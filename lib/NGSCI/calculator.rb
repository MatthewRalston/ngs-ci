require 'bio'
require 'parallel'
require 'bio-samtools'
require 'ruby-prof'

module NGSCI

  # A calculator calculates the sequencing complexity index.
  #
  # @authoer Matthew Ralston
  # @abstract A class for calculating the complexity index on next generation sequencing reads
  # @attr_reader [Integer] block_size The block size for parallelizing disk access
  # @attr_reader [Hash<Symbol,Integer>] chroms A hash of chromosomes and their sizes
  # @attr_reader [Integer] read_length The read length obtained from a bam file
  # @attr_reader [Integer] denominator The denominator and normalization factors calculated from the read length
  class Calculator
    attr_reader :block_size, :chroms, :read_length, :denominator

    # A new calculator to compute the sequencing complexity index given
    # a loaded Bio::DB::Sam object and optional thread argument.
    #
    # @param [Bio::DB::Sam] bam Opened bam file with loaded reference.
    # @param [Int] threads The number of threads used to compute NGSCI.
    # @param [String] strand One of [FR RF F] or nil for strandedness.
    def initialize(bam, reference, strand: nil, threads: 1)
      @block_size = 1600
      @results = nil
      @reference=reference
      @bam = Bio::DB::Sam.new(:bam=>bam,:fasta=>reference)
      unless @bam.indexed?
        @bam.index
      end
      @bam.open
      @threads = threads
      @chroms = reference_sequences(reference)
      @read_length = NGSCI::Calculator.read_length_calc(@bam,@block_size)
      @denominator = denominator_calc(@read_length)
      if strand
        unless %w(FR RF F).include?(strand)
          raise NGSCI::NGSCIError.new "Strand specific option #{opts.strand} is invalid." +
      " It must be one of: [FR, RF, F]"
        end
        @strand = strand.downcase
      else
        @strand = nil
      end
    end

    # Calculation of the sequencing complexity index
    #
    def run(runtime: false)
      RubyProf.start if runtime
      # Convert each aligned read to Read clas
      chroms={}
      @chroms.each do |chrom,size|
        chroms[chrom] = @strand ? {"+"=>[],"-"=>[]} : {nil=>[]}
        disk_accesses = (size/@block_size.to_f).ceil
=begin
        # N O N - P A R A L L E L
        i=0
        while i < disk_accesses

          readblock(chrom,i).each do |key,val|
            chroms[chrom][key] += val
          end
          i+=1
        end
=end

        data = Parallel.map((0...disk_accesses).to_a,:in_processes => @threads) do |i|
          readblock(chrom,i)
        end
        chroms[chrom].keys.each do |key|
          chroms[chrom][key] = data.map{|x| x[key]}.flatten(1)
        end

      end
      
      # Printing runtime information for optimization
      if runtime
        runtime=RubyProf.stop
        printer=RubyProf::FlatPrinter.new(runtime)
        printer.print(STDOUT)
      end
      @results = chroms
    end

    # Reads a single block from the disk and calculates the NGSCI
    #
    # @param chrom [String] The chromosome from the bam file
    # @param i [Integer] The number of blocks that have been read
    # @return localNGSCI [Hash<Symbol,Array>]
    #   * :+ (Array[Integer]) The NGSCI for the + strand
    #   * :- (Array[Integer]) The NGSCI for the - strand
    def readblock(chrom,i)
      reads=[]
      results = @strand ? {"+" => [],"-" => []}: {nil => []}
      start = [0,(i * @block_size) - @read_length].max
      stop = [(i + 1) * @block_size, self.chroms[chrom]].min
      @bam.fetch(chrom,start,stop) {|read| reads << convert(read)}
      start += @read_length unless start == 0
      reads.compact!
      reads.sort_by!(&:start) unless reads.empty?
      x=0
      bases = (start...stop).to_a
      block = stop - start
      while x < block
        b = bases[x]
        aligned = reads.select{|r| r.start <= b && r.stop - 1 >= b}.group_by &:strand
        results.keys.each do|key|
          results[key] << [b,*sci(aligned[key] || [])]
        end
        x+=1
      end
      return results
    end

    # Calculates sequencing complexity index for a single base
    # 
    # @param reads [Array<NGSCI::Read>] A group of reads aligned to a single base.
    # @return localNGSCI [Array<Integer,Integer,Float,Float>]
    def sci(reads)
      numreads=reads.size
      # Groups reads by start site
      # selects the largest read length from the groups
      reads = reads.group_by(&:start).map{|k,v| v.max{|x,y| x.length <=> y.length}}
      d = summed_dissimilarity(reads)
      uniquereads = reads.size
      return [numreads,uniquereads,(d.to_f/@read_length).round(4),(uniquereads*d/@denominator).round(4)]
    end


    # Calculates summed dissimilarity between a group of reads
    #
    # @param reads [Array<NGSCI::Read>] Array of reads
    # @return sum_dissimilarity [Integer] Summed dissimilarity between reads
    def summed_dissimilarity(reads)
      numreads = reads.size
      sum=0
      unless numreads <= 1
        i = 0
        while i < numreads
          r1 = reads[i] # for each of n reads
          sum+=reads.
                reject{|r| r == r1}. # select the n-1 other reads
                map{|r| dissimilarity(r,r1)}. # calculate their overlap to r1
                reduce(:+)
          i+=1
        end
      end
      return sum
    end

    # Calculation of the dissimilarity between two reads
    #
    # @param read1 [NGSCI::Read]  First read to be compared
    # @param read2 [NGSCI::Read]  Second read to be compared
    # @return unique_bases [Integer] Length of non-overlapping bases
    def dissimilarity(read1,read2)
      if read1.start > read2.start
        if read1.stop < read2.stop # Read 1 is inside read 2
          (read1.start - read2.start) + (read2.stop - read1.stop)
        else # Normal overlap
          read1.start - read2.start
        end
      else
        if read1.stop > read2.stop # Read 2 is inside read 1
          (read2.start - read1.start) + (read1.stop - read2.stop)
        else # Normal overlap
          read2.start - read1.start
        end
      end
    end

    # Calculates the read length of a bam file by sampling at least on full block of reads
    #
    # @param [Bio::DB::Sam] bam A bam reader object
    # @param [Integer] block_size The number of reads to read from a bam file
    # @return [Integer] read_length The read length acquired from reading a block at a time until at least 100 reads are acquired
    def self.read_length_calc(bam,block_size)
      stats=bam.index_stats.select {|k,v| k != "*" && v[:mapped_reads] > 0}
      if stats.empty?
        raise NGSCIIOError.new "BAM file is empty! Check samtools idxstats."
      end
      i=0
      lengths=[]
      test = block_size
      while i <= test 
        bam.view do |read|
          lengths << read.seq.size
          i +=1
        end
        if i == test && lengths.size < 100
          test += block_size
        end
      end
      lengths.max
    end

    # Calculates the denominator for the complexity index from the read length, assuming maximum saturation (i.e. number of unique reads == read_length)
    #
    # @param [Integer] read_length The read length
    # @return [Float] denominator The denominator including normalization factors for the complexity index
    def denominator_calc(read_length)
      read_length*3*max_summed_dissimilarity(read_length)/(read_length - 1).to_f
    end

    # Calculates the average summed dissimilarity (per read) of that read to all other reads
    #
    # @param [Integer] read_length The read length
    # @return [Integer] avg_summed_dissimilarity
    def max_summed_dissimilarity(read_length)
      # For each unique read under maximum saturation, calculate the sum of dissimilarities for that read to all other reads
      summed_dissimilarities = (1..read_length).to_a.map { |r| 
        (read_length ** 2) / 2 - read_length*r + read_length/2 + r**2 - r }.reduce(:+)
    end

    # Converts strand specific BAM read into a sequence object format
    # Uses the @strand instance variable to determine the strand of conversion
    # 
    # @param [Bio::DB::Alignment] read Read to be converted.
    # @return [NGSCI::Read] read Converted Read object
    def convert(read)
      unless read.query_unmapped
        if @strand
          return self.send(@strand.to_sym,read)
        else
          return newread(read)
        end
      end
      return nil
    end


    # Converts strand specific BAM read into a sequence object format
    # Assumes paired-end strand-specific sequencing with "fr" chemistry
    # 
    # @param read [Bio::DB::Alignment] Read to be converted.
    # @return read [NGSCI::Read] Converted Read object
    def fr(read)
      if read.first_in_pair
        read.query_strand ? newread(read,strand:"+") : newread(read,strand:"-")
      else
        read.query_strand ? newread(read,strand:"-") : newread(read,strand:"+")
      end
    end


    # Converts strand specific BAM read into a sequence object format
    # Assumes paired-end strand-specific sequencing with "rf" chemistry
    # 
    # @param read [Bio::DB::Alignment] Read to be converted.
    # @return read [NGSCI::Read] Converted Read object
    def rf(read)
      if read.first_in_pair
        read.query_strand ? newread(read,strand:"-") : newread(read,strand:"+")
      else
        read.query_strand ? newread(read,strand:"+") : newread(read,strand:"-")
      end      
    end


    # Converts strand specific BAM read into a sequence object format
    # Assumes single-end strand-specific sequencing with "f" chemistry
    # 
    # @param read [Bio::DB::Alignment] Read to be converted.
    # @return read [NGSCI::Read] Converted Read object
    def f(read)
      read.query_strand ? newread(read,strand:"+") : newread(read,strand:"-")
    end

    # Creates a new read with optional strand argument
    #
    # @param read [Bio::DB::Alignment] Aligned read to be converted
    # @param strand [String] Strand of read
    # @return read [NGSCI::Read] Converted Read object
    def newread(read,strand: nil)
      Read.new(read.pos,read.pos+read.seq.size,strand: strand)
    end

    # Acquires names and sizes of reference sequences included in the bam file
    #
    # @param reference [String] Path to reference fasta file.
    # @return chromosomes [Hash<Symbol,Integer>] A dictionary of chromosome sizes
    def reference_sequences(reference)
      chromosomes={}
      Bio::FastaFormat.open(@reference).each_entry do |f| 
        chromosomes[f.entry_id]=f.seq.size
      end      
      chromosomes.select {|chrom| @bam.index_stats.keys.include?(chrom)}
    end

    # Exports the results to outfile
    #
    # @param outfile [String] Path to outfile
    def export(outfile)
      if @results
        File.open(outfile,'w') do |file|
          file.puts("Chrom,Base,Strand,Depth,Unique_Reads,Overlap,NGS-CI")
          @results.each do |chrom,results|
            results.each do |strand,val|
              val.each do |x|
                file.puts([chrom,x[0],strand,*x[1..-1]].join(","))
              end
            end
          end
        end
        return outfile
      else
        return nil
      end
    end
  end # End calculator class
end
