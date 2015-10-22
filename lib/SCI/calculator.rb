require 'bio'
require 'ruby-prof'
require 'parallel'
require 'bio-samtools'

module SCI

  # A calculator calculates the sequencing complexity index.
  #
  # @!attribute [r] sci
  class Calculator
    attr_reader :sci, :block_size, :buffer, :chroms

    # A new calculator to compute the sequencing complexity index given
    # a loaded Bio::DB::Sam object and optional thread argument.
    #
    # @param bam [Bio::DB::Sam] Opened bam file with loaded reference.
    # @param threads [Int] The number of threads used to compute SCI.
    # @param strand [String] One of [FR RF F] or nil for strandedness.
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
      read_length
      if strand
        unless %w(FR RF F).include?(strand)
          raise SCI::SCIError.new "Strand specific option #{opts.strand} is invalid." +
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

    # Reads a single block from the disk and calculates the SCI
    #
    # @param chrom [String] The chromosome from the bam file
    # @param i [Integer] The number of blocks that have been read
    # @return localSCI [Hash<Symbol,Array>]
    #   * :+ (Array[Integer]) The SCI for the + strand
    #   * :- (Array[Integer]) The SCI for the - strand
    def readblock(chrom,i)
      reads=[]
      results = @strand ? {"+" => [],"-" => []}: {nil => []}
      start = [0,(i * @block_size) - @buffer].max
      stop = [(i + 1) * @block_size, self.chroms[chrom]].min
      @bam.fetch(chrom,start,stop) {|read| reads << convert(read)}
      start += @buffer unless start == 0
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
    # @param reads [Array<SCI::Read>] A group of reads aligned to a single base.
    # @return sci [Float] 
    def sci(reads)
      numreads=reads.size
      # Groups reads by start site
      # selects the largest read length from the groups
      reads = reads.group_by(&:start).map{|k,v| v.max{|x,y| (x.stop-x.start).abs <=> (y.stop-y.start).abs}}
      o = summed_overlaps(reads)
      uniquereads = reads.size
      return [numreads,uniquereads,(@buffer*o.to_f/@denom).round(4),(300*uniquereads*o/(2*@denom)).round(4)]
    end

    # Calculates summed overlap between a group of reads
    #
    # @param reads [Array<SCI::Read>] Array of reads
    # @return avg_overlap [Integer] Summed overlap between reads
    def summed_overlaps(reads)
      numreads = reads.size
      sum=0
      unless numreads == 1
        i = 0
        while i < numreads
          r1 = reads[i] # for each of n reads
          sum+=reads.
                reject{|r| r == r1}. # select the n-1 other reads
                map{|r| overlap(r,r1)}. # calculate their overlap to r1
                reduce(:+)
          i+=1
        end
      end
      return sum
    end    

    # Calculation of the overlap between two reads
    # 
    # @param read1 [SCI::Read] First read to be compared
    # @param read2 [SCI::Read] First read to be compared
    # @return overlap_length [Integer] Length of overlap
    def overlap(read1,read2)
      if read1.start > read2.start
        if read1.stop < read2.stop # Read 1 is inside read 2
          read1.stop - read1.start
        else # Normal overlap
          read2.stop - read1.start
        end
      else
        if read1.stop > read2.stop # Read 2 is inside read 1
          read2.stop - read2.start
        else # Normal overlap
          read1.stop - read2.start
        end
      end
    end

    # Loads the read length from a bam file into the @buffer variable
    #
    def read_length
      buffer=0
      stats=@bam.index_stats.select {|k,v| k != "*" && v[:mapped_reads] > 0}
      if stats.empty?
        raise SCIIOError.new "BAM file is empty! Check samtools idxstats."
      else
        i=0
        lengths=[]
        test = @block_size
        while i <= test 
          @bam.view do |read|
            lengths << read.seq.size
            i +=1
          end
          if i == test && lengths.size < 100
            test += @block_size
          end
        end
        @buffer = lengths.max
        @denom = @buffer**2 * (@buffer - 1)**2
      end
    end

    # Converts strand specific BAM read into a sequence object format
    # Uses the @strand instance variable to determine the strand of conversion
    # 
    # @param read [Bio::DB::Alignment] Read to be converted.
    # @return read [SCI::Read] Converted Read object
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
    # @return read [SCI::Read] Converted Read object
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
    # @return read [SCI::Read] Converted Read object
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
    # @return read [SCI::Read] Converted Read object
    def f(read)
      read.query_strand ? newread(read,strand:"+") : newread(read,strand:"-")
    end

    # Creates a new read with optional strand argument
    #
    # @param read [Bio::DB::Alignment] Aligned read to be converted
    # @param strand [String] Strand of read
    # @return read [SCI::Read] Converted Read object
    def newread(read,strand: nil)
      Read.new(read.pos,read.pos+read.seq.size,strand: strand)
    end

    # Acquires names and sizes of reference sequences included in the bam file
    #
    # @param reference [String] Path to reference fasta file.
    # @return chromosomes [Hash<Symbol,Object>] A dictionary of chromosome sizes
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
          file.puts("Chrom,Base,Strand,Depth,Unique_Reads,Overlap,SCI")
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
