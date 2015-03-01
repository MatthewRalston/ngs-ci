require 'bio'
require 'ruby-prof'
require 'parallel'
require 'bio-samtools'

module LCI

  # A calculator calculates the library complexity index.
  #
  # @!attribute [r] lci
  class Calculator
    attr_reader :lci, :block_size, :buffer, :chroms

    # A new calculator to compute the library complexity index given
    # a loaded Bio::DB::Sam object and optional thread argument.
    #
    # @param bam [Bio::DB::Sam] Opened bam file with loaded reference.
    # @param threads [Int] The number of threads used to compute LCI.
    # @param strand [String] One of [FR RF F] or nil for strandedness.
    def initialize(bam, reference, strand: nil, threads: 1)
      @block_size = 1000
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
          raise LCIError.new "Strand specific option #{opts.strand} is invalid." +
      " It must be one of: [FR, RF, F]"
        end
        @strand = strand.downcase
      else
        @strand = nil
      end
    end

    # Calculation of the library complexity index
    #
    def run(runtime: false)
      RubyProf.start if runtime
      # Convert each aligned read to Read clas
      chroms={}
      @chroms.each do |chrom,size|
        chroms[chrom]={}
        disk_accesses = (size/@block_size).ceil
        #disk_accesses.times do |i|
        Parallel.each((0...disk_accesses).to_a,:in_threads => @threads) do |i|
          readblock(chrom,i).each do |key,val|
            if chroms[chrom][key]
              chroms[chrom][key] += val
            else
              chroms[chrom][key] = val
            end
          end
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

    # Reads a single block from the disk and calculates the LCI
    #
    # @param chrom [String] The chromosome from the bam file
    # @param i [Integer] The number of blocks that have been read
    # @return localLCI [Hash<Symbol,Array>]
    #   * :+ (Array[Integer]) The LCI for the + strand
    #   * :- (Array[Integer]) The LCI for the - strand
    def readblock(chrom,i)
      reads=[]
      results = @strand ? {"+" => [],"-" => []}: {nil => []}
      start = [0,(i * @block_size) - @buffer].max
      stop = (i + 1) * @block_size
      @bam.fetch(chrom,start,stop) {|read| reads << convert(read)}
      start += @buffer unless start == 0
      reads.compact!
      reads.sort! {|r| r.start} unless reads.empty?
      (start...stop).to_a.each do |b|
      #Parallel.each((start...stop).to_a, :in_threads => @threads) do |b|
      #Parallel.each((start...stop).to_a, :in_processes => @threads) do |b|
        aligned = reads.select{|r| r.start <= b && r.stop >= b}.
          group_by{|r|r.strand} # Alternative selection of reads
        results.keys.each do|key|
          val = aligned[key] ? aligned[key] : []
          results[key] << [b,lci(val)]
        end
        
      end
      results
    end

    # Calculates library complexity index for a single base
    # 
    # @param reads [Array<LCI::Read>] A group of reads aligned to a single base.
    # @return lci [Float] 
    def lci(reads)
      #x=(unique_start_sites(reads)*average_overlap(reads)/@buffer.to_f).round(4)
      x=(unique_start_sites(reads)*average_overlap(unique_reads(reads))/@buffer.to_f).round(4)
      raise LCIError.new "lci is < 0: #{x}" unless x >= 0
      x
    end

    # Calculates the unique reads from a group of reads
    #
    # @param reads [Array<LCI::Read>] A group of reads aligned to a single base.
    # @return unique [Array<LCI::Read>] A group of reads with unique start sites.
    def unique_reads(reads)
      unique=[]; ustarts=[]
      reads.each do |read|
        unless ustarts.include? read.start
          unique << read
          ustarts << read.start
        end
      end
      unique
    end

    # Calculates the number of unique start sites
    #
    # @param reads [Array<LCI::Read>] Array of reads
    # @return [Integer] Number of unique start sites
    def unique_start_sites(reads)
      reads.map(&:start).uniq.size
    end

    # Calculates average overlap between a group of reads
    #
    # @param reads [Array<LCI::Read>] Array of reads
    # @return avg_overlap [Integer] Average overlap between reads
    def average_overlap(reads)
      avg=0
      reads.each do |r1|
        sum=0
        reads.each do |r2|
          unless r1 == r2
            sum+=overlap(r1,r2).to_f
          end
        end
        avg+=(sum/(reads.size-1))/reads.size unless reads.size == 1
      end
      return avg
    end

    # Calculation of the overlap between two reads
    # 
    # @param read1 [LCI::Read] First read to be compared
    # @param read2 [LCI::Read] First read to be compared
    # @return overlap_length [Integer] Length of overlap
    def overlap(read1,read2)
      read1.start > read2.start ? (read2.stop - read1.start) : (read1.stop - read2.start)
    end

    # Loads the read length from a bam file into the @buffer variable
    #
    def read_length
      buffer=0
      stats=@bam.index_stats.select {|k,v| k != "*" && v[:mapped_reads] > 0}
      if stats.empty?
        raise LCIIOError.new "BAM file is empty! Check samtools idxstats."
      else
        @bam.view do |read|
          @buffer=read.seq.size
          break
        end
      end
    end

    # Converts strand specific BAM read into a sequence object format
    # Uses the @strand instance variable to determine the strand of conversion
    # 
    # @param read [Bio::DB::Alignment] Read to be converted.
    # @return read [LCI::Read] Converted Read object
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
    # @return read [LCI::Read] Converted Read object
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
    # @return read [LCI::Read] Converted Read object
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
    # @return read [LCI::Read] Converted Read object
    def f(read)
      read.query_strand ? newread(read,strand:"+") : newread(read,strand:"-")
    end

    # Creates a new read with optional strand argument
    #
    # @param read [Bio::DB::Alignment] Aligned read to be converted
    # @param strand [String] Strand of read
    # @return read [LCI::Read] Converted Read object
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
          file.puts("Chrom,Base,Strand,LCI")
          @results.each do |chrom,results|
            results.each do |strand,val|
              val.each do |x|
                file.puts([chrom,x[0],strand,x[1]].join(","))
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
