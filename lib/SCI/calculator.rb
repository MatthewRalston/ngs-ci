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
          raise SCIError.new "Strand specific option #{opts.strand} is invalid." +
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
        reads = []
        @bam.fetch(chrom,0,size) {|read| reads << convert(read)}
        chroms[chrom] = @strand ? {"+"=>[],"-"=>[]} : {nil=>[]}
        process(reads,size).each do |hash|
          hash.each do |key,val|
            chroms[chrom][key] << val
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

    # Calculates the SCI in parallel for each base in a chromosome
    #
    # @param reads [Array<SCI::Read>] A list of Read objects for processing
    # @param size [Integer] The size of the chromosome
    # @return [Array<Hash>]
    #   * :+ (Array[Integer]) The position and SCI for the + strand
    #   * :- (Array[Integer]) The position and SCI for the - strand
    def process(reads,size)
      strands = @strand ? {"+" => [],"-" => []}: {nil => []}
      reads.compact!
      reads.sort_by!(&:start) unless reads.empty?
      results = Parallel.map((0...size).to_a, :in_processes => @threads) do |b|
        aligned = reads.select{|r| r.start <= b && r.stop >= b}.uniq(&:start).group_by &:strand
        temp = {}
        strands.keys.each do |key|
          temp[key] = [b,sci(aligned[key] || [])]
        end
        temp
      end
      return results
    end

    # Calculates sequencing complexity index for a single base
    # 
    # @param reads [Array<SCI::Read>] A group of reads aligned to a single base.
    # @return sci [Float] 
    def sci(reads)
      #x=(unique_start_sites(reads)*average_overlap(reads)/@buffer.to_f).round(4)
      #raise SCIError.new "sci is < 0: #{x}" unless x >= 0
      return (reads.map(&:start).size*average_overlap(reads)/@buffer.to_f).round(4)
    end

    # Calculates average overlap between a group of reads
    #
    # @param reads [Array<SCI::Read>] Array of reads
    # @return avg_overlap [Integer] Average overlap between reads
    def average_overlap(reads)
      avg=0
      for r1 in reads
        avg+=(reads.
              reject{|r| r == r1}.map{|r| overlap(r,r1)}.
              reduce(:+)/(reads.size - 1)/reads.size.to_f) unless reads.size == 1
      end
      return avg
    end

    # Calculation of the overlap between two reads
    # 
    # @param read1 [SCI::Read] First read to be compared
    # @param read2 [SCI::Read] First read to be compared
    # @return overlap_length [Integer] Length of overlap
    def overlap(read1,read2)
      return read1.start > read2.start ? (read2.stop - read1.start) : (read1.stop - read2.start)
    end

    # Loads the read length from a bam file into the @buffer variable
    #
    def read_length
      buffer=0
      stats=@bam.index_stats.select {|k,v| k != "*" && v[:mapped_reads] > 0}
      if stats.empty?
        raise SCIIOError.new "BAM file is empty! Check samtools idxstats."
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
          file.puts("Chrom,Base,Strand,SCI")
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
