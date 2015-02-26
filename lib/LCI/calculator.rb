require 'bio'
require 'ruby-prof'

module LCI

  # A calculator calculates the library complexity index.
  #
  # @!attribute [r] lci
  class Calculator
    attr_reader :lci, :block_size

    # A new calculator to compute the library complexity index given
    # a loaded Bio::DB::Sam object and optional thread argument.
    #
    # @param bam [Bio::DB::Sam] Opened bam file with loaded reference.
    # @param threads [Int] The number of threads used to compute LCI.
    # @param strand [String] One of [FR RF F] or nil for strandedness.
    def initialize(bam, reference, strand: nil, threads: 1)
      @block_size = 1000
      
      @reference=reference
      @bam = Bio::DB::Sam.new(:bam=>bam,:fasta=>reference)
      @bam.open
      @threads = threads
      @chroms = reference_sequences(reference)
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
    def librarycomplexityindex
      RubyProf.start
      # Convert each aligned read to Read class
      @chroms.each do |chrom,size|
        disk_accesses = size/Calculator.block_size.ceil
        disk_accesses.times {|x| iteration}
      end
      
      # Printing runtime information for optimization
      runtime=RubyProf.stop
      printer=RubyProf::FlatPrinter.new(runtime)
      printer.print(STDOUT)
    end



    # Converts strand specific BAM read into a sequence object format
    # 
    # @param read [Bio::DB::Alignment] Read to be converted.
    def convert(read)
      if @strand
        self.send(@strand.to_sym,read)
      else
        newread(read)
      end
    end


    # Converts strand specific BAM read into a sequence object format
    # 
    # @param read [Bio::DB::Alignment] Read to be converted.
    def fr(read)
      if read.first_in_pair
        read.query_strand ? newread(read,strand:"+") : newread(read,strand:"-")
      else
        read.query_strand ? newread(read,strand:"-") : newread(read,strand:"+")
      end
    end


    # Converts strand specific BAM read into a sequence object format
    # 
    # @param read [Bio::DB::Alignment] Read to be converted.
    def rf(read)
      if read.first_in_pair
        read.query_strand ? newread(read,strand:"-") : newread(read,strand:"+")
      else
        read.query_strand ? newread(read,strand:"+") : newread(read,strand:"-")
      end      
    end


    # Converts strand specific BAM read into a sequence object format
    # 
    # @param read [Bio::DB::Alignment] Read to be converted.
    def f(read)
      read.query_strand ? newread(read,strand:"+") : newread(read,strand:"+")
    end

    # Creates a new read with optional strand argument
    #
    # @param
    # @param strand [String] Strand of read
    def newread(read,strand: nil)
      Read.new(read.pos,read.pos+read.qlen,strand: strand)
    end

    # Acquires names and sizes of reference sequences
    #
    # @param reference [String] Path to reference fasta file.
    def reference_sequences(reference)
      chromosomes={}
      Bio::FastaFormat.open(@reference).each_entry do |f| 
        chromosomes[f.entry_id]=f.seq.size
      end
      chromosomes
    end
    
  end # End calculator class
end
