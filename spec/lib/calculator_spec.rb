require 'spec_helper'
require 'bio-samtools'

testbam="spec/test_files/test.bam"
emptybam="spec/test_files/empty.bam"
testfasta="spec/test_files/test.fa"
testout="spec/test_files/testfile.txt"



describe "#run" do
  context "during a strand specific run" do
    before(:each) do
      @calc=NGSCI::Calculator.new(testbam,testfasta,strand:"FR")
      @testchrom=@calc.chroms.keys[0]
    end
    it "returns a hash" do
      expect(@calc.run(runtime: false)).to be_instance_of Hash
    end
    it "returns the hash with keys of the chromosomes names" do
      expect(@calc.run.keys).to eq(@calc.chroms.keys)
    end
    it "returns the hash with keys for each strand" do
      expect(@calc.run[@testchrom].keys).to eq(%w(+ -))
    end
    it "returns NGSCI for each base of the genome" do
      expect(@calc.run[@testchrom]["+"].size).to eq(@calc.chroms[@testchrom])
    end
  end

  context "during an unstranded run" do
    before(:each) do
      @calc=NGSCI::Calculator.new(testbam,testfasta)
      @testchrom=@calc.chroms.keys[0]
    end
    it "returns the hash with keys of the chromosomes names" do
      expect(@calc.run.keys).to eq(@calc.chroms.keys)
    end
    it "returns the hash with a nil strand key" do
      expect(@calc.run[@testchrom].keys[0]).to be nil
    end
    it "returns NGSCI for each base of the genome" do
      expect(@calc.run[@testchrom][nil].size).to eq(@calc.chroms[@testchrom])
    end
  end
end

describe "#readblock" do
  context "when reading the first block" do
    it "returns a hash with an array of length @block_size" do
      @calc=NGSCI::Calculator.new(testbam,testfasta)
      results=@calc.readblock(@calc.chroms.keys[0],0)
      result_length=results[results.keys[0]].size
      expect(result_length).to eq(@calc.block_size)
    end
  end
  context "when reading any other block" do
    it "returns a hash with and array of length @block_size" do
      @calc=NGSCI::Calculator.new(testbam,testfasta)
      results=@calc.readblock(@calc.chroms.keys[0],1)
      result_length=results[results.keys[0]].size
      expect(result_length).to eq(@calc.block_size)      
    end
  end
end

describe "#sci" do
  context "when passed an array of read objects" do
    before(:each) do
      @calc = NGSCI::Calculator.new(testbam,testfasta)
      @bam = Bio::DB::Sam.new(:bam=>testbam,:fasta=>testfasta)
      @bam.open
      @reads = []
      @bam.fetch("NC_001988.2",75,75){|x| read = @calc.convert(x); @reads << read unless read.nil?}
      @reads = @reads.uniq{|r|r.start}
    end
    it "returns an array" do
      expect(@calc.sci(@reads)).to be_kind_of(Array)
    end
    it "returns the sequencing complexity index" do
      expect(@calc.sci(@reads)[-1]).to eq(0.0)
    end
  end
  context "when passed an empty array" do
    it "returns nil" do
      @calc = NGSCI::Calculator.new(testbam,testfasta)
      empty_sci = @calc.sci([])[-1]
      expect(empty_sci).to be_zero
    end
  end
end

describe "#read_length" do
  it "calculates the read length" do
    @calc=NGSCI::Calculator.new(testbam,testfasta)
    expect(@calc.buffer).to eq(76)
  end
  it "fails on an empty bam file" do
    expect{NGSCI::Calculator.new(emptybam,testfasta)}.to raise_error(NGSCI::NGSCIIOError)
    `rm #{emptybam}.bai`
  end
end

describe "#summed_overlaps" do
  it "returns an int" do
    @bam=Bio::DB::Sam.new(:bam=>testbam,:fasta=>testfasta)
    @bam.open
    @reads = []
    @calc=NGSCI::Calculator.new(testbam,testfasta)
    @bam.fetch("NC_001988.2",8,75) {|x| read=@calc.convert(x); @reads << read if read}
    @reads = @reads.uniq{|r| r.start}
    expect(@calc.summed_overlaps(@reads)).to be_an(Integer)
  end
  context "when passed an array of read objects" do
    before(:each) do
      @bam=Bio::DB::Sam.new(:bam=>testbam,:fasta=>testfasta)
      @bam.open
      @reads = []
      @calc=NGSCI::Calculator.new(testbam,testfasta)
      @bam.fetch("NC_001988.2",8,75) {|x| read=@calc.convert(x); @reads << read if read}
      @reads = @reads.uniq{|r| r.start}
    end
    it "returns the #overlap of two reads" do
      summed_overlap = 2*@calc.overlap(@reads[0],@reads[1])
      expect(@calc.summed_overlaps(@reads[0..1])).to eq(summed_overlap)
    end

    it "calculates the average overlap between a group of reads" do
      expect(@calc.summed_overlaps(@reads[0..7]).round(4)).to eq(380.0)
    end
  end
  context "when passed an array with a single read object" do
    it "returns zero" do
      @bam=Bio::DB::Sam.new(:bam=>testbam,:fasta=>testfasta)
      @bam.open
      @reads=[]
      @calc=NGSCI::Calculator.new(testbam,testfasta)
      @bam.fetch("NC_001988.2",8,75) {|x| read=@calc.convert(x); @reads << read if read}
      expect(@calc.summed_overlaps([@reads[0]])).to be_zero
    end
  end
  context "when passed an empty array" do
    it "returns zero" do
      @calc=NGSCI::Calculator.new(testbam,testfasta)
      expect(@calc.summed_overlaps([])).to be_zero
    end
  end
end

describe "#overlap" do
  before(:each) do
    @bam=Bio::DB::Sam.new(:bam=>testbam,:fasta=>testfasta)
    @bam.open
    @reads=[]
    @bam.fetch("NC_001988.2",0,200) {|x| @reads << x}
    @calc=NGSCI::Calculator.new(testbam,testfasta)
    @read1=@calc.convert(@reads[2])
    @read2=@calc.convert(@reads[3])
  end
  it "calculates the overlap between two reads" do
    expect(@calc.overlap(@read1,@read2)).to eq(14)
  end

  it "calculates the overlap regardless of order" do
    expect(@calc.overlap(@read2,@read1)).to eq(14)
  end
end

describe '#reference_sequences' do
  before(:each) do
    @calc=NGSCI::Calculator.new(testbam,testfasta)
  end

  it "retrieves reference sequences" do
    expect(@calc.reference_sequences(testfasta).keys).to include "NC_001988.2"
  end
end

describe '#newread' do
  before(:each) do
    @bam=Bio::DB::Sam.new(:bam=>testbam,:fasta=>testfasta)
    @bam.open
    @reads=[]
    @bam.fetch("NC_001988.2",0,200) {|x| @reads << x}
    @calc=NGSCI::Calculator.new(testbam,testfasta)
  end

  it "converts an alignment object to a read object" do
    expect(@calc.newread(@reads[0])).to be_instance_of NGSCI::Read
  end  
end

describe "#fr" do
  before(:each) do
    @bam=Bio::DB::Sam.new(:bam=>testbam,:fasta=>testfasta)
    @bam.open
    @reads=[]
    @bam.fetch("NC_001988.2",0,200) {|x| @reads << x}
    @calc=NGSCI::Calculator.new(testbam,testfasta)
    @testpair=[@reads[5],@reads[10]]
  end

  it "converts the first read with FR chemistry" do
    first=@calc.fr(@testpair[0])
    expect(first.strand).to eq("+")
  end
  it "converts the second read with FR chemistry" do
    second=@calc.fr(@testpair[1])
    expect(second.strand).to eq("+")
  end
end

describe "#rf" do
  before(:each) do
    @bam=Bio::DB::Sam.new(:bam=>testbam,:fasta=>testfasta)
    @bam.open
    @reads=[]
    @bam.fetch("NC_001988.2",0,200) {|x| @reads << x}
    @calc=NGSCI::Calculator.new(testbam,testfasta)
    @testpair=[@reads[5],@reads[10]]
  end
  it "converts the first read with RF chemistry" do
    first=@calc.rf(@testpair[0])
    expect(first.strand).to eq("-")
  end

  it "converts the second read with RF chemistry" do
    second=@calc.rf(@testpair[1])
    expect(second.strand).to eq("-")
  end
end

describe "#f" do
  before(:each) do
    @bam=Bio::DB::Sam.new(:bam=>testbam,:fasta=>testfasta)
    @bam.open
    @reads=[]
    @bam.fetch("NC_001988.2",0,200) {|x| @reads << x}
    @calc=NGSCI::Calculator.new(testbam,testfasta)
    @testpair=[@reads[5],@reads[10]]
  end
  it "converts a read with F chemistry on the + strand" do
    first=@calc.f(@testpair[0])
    expect(first.strand).to eq("-")
  end
  it "converts a read with F chemistry on the - strand" do
    second=@calc.f(@testpair[1])
    expect(second.strand).to eq("+")
  end

end

describe '#convert' do
  before(:each) do
    @bam=Bio::DB::Sam.new(:bam=>testbam,:fasta=>testfasta)
    @bam.open
    @reads=[]
    @bam.fetch("NC_001988.2",0,200) {|x| @reads << x}
  end
  it "converts an alignment object to a read object" do
    calc=NGSCI::Calculator.new(testbam,testfasta)
    expect(calc.convert(@reads[2])).to be_instance_of NGSCI::Read
  end

  it "returns nil for an unmapped read" do
    calc=NGSCI::Calculator.new(testbam,testfasta)
    expect(calc.convert(@reads[1])).to be_nil
  end

  it "converts the first read in FR chemistry aligned to the + strand" do
    calc=NGSCI::Calculator.new(testbam,testfasta,strand:"FR")
    testpair=[@reads[5],@reads[10]]
    first=calc.convert(testpair[1])
    expect(first.strand).to eq("+")
  end
  it "converts the second read in FR chemistry aligned to the - strand" do
    calc=NGSCI::Calculator.new(testbam,testfasta,strand:"FR")
    testpair=[@reads[5],@reads[10]]
    second=calc.convert(testpair[0])
    expect(second.strand).to eq("+")
  end
  it "converts the first read in RF chemistry aligned to the + strand" do
    calc=NGSCI::Calculator.new(testbam,testfasta,strand:"RF")
    testpair=[@reads[5],@reads[10]]
    first=calc.convert(testpair[1])
    expect(first.strand).to eq("-")
  end
  it "converts the second read in FR chemistry aligned to the - strand" do
    calc=NGSCI::Calculator.new(testbam,testfasta,strand:"RF")
    testpair=[@reads[5],@reads[10]]
    second=calc.convert(testpair[0])
    expect(second.strand).to eq("-")
  end

end

describe "#export" do
  context "calculator has not run" do
    it "returns nil" do
      calc=NGSCI::Calculator.new(testbam,testfasta)
      expect(calc.export(testout)).to be nil
    end
  end
  context "calculator has run" do
    after(:all) do
      `rm #{testout}`
    end
    it "returns outfile" do
      calc=NGSCI::Calculator.new(testbam,testfasta)
      calc.run
      expect(calc.export(testout)).to eq(testout)
    end
  end
end


