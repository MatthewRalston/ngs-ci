require 'spec_helper'

testbam="test_files/test.bam"
testfasta="test_files/NC_001988.ffn"

before(:all) do
  @bam=Bio::DB::Sam.new(:bam=>testbam,:fasta=>testfasta)
  @bam.open
  @reads=[]
  @bam.fetch("NC_001988.2",0,200) {|x| @reads << x}
  @calc=LCI::Calculator.new(testbam,testfasta)
end

describe '#reference_sequences' do
  it "retrieves reference sequences" do
    expect(@calc.reference_sequences(testfasta).keys).to include "NC_001988.2"
  end
end

describe '#newread' do
  it "converts an alignment object to a read object" do
    expect(newread(@reads[0])).to be_instance_of LCI::Read
  end  
end


describe '#convert' do


end



