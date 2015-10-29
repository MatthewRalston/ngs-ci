require 'spec_helper'


describe NGSCI::Read do
  context "before created" do
    it "fails to instantiate on a string start site" do
      expect{NGSCI::Read.new("foo",3)}.to raise_error(NGSCI::NGSCIError)
    end

    it "fails to instantiate on a string stop site" do
      expect{NGSCI::Read.new(1,"foo")}.to raise_error(NGSCI::NGSCIError)
    end

    it "fails to instantiate when the stop site is greater than the start site" do
      expect{NGSCI::Read.new(3,1)}.to raise_error(NGSCI::NGSCIError)
    end

    it "fails to instantiate on an improper strand argument" do
      expect{NGSCI::Read.new(1,3,strand:"foo")}.to raise_error(NGSCI::NGSCIError)
    end

    it "fails to instantiate without the three necessary arguments" do
      expect{NGSCI::Read.new(1)}.to raise_error(ArgumentError)
    end

    it "instantiates a new read with proper unstranded arguments" do
      expect{NGSCI::Read.new(1,3)}.to_not raise_error
    end

    it "instantiates a new read with proper stranded arguments" do
      expect{NGSCI::Read.new(1,3,strand:"+")}.to_not raise_error
    end
  end

  context "after created" do
    before(:each) do
      @read = NGSCI::Read.new(1,3,strand:"+")
    end

    it "has a start attribute" do
      expect(@read.methods).to include(:start)
    end
    
    it "has a stop attribute" do
      expect(@read.methods).to include(:stop)
    end
  
    it "has a strand attribute" do
      expect(@read.methods).to include(:strand)
    end

  end
end
