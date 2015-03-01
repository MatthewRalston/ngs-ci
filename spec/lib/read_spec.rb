require 'spec_helper'


describe "reads" do

  it "fails to instantiate on a string start site" do
    expect{LCI::Read.new("foo",3)}.to raise_error(LCI::LCIError)
  end

  it "fails to instantiate on a string stop site" do
    expect{LCI::Read.new(1,"foo")}.to raise_error(LCI::LCIError)
  end

  it "fails to instantiate when the stop site is greater than the start site" do
    expect{LCI::Read.new(3,1)}.to raise_error(LCI::LCIError)
  end

  it "fails to instantiate on an improper strand argument" do
    expect{LCI::Read.new(1,3,strand:"foo")}.to raise_error(LCI::LCIError)
  end

  it "fails to instantiate without the three necessary arguments" do
    expect{LCI::Read.new(1)}.to raise_error(ArgumentError)
  end

  it "instantiates a new read with proper unstranded arguments" do
    expect{LCI::Read.new(1,3)}.to_not raise_error
  end

  it "instantiates a new read with proper stranded arguments" do
    expect{LCI::Read.new(1,3,strand:"+")}.to_not raise_error
  end


end
