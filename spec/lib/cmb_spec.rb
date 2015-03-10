require 'spec_helper'


describe "command" do

  it "runs commands" do
    c=SCI::Cmd.new("echo success")
    c.run
    expect(c.stdout.chomp).to eq("success")
  end

  it "receives commands" do
    c=SCI::Cmd.new("echo success")
    expect(c.to_s).to eq("echo success")    
  end

end
