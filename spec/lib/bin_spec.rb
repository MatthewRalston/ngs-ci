require 'spec_helper'

testfasta="spec/test_files/test.fa"
testbam="spec/test_files/test.bam"

describe "bin executable" do

  it "runs the --help option" do
    c=SCI::Cmd.new("bundle exec bin/sci --help")
    c.run
    expect(c.stdout).to match /DESCRIPTION/
    expect(c.status).to be_success
  end
  context "non-existent input files" do
    it "fails on a non-existent bam file" do
      c=SCI::Cmd.new("bundle exec bin/sci --bam foo.bam --reference #{testfasta}")
      c.run
      expect(c.status).to_not be_success
    end

    it "fails on a non-existent fasta file" do
      c=SCI::Cmd.new("bundle exec bin/sci --bam #{testbam} --reference foo.fasta")
      c.run
      expect(c.status).to_not be_success
    end
  end
  context "improper options" do
    it "fails on improper strand-specific argument" do
      c=SCI::Cmd.new("bundle exec bin/sci --bam #{testbam} --reference #{testfasta} --strand G")
      c.run
      expect(c.status).to_not be_success
    end
    it "fails on improper loglevel argument" do
      c=SCI::Cmd.new("bundle exec bin/sci --bam #{testbam} --reference #{testfasta} --loglevel foo")
      c.run
      expect(c.status).to_not be_success
    end    
  end
  it "runs on test data" do
    c=SCI::Cmd.new("bundle exec bin/sci --bam #{testbam} --reference #{testfasta} --strand F")
    c.run
    expect(c.status).to be_success
  end

  it "produces speed information" do
    c=SCI::Cmd.new("bundle exec bin/sci --bam #{testbam} --reference #{testfasta} --loglevel debug")
    c.run
    expect(c.stdout).to include("Measure Mode:")
    
  end
end
