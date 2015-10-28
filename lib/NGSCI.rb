require 'yell'

# NGSCI stands for Sequencing Complexity Index
# This program calculates a sequencing complexity index for each base and/or strand in a genome.
# This program calculates this by averaging average overlaps of reads aligned to that base.
module NGSCI
  # For custom error handling in the future, unimplemented
  class NGSCIError < StandardError; end
  class NGSCIIOError < NGSCIError; end
  class NGSCIArgError < NGSCIError; end


  # Create the universal logger and include it in Object
  # making the logger object available everywhere
  format = Yell::Formatter.new("[%5L] %d : %m", "%Y-%m-%d %H:%M:%S")
  # http://xkcd.com/1179/
  Yell.new(:format => format) do |l|
    l.level = :info
    l.name = Object
    l.adapter STDOUT, level: [:debug, :info, :warn]
    l.adapter STDERR, level: [:error, :fatal]
  end
  Object.send :include, Yell::Loggable

end # NGSCI

# Integrate modules
require 'NGSCI/cmd'
require 'NGSCI/version'
require 'NGSCI/calculator'
require 'NGSCI/read'
