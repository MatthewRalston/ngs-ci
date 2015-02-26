require 'yell'

# LCI stands for Library Complexity Index
# This program calculates a library complexity index for each base and/or strand in a genome.
# This program calculates this by averaging average overlaps of reads aligned to that base.
module LCI
  # For custom error handling in the future, unimplemented
  class LCIError < StandardError; end
  class LCIIOError < LCIError; end
  class LCIArgError < LCIError; end


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

end # LCI

# Integrate modules
require 'LCI/cmd'
require 'LCI/version'
require 'LCI/calculator'
require 'LCI/read'
