module SCI

  # A simple read class
  #
  # @!attribute [r] start
  # @!attribute [r] stop
  # @!attribute [r] strand
  class Read
    attr_reader :start, :stop, :strand
    def initialize(start,stop,strand: nil)
=begin DEPRECATED chromosome variable
      unless chr.is_a?(String)
        raise SCIError.new "Invalid chromosome argument:\n"
        "chr:#{chr}\tstart:#{start}\tstop:#{stop}\tstrand:#{strand}"
      end
=end
      unless start.is_a?(Integer) && stop.is_a?(Integer) && stop > start
        raise SCIError.new "Invalid coordinate arguments:\n"
        "chr:#{chr}\tstart:#{start}\tstop:#{stop}\tstrand:#{strand}"
      end
      if strand && !%w(+ -).include?(strand)
        raise SCIError.new "Invalid strand argument:\n"
        "chr:#{chr}\tstart:#{start}\tstop:#{stop}\tstrand:#{strand}"
      end
      @start=start
      @stop=stop
      @strand=strand
    end

  end
end
