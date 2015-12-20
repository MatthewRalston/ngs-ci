
# coding: utf-8
lib = File.expand_path('../lib', __FILE__)
$LOAD_PATH.unshift(lib) unless $LOAD_PATH.include?(lib)
require 'NGSCI/version'

Gem::Specification.new do |spec|
  spec.name          = "ngs-ci"
  spec.version       = NGSCI::VERSION
  spec.authors       = ["Matthew Ralston"]
  spec.email         = ["mrals89@gmail.com"]
  spec.summary       = %q{Next Generation Sequencing Complexity Index.}
  spec.description   = %q{Calculated a metric that estimates read complexity at each base for RNA-seq BAM files. Alternative to pileup format.}
  spec.homepage      = ""
  spec.license       = "GPL v3"

  spec.files         = `git ls-files -z`.split("\x0")
  spec.executables   = spec.files.grep(%r{^bin/}) { |f| File.basename(f) }
  spec.test_files    = spec.files.grep(%r{^(test|spec|features)/})
  spec.require_paths = ["lib"]

  spec.add_dependency 'trollop','~> 2.1.2'
  spec.add_dependency 'bio-samtools', '= 2.3.2'
  spec.add_dependency 'parallel', '~> 1.4'
  spec.add_dependency 'yell', '~> 0'
  spec.add_dependency "ruby-prof", "~> 0.15"
  spec.has_rdoc = 'yard', '~> 0'

  spec.add_development_dependency "bundler", "~> 0"
  spec.add_development_dependency "rake", "~> 10.0"
  spec.add_development_dependency "rspec", "~> 3.1"
  #spec.add_development_dependency "guard", "~> 2.12"
  spec.add_development_dependency "coveralls", "~> 0"
  #spec.add_development_dependency "cucumber",  "~> 1.3"
end
