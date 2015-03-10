
# coding: utf-8
lib = File.expand_path('../lib', __FILE__)
$LOAD_PATH.unshift(lib) unless $LOAD_PATH.include?(lib)
require 'SCI/version'

Gem::Specification.new do |spec|
  spec.name          = "SCI"
  spec.version       = SCI::VERSION
  spec.authors       = ["Matthew Ralston"]
  spec.email         = ["mrals89@gmail.com"]
  spec.summary       = %q{Sequencing Complexity Index.}
  spec.description   = %q{Calculated a metric that estimates read complexity at each base for RNA-seq BAM files. Alternative to pileup format.}
  spec.homepage      = ""
  spec.license       = "GPL v3"

  spec.files         = `git ls-files -z`.split("\x0")
  spec.executables   = spec.files.grep(%r{^bin/}) { |f| File.basename(f) }
  spec.test_files    = spec.files.grep(%r{^(test|spec|features)/})
  spec.require_paths = ["lib"]

  spec.add_dependency 'trollop','~> 2.0'
  spec.add_dependency 'bio-samtools', '~> 2.3'
  spec.add_dependency 'parallel', '~> 1.4'
  spec.add_dependency 'yell'

  spec.has_rdoc = 'yard'

  spec.add_development_dependency "bundler"
  spec.add_development_dependency "rake", "~> 10.0"
  spec.add_development_dependency "rspec", "~> 3.1"
  spec.add_development_dependency "guard", "~> 2.12"
  spec.add_development_dependency "coveralls"
  spec.add_development_dependency "ruby-prof", "~> 0.15"
  spec.add_development_dependency "cucumber",  "~> 1.3"
end
