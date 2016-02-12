[![Build Status](https://travis-ci.org/MatthewRalston/ngs-ci.png?branch=master)](https://travis-ci.org/MatthewRalston/ngs-ci)

[![Gem Version](https://badge.fury.io/rb/ngs-ci.png)](http://badge.fury.io/rb/ngs-ci)

[![Coverage Status](https://coveralls.io/repos/MatthewRalston/ngs-ci/badge.svg?branch=master&service=github)](https://coveralls.io/github/MatthewRalston/ngs-ci?branch=master)


# Todo

The inconsistency between the max summed dissimilarity and the denominator calculation is likely:
1. Add option to calculate the denominator as the sum of all read lengths in the group instead of the number of reads times the read length if a flag is provided
2. Issue in the complexity index (when present is max, 
    a. ( present - missing ) / max = present/max_similarity - missing/max_dissim

# NGS Complexity Index

NOTE: This is a project in progress. 
This gem will calculate a sequencing complexity index for BAM files.

## Installation

Add this line to your application's Gemfile:

```ruby
gem 'NGSCI'
```

And then execute:

    $ bundle

Or install it yourself as:

    $ gem install ngs-ci --pre

Or install manually:

    $ git clone https://github.com/MatthewRalston/ngs-ci.git
    $ cd ngs-ci
    $ gem build ngs-ci.gemspec
    $ gem install ngs-ci-[Version].gem

## Usage

* See ```--help``` for details. More to come.

## Contributing

1. Fork it ( https://github.com/MatthewRalston/SCI/fork )
2. Create your feature branch (`git checkout -b my-new-feature`)
3. Commit your changes (`git commit -am 'Add some feature'`)
4. Push to the branch (`git push origin my-new-feature`)
5. Create a new Pull Request

## License
GPL v3. See LICENSE.txt for details.
