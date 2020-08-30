# README - NGS-CI
> A Python CLI and module for calculating per-base sequencing complexity.

NOTE: This project is pre-alpha, all of the badge links are broken and are just placeholders at the moment. Development is ongoing. But feel free to clone the repository and play with the code for yourself!

## Development Status

[![PyPI version](https://img.shields.io/pypi/v/ngsci.svg)][pip]
[![Python versions](https://img.shields.io/pypi/pyversions/ngsci.svg)][Pythons]
[![Travis Build Status](https://travis-ci.org/MatthewRalston/ngs-ci.svg?branch=master)](https://travis-ci.org/MatthewRalston/ngsci)
[![Coveralls code coverage](https://img.shields.io/coveralls/MatthewRalston/ngs-ci/master.svg)][Coveralls]
[![ReadTheDocs status](https://readthedocs.org/projects/ngs-ci/badge/?version=latest&style=flat)][RTD]


[pip]: https://pypi.org/project/ngsci/
[Pythons]: https://pypi.org/project/ngsci/
[Coveralls]: https://coveralls.io/r/MatthewRalston/ngs-ci?branch=master
[RTD]: https://ngs-ci.readthedocs.io/en/latest/

## Summary 

The Next Generation Sequencing Complexity Index (NGS-CI) is a Python command-line tool designed for bioinformatics applications. It creates a complexity score on a scale of 0-100 that represents the theoretical percent of possible complexity of reads aligned to a given base.

The principle goal of the library is to generate a single per-base metric to summarize sequencing complexity as a percentage of the total possible unique reads aligned to a given base.



## Installation

OS X and Linux release:

```sh
pip install ngsci
```

Development installation:

```sh
git clone https://github.com/MatthewRalston/ngsci.git
cd ngsci
pip install -r requirements.txt
pip install -r requirements-dev.txt
PYTHONPATH=$(pwd):$PYTHONPATH
```

## Usage Example

CLI Usage

```bash
./bin/ngsci --help
./bin/ngsci -vv input.bam > output.txt
```



## Documentation

Check out the [Readthedocs documentation](https://ngsci.readthedocs.io/en/latest/), with examples and descriptions of the module usage.

## Development

```bash
pytest --cov=ngsci
```

## License

Created by Matthew Ralston - [Scientist, Programmer, Musician](http://matthewralston.github.io) - [Email](mailto:mrals89@gmail.com)

Distributed under the GPL v3.0 license. See `LICENSE.txt` for the copy distributed with this project. Open source software is not for everyone, but for those of us starting out and trying to put the ecosystem ahead of ego, we march into the information age with this ethos.

## Contributing

1. Fork it (<https://github.com/MatthewRalston/ngs-ci/fork>)
2. Create your feature branch (`git checkout -b feature/fooBar`)
3. Commit your changes (`git commit -am 'Add some fooBar'`)
4. Push to the branch (`git push origin feature/fooBar`)
5. Create a new Pull Request

