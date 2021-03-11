# AxiRunner
Make synthetic seismograms using [Axitra]. The easy way™.

(c) 2020-2021 Claudio Satriano

[axitra]:https://github.com/coutanto/axitra

## Installation
Clone or download this repository, then from within the main repository directory, run:

    pip install .

You can also install in editable mode (for developers), with:

    pip install -e .

## Running
To generate a sample set of config files:

    axirunner -s

The sample config files should be (hopefully) self-explanatory.

To run, using a config file:

    axirunner -c <CONFIGFILE>

To get help:

    axirunner -h

## Limitations
Only the "moment" version of Axitra is currently supported and only for
double-couple sources.

## Getting Axitra
Take a look [here](https://github.com/coutanto/axitra).
