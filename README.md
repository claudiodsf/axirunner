# PyAxitra
> Python bindings to Axitra
>
> (c) 2020 Claudio Satriano

## Installation
Clone or donwload this repository, then from within the main repository directory, run:

    pip install .

You can also install in editable mode (for developers), with:

    pip install -e .

## Running
To generate a sample config file:

    run_axitra -s

The sample config file should be (hopefully) self-explanatory.

To run, using a config file:

    run_axitra -c <CONFIGFILE>

To get help:

    run_axitra -h

## Limitations
Only the "moment" version of Axitra is currently supported and only for
double-couple sources.

## Getting Axitra
Take a look [here](https://github.com/coutanto/axitra).
