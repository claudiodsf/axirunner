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

## Citing Axitra
If you used this package or the original Axitra for a scientific publication,
please cite the following paper:

> Cotton, F., & Coutant, O. (1997). Dynamic stress variations due to shear
> faults in a plane-layered medium. Geophysical Journal International,
> 128(3), 676–688, doi:
> [10.1111/j.1365-246x.1997.tb05328.x](https://doi.org/10.1111/j.1365-246x.1997.tb05328.x).
