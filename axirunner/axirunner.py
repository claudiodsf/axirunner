#!/usr/bin/env python
# -*- coding: utf8 -*-
"""
Wrapper for running Axitra.

:copyright:
    2020-2021 Claudio Satriano <satriano@ipgp.fr>
:license:
    CeCILL Free Software License Agreement, Version 2.1
    (http://www.cecill.info/index.en.html)
"""
import sys
import os
import argparse
import numpy as np
from glob import glob
from shutil import copyfile
from obspy import read, UTCDateTime
from obspy.io.sac import SACTrace
from pyproj import Proj
from ._version import get_versions
from .configobj import ConfigObj
from .configobj.validate import Validator


def err_exit(msg):
    msg = str(msg)
    sys.stderr.write(msg + '\n')
    sys.exit(1)


class Source():
    """An axitra source."""

    def __init__(self, line):
        word = line.split()
        # we don't know, at this stage, wether source is in cartesian...
        self.x = float(word[0])
        self.y = float(word[1])
        # ...or geographical coordinates
        self.lat = float(word[0])
        self.lon = float(word[1])
        self.z = float(word[2])
        self.moment = float(word[3])
        self.strike = float(word[4])
        self.dip = float(word[5])
        self.rake = float(word[6])
        self.origin_time = UTCDateTime(word[7])
        self.offset = 0.


class Station():
    """An axitra station."""

    def __init__(self, line):
        word = line.split()
        self.code = word[0]
        # we don't know, at this stage, wether source is in cartesian...
        self.x = float(word[1])
        self.y = float(word[2])
        # ...or geographical coordinates
        self.lat = float(word[1])
        self.lon = float(word[2])
        self.z = float(word[3])


class Axirunner():
    """Wrapper class for configuring and running axitra and convm."""

    def __init__(self, config):
        for key in config:
            setattr(self, key, config[key])
        os.makedirs(self.run_name, exist_ok=True)
        self.compute_nfreqs()
        self.read_model_file()
        self.read_station_file()
        self.read_source_file()
        self.coordinate_projection()
        self.write_axi_data()
        self.write_stations()
        self.write_sources()

    def compute_nfreqs(self):
        npts = self.time_length * self.sampling_rate
        # round to closest power of 2
        npts = int(2**np.around(np.log2(npts)))
        # recompute time_length
        self.time_length = npts/self.sampling_rate
        self.number_of_frequencies = int(npts/2)

    def read_model_file(self):
        try:
            fp = open(self.velocity_model_file, 'r')
        except Exception as err:
            err_exit(err)
        layers = []
        skip = False
        for line in fp.readlines():
            words = line.split()
            if line.strip().startswith('//'):
                continue
            if line.strip().startswith('#'):
                continue
            if line.strip().startswith('/*'):
                skip = True
            if line.strip().endswith('*/'):
                skip = False
                continue
            if skip:
                continue
            if len(words) < 6:
                continue
            try:
                layer = [float(w) for w in words]
            except Exception:
                continue
            layers.append(layer)
        fp.close()
        self.layers = layers

    def read_station_file(self):
        try:
            fp = open(self.stations_file, 'r')
        except Exception as err:
            err_exit(err)
        stations = []
        for line in fp.readlines():
            if line.strip().startswith('#'):
                continue
            stations.append(Station(line))
        fp.close()
        self.stations = stations

    def read_source_file(self):
        try:
            fp = open(self.sources_file, 'r')
        except Exception as err:
            err_exit(err)
        sources = []
        for line in fp.readlines():
            if line.strip().startswith('#'):
                continue
            sources.append(Source(line))
        fp.close()
        self.min_origin_time = min(s.origin_time for s in sources)
        for s in sources:
            s.offset = s.origin_time - self.min_origin_time
        self.sources = sources

    def coordinate_projection(self):
        if not self.geographical_coordinates:
            return
        lats_sta = np.array([s.lat for s in self.stations])
        lons_sta = np.array([s.lon for s in self.stations])
        lats_src = np.array([s.lat for s in self.sources])
        lons_src = np.array([s.lon for s in self.sources])
        lats = np.hstack((lats_sta, lats_src))
        lons = np.hstack((lons_sta, lons_src))
        lon0 = np.mean(lons)
        lat0 = np.mean(lats)
        lat1 = np.floor(np.min(lats))
        lat2 = np.ceil(np.max(lats))
        p = Proj(proj='lcc', lat_0=lat0, lon_0=lon0, lat_1=lat1, lat_2=lat2)
        for station in self.stations:
            # projection output is in meters
            # AXITRA convention is x=north and y=east
            station.y, station.x = p(station.lon, station.lat)
        for source in self.sources:
            # projection output is in meters
            # AXITRA convention is x=north and y=east
            source.y, source.x = p(source.lon, source.lat)
        self.projection = p

    def write_axi_data(self):
        axi_data = '&input\n'
        axi_data +=\
            'nfreq={},tl={},aw={},xl={},ikmax={},\n'.format(
                self.number_of_frequencies,
                self.time_length,
                self.imaginary_freq_coefficient,
                self.medium_periodicity,
                self.max_iterations
            )
        free_surface = '.true.' if self.free_surface else '.false'
        axi_data += 'latlon=.false.,freesurface={},'.format(free_surface)
        axi_data += 'sourcefile="source.xyz",statfile="station.xyz"\n//\n'
        # From AXITRA README:
        # (2) in free format and for each layer:
        # - thickness (or depth of the upper interface), vp, vs, rho, qp,qs
        #   Units are [m, m/s and Kg/m^3], or [km, km/s and Kg/km^3]
        #   If you specify the upper interface depth, it must be set to 0. for
        #   the free surface. !!!! rho unit must be consistent with the length
        #   unit (m or km) !!!!
        for layer in self.layers:
            layer[0] *= 1e3  # km --> m;
            layer[1] *= 1e3  # km/s --> m/s;
            layer[2] *= 1e3  # km/s --> m/s;
            layer[3] *= 1e3  # g/cm3 --> kg/m3
            line = '{} {} {} {} {} {}'.format(*layer[:6])
            axi_data += line + '\n'
        outfile = os.path.join(self.run_name, 'axi.data')
        with open(outfile, 'w') as fp:
            fp.write(axi_data)

    def write_stations(self):
        outfile = os.path.join(self.run_name, 'station.xyz')
        with open(outfile, 'w') as fp:
            for n, station in enumerate(self.stations):
                line = '{:10d} {:14.3f} {:14.3f} {:14.3f}'.format(
                    n+1, station.x, station.y, station.z)
                fp.write(line + '\n')

    def write_sources(self):
        xyz_file = os.path.join(self.run_name, 'source.xyz')
        axi_hist = os.path.join(self.run_name, 'axi.hist')
        fp_xyz = open(xyz_file, 'w')
        fp_hist = open(axi_hist, 'w')
        for n, source in enumerate(self.sources):
            line = '{:10d} {:14.3f} {:14.3f} {:14.3f}'.format(
                n+1, source.x, source.y, source.z)
            fp_xyz.write(line + '\n')
            line = '{} {:.1e} {} {} {} 0. 0. {}'.format(
                n+1, source.moment, source.strike, source.dip, source.rake,
                source.offset-self.trace_start_offset)
            fp_hist.write(line + '\n')
        fp_xyz.close()
        fp_hist.close()

    def run_axitra(self):
        if self.axitra_path != 'None':
            axitra = self.axitra_path
        else:
            axitra = 'axitra'
        cmd = 'cd {} && {} && cd ..'.format(self.run_name, axitra)
        os.system(cmd)

    def run_convms(self):
        if self.convms_path != 'None':
            convms = self.convms_path
        else:
            convms = 'convms'
        sf = {
            'dirac': 0,
            'ricker': 1,
            'step': 2,
            'triangle': 4,
            'ramp': 5,
            'true_step': 7,
            'trapezoid': 8
        }
        sf_with_duration = ['ricker', 'step', 'triangle', 'ramp', 'trapezoid']
        output = {
            'displacement': 1,
            'velocity': 2,
            'acceleration': 3
        }
        convms_file = os.path.join(self.run_name, 'convms.in')
        with open(convms_file, 'w') as fp:
            fp.write('{}\n'.format(sf[self.source_function]))
            if self.source_function in sf_with_duration:
                fp.write('{}\n'.format(self.source_duration))
            fp.write('{}\n'.format(output[self.output]))
        cmd = 'cd {} && {} < convms.in && cd ..'.format(self.run_name, convms)
        os.system(cmd)

    def update_traces(self):
        # AXITRA convention is x=north and y=east
        cmp = {'X': 'N', 'Y': 'E', 'Z': 'Z'}
        instr = {'velocity': 'H', 'displacement': 'H', 'acceleration': 'N'}
        idep = {'velocity': 7, 'displacement': 6, 'acceleration': 8}
        for n, station in enumerate(self.stations):
            stid = station.code.split('.')
            if len(stid) == 3:
                net, sta, loc = stid
            elif len(stid) == 2:
                net, sta = stid
                loc = None
            else:
                sta = '_'.join(stid)
                net = None
                loc = None
            filenames = 'axi{:03d}.?.sac'.format(n+1)
            filenames = os.path.join(self.run_name, filenames)
            for fname in glob(filenames):
                tr = read(fname)[0]
                tr.stats.network = net
                tr.stats.station = sta
                tr.stats.location = loc
                sr = tr.stats.sampling_rate
                if sr < 10:
                    band = 'L'
                elif 10 <= sr < 80:
                    band = 'B'
                elif sr >= 80:
                    band = 'H'
                channel = band + instr[self.output] + cmp[tr.stats.channel]
                tr.stats.channel = channel
                # Set SAC data type. Note that SAC data amplitude is in nm,
                # nm/s or nm/s/s, so we need to multiply by 1e9
                # https://ds.iris.edu/files/sac-manual/manual/file_format.html
                tr.data *= 1e9
                tr.stats.sac.idep = idep[self.output]
                tr.stats.sac.stla = station.lat
                tr.stats.sac.stlo = station.lon
                tr.stats.sac.stel = -station.z
                sac = SACTrace.from_obspy_trace(tr)
                sac.reftime = self.min_origin_time+self.trace_start_offset
                sac.b = 0
                outfile = tr.id + '.sac'
                outfile = os.path.join(self.run_name, outfile)
                sac.write(outfile)
                os.remove(fname)


# USER INTERFACE --------------------------------------------------------------
def parse_configspec():
    curdir = os.path.dirname(__file__)
    configspec_file = os.path.join(curdir, 'conf', 'configspec.conf')
    configspec = read_config(configspec_file)
    return configspec


def _write_ok(filepath):
    if os.path.exists(filepath):
        ans = input(
            '"{}" already exists. '
            'Do you want to overwrite it? [y/N] '.format(filepath))
        if ans in ['y', 'Y']:
            return True
        else:
            return False
    return True


def write_sample_config(configspec, progname):
    c = ConfigObj(configspec=configspec, default_encoding='utf8')
    val = Validator()
    c.validate(val)
    c.defaults = []
    c.initial_comment = configspec.initial_comment
    c.comments = configspec.comments
    configfile = progname + '.conf'
    if _write_ok(configfile):
        with open(configfile, 'wb') as fp:
            c.write(fp)
        print('Sample config file written to: "{}"'.format(configfile))
    curdir = os.path.dirname(__file__)
    confdir = os.path.join(curdir, 'conf')
    for file in ['sources.conf', 'stations.conf', 'velocity.nd']:
        if _write_ok(file):
            copyfile(os.path.join(confdir, file), file)
            file_type = os.path.splitext(file)[0].rstrip('s')
            print('Sample {} file written to: "{}"'.format(file_type, file))


def read_config(config_file, configspec=None):
    kwargs = dict(
        configspec=configspec, file_error=True, default_encoding='utf8')
    if configspec is None:
        kwargs.update(
            dict(interpolation=False, list_values=False, _inspec=True))
    try:
        config_obj = ConfigObj(config_file, **kwargs)
    except IOError as err:
        err_exit(err)
    except Exception as err:
        msg = 'Unable to read "{}": {}'.format(config_file, err)
        err_exit(msg)
    return config_obj


def validate_config(config_obj):
    val = Validator()
    test = config_obj.validate(val)
    if isinstance(test, dict):
        for entry in test:
            if not test[entry]:
                sys.stderr.write('Invalid value for "%s": "%s"\n' %
                                 (entry, config_obj[entry]))
        sys.exit(1)
    if not test:
        err_exit('No configuration value present!')


def parse_arguments():
    """Parse command line arguments."""
    parser = argparse.ArgumentParser(
        description='Run Axitra.')
    group = parser.add_mutually_exclusive_group(required=True)
    group.add_argument(
        '-c', '--configfile', type=str,
        help='config file for data sources and processing params')
    group.add_argument(
        '-s', '--sampleconfig', default=False, action='store_true',
        required=False,
        help='write sample config file to current directory and exit')
    parser.add_argument(
        '-v', '--version', action='version',
        version='%(prog)s {}'.format(get_versions()['version']))
    args = parser.parse_args()
    return args
# END: USER INTERFACE ---------------------------------------------------------


def main():
    args = parse_arguments()
    configspec = parse_configspec()
    if args.sampleconfig:
        write_sample_config(configspec, 'axirunner')
        sys.exit(0)
    config = read_config(args.configfile, configspec)
    validate_config(config)
    axirunner = Axirunner(config)
    axirunner.run_axitra()
    axirunner.run_convms()
    axirunner.update_traces()
