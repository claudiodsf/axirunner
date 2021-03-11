#!/usr/bin/env python
# -*- coding: utf8 -*-
"""
Wrapper for running Axitra

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
from .version import get_git_version
from .configobj import ConfigObj
from .configobj.validate import Validator


def err_exit(msg):
    msg = str(msg)
    sys.stderr.write(msg + '\n')
    sys.exit(1)


def compute_nfreqs(config):
    npts = config.time_length * config.sampling_rate
    # round to closest power of 2
    npts = int(2**np.around(np.log2(npts)))
    # recompute time_length
    config.time_length = npts/config.sampling_rate
    config.number_of_frequencies = int(npts/2)


def write_axi_data(config):
    axi_data = '&input\n'
    axi_data +=\
        'nc={},nfreq={},tl={},aw={},nr={},ns={},xl={},ikmax={},\n'.format(
            len(config.layers),
            config.number_of_frequencies,
            config.time_length,
            config.imaginary_freq_coefficient,
            len(config.stations),
            len(config.sources),
            config.medium_periodicity,
            config.max_iterations
        )
    free_surface = '.true.' if config.free_surface else '.false'
    axi_data += 'latlon=.false.,freesurface={},'.format(free_surface)
    axi_data += 'sourcefile="source.xyz",statfile="station.xyz"\n//\n'
    # From AXITRA README:
    # (2) in free format and for each layer:
    # - thickness (or depth of the upper interface), vp, vs, rho, qp,qs
    #   Units are [m, m/s and Kg/m^3], or [km, km/s and Kg/km^3]
    #   If you specify the upper interface depth, it must be set to 0. for the
    #   free surface. !!!! rho unit must be consistent with the length unit
    #   (m or km) !!!!
    for layer in config.layers:
        layer[0] *= 1e3  # km --> m;
        layer[1] *= 1e3  # km/s --> m/s;
        layer[2] *= 1e3  # km/s --> m/s;
        layer[3] *= 1e3  # g/cm3 --> kg/m3
        line = '{} {} {} {} {} {}'.format(*layer[:6])
        axi_data += line + '\n'
    outfile = os.path.join(config.run_name, 'axi.data')
    with open(outfile, 'w') as fp:
        fp.write(axi_data)


def write_stations(config):
    outfile = os.path.join(config.run_name, 'station.xyz')
    with open(outfile, 'w') as fp:
        for n, station in enumerate(config.stations_xyz):
            line = '{:10d} {:14.3f} {:14.3f} {:14.3f}'.format(
                n+1, *station[1:4])
            fp.write(line + '\n')


def write_sources(config):
    xyz_file = os.path.join(config.run_name, 'source.xyz')
    axi_hist = os.path.join(config.run_name, 'axi.hist')
    fp_xyz = open(xyz_file, 'w')
    fp_hist = open(axi_hist, 'w')
    for n, source in enumerate(config.sources):
        line = '{:10d} {:14.3f} {:14.3f} {:14.3f}'.format(
            n+1, source.x, source.y, source.z)
        fp_xyz.write(line + '\n')
        line = '{} {:.1e} {} {} {} 0. 0. {}'.format(
            n+1, source.moment, source.strike, source.dip, source.rake,
            source.offset-config.trace_start_offset)
        fp_hist.write(line + '\n')
    fp_xyz.close()
    fp_hist.close()


def read_model_file(config):
    try:
        fp = open(config.velocity_model_file, 'r')
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
    config.layers = layers


def read_station_file(config):
    try:
        fp = open(config.stations_file, 'r')
    except Exception as err:
        err_exit(err)
    stations = []
    for line in fp.readlines():
        words = line.split()
        if words[0] == '#':
            continue
        station = [words[0]] + [float(w) for w in words[1:4]]
        stations.append(station)
    fp.close()
    config.stations = stations


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


def read_source_file(config):
    try:
        fp = open(config.sources_file, 'r')
    except Exception as err:
        err_exit(err)
    sources = []
    for line in fp.readlines():
        if line.strip().startswith('#'):
            continue
        sources.append(Source(line))
    fp.close()
    config.min_origin_time = min(s.origin_time for s in sources)
    for s in sources:
        s.offset = s.origin_time - config.min_origin_time
    config.sources = sources


def project(config):
    if not config.geographical_coordinates:
        config.stations_xyz = config.stations
        return
    lats_sta = np.array([s[1] for s in config.stations])
    lons_sta = np.array([s[2] for s in config.stations])
    lats_src = np.array([s.lat for s in config.sources])
    lons_src = np.array([s.lon for s in config.sources])
    lats = np.hstack((lats_sta, lats_src))
    lons = np.hstack((lons_sta, lons_src))
    lon0 = np.mean(lons)
    lat0 = np.mean(lats)
    lat1 = np.floor(np.min(lats))
    lat2 = np.ceil(np.max(lats))
    p = Proj(proj='lcc', lat_0=lat0, lon_0=lon0, lat_1=lat1, lat_2=lat2)
    stations_xyz = []
    for station in config.stations:
        # projection output is in meters
        # AXITRA convention is x=north and y=east
        y, x = p(station[2], station[1])
        # z is in meters
        z = station[3]
        stations_xyz.append([station[0], x, y, z])
    config.stations_xyz = stations_xyz
    # sources_xyz = []
    for source in config.sources:
        # projection output is in meters
        # AXITRA convention is x=north and y=east
        source.y, source.x = p(source.lon, source.lat)
    # config.sources_xyz = sources_xyz
    config.projection = p


def run_axitra(config):
    if config.axitra_path != 'None':
        axitra = config.axitra_path
    else:
        axitra = 'axitra'
    cmd = 'cd {} && {} && cd ..'.format(config.run_name, axitra)
    os.system(cmd)


def run_convms(config):
    if config.convms_path != 'None':
        convms = config.convms_path
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
    convms_file = os.path.join(config.run_name, 'convms.in')
    with open(convms_file, 'w') as fp:
        fp.write('{}\n'.format(sf[config.source_function]))
        if config.source_function in sf_with_duration:
            fp.write('{}\n'.format(config.source_duration))
        fp.write('{}\n'.format(output[config.output]))
    cmd = 'cd {} && {} < convms.in && cd ..'.format(config.run_name, convms)
    os.system(cmd)


def update_traces(config):
    # AXITRA convention is x=north and y=east
    cmp = {'X': 'N', 'Y': 'E', 'Z': 'Z'}
    instr = {'velocity': 'H', 'displacement': 'H', 'acceleration': 'N'}
    for n, station in enumerate(config.stations):
        stid = station[0].split('.')
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
        filenames = os.path.join(config.run_name, filenames)
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
            channel = band + instr[config.output] + cmp[tr.stats.channel]
            tr.stats.channel = channel
            tr.stats.sac.stla = station[1]
            tr.stats.sac.stlo = station[2]
            tr.stats.sac.stel = station[3]
            sac = SACTrace.from_obspy_trace(tr)
            sac.reftime = config.min_origin_time+config.trace_start_offset
            sac.b = 0
            outfile = tr.id + '.sac'
            outfile = os.path.join(config.run_name, outfile)
            sac.write(outfile)
            os.remove(fname)


# USER INTERFACE --------------------------------------------------------------
class dotdict(dict):
    """dot.notation access to dictionary attributes."""
    __getattr__ = dict.get
    __setattr__ = dict.__setitem__
    __delattr__ = dict.__delitem__


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
        version='%(prog)s {}'.format(get_git_version()))
    args = parser.parse_args()
    return args
# END: USER INTERFACE ---------------------------------------------------------


def main():
    args = parse_arguments()
    configspec = parse_configspec()
    if args.sampleconfig:
        write_sample_config(configspec, 'run_axitra')
        sys.exit(0)
    config = read_config(args.configfile, configspec)
    validate_config(config)
    config = dotdict(config)
    config.args = args
    os.makedirs(config.run_name, exist_ok=True)
    compute_nfreqs(config)
    read_model_file(config)
    read_station_file(config)
    read_source_file(config)
    project(config)
    write_axi_data(config)
    write_stations(config)
    write_sources(config)
    run_axitra(config)
    run_convms(config)
    update_traces(config)
