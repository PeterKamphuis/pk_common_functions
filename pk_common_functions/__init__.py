# -*- coding: future_fstrings -*-

import pkg_resources
import os
import subprocess

def report_version():
    # Distutils standard  way to do version numbering
    try:
        __version__ = pkg_resources.require("pk_common_functions")[0].version
    except pkg_resources.DistributionNotFound:
        __version__ = "dev"
    # perhaps we are in a github with tags; in that case return describe
    path = os.path.dirname(os.path.abspath(__file__))
    try:
        # work round possible unavailability of git -C
        result = subprocess.check_output(
            'cd %s; git describe --tags' % path, shell=True, stderr=subprocess.STDOUT).rstrip().decode()
    except subprocess.CalledProcessError:
        result = None
    if result != None and 'fatal' not in result:
        # will succeed if tags exist
        return result
    else:
        # perhaps we are in a github without tags? Cook something up if so
        try:
            result = subprocess.check_output(
                'cd %s; git rev-parse --short HEAD' % path, shell=True, stderr=subprocess.STDOUT).rstrip().decode()
        except subprocess.CalledProcessError:
            result = None
        if result != None and 'fatal' not in result:
            return __version__+'-'+result
        else:
            # we are probably in an installed version
            return __version__

def report_branch():
    path = os.path.dirname(os.path.abspath(__file__))
    result = None
    try:
        branches = subprocess.check_output(\
            'cd %s; git branch' % path, shell=True, stderr=subprocess.STDOUT).rstrip().decode().split()
        for i,line in enumerate(branches):
            if line == '*':
                result=branches[i+1]
                break
    except subprocess.CalledProcessError:
        result= None
    return result

__version__ = report_version()
__branch__ = report_branch()

#set_logger_values()
