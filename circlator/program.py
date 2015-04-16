import shutil
import os
import re
import subprocess
from distutils.version import LooseVersion
from circlator import common


class Program:
    def __init__(self, name, version_cmd, version_regex, environment_var=None):
        self.name = name
        if environment_var is not None and environment_var in os.environ:
            self.path = os.environ[environment_var]
        else:
            self.path = name


        self.version_cmd = version_cmd
        self.version_regex = version_regex


    def in_path(self):
        '''Returns true iff it is in the path'''
        return shutil.which(self.path) is not None


    def version(self):
        '''Returns version. If not in path, or in path but can't get the version, returns None'''
        if not self.in_path():
            return None

        cmd = self.exe() + ' ' + self.version_cmd
        cmd_output = subprocess.Popen(cmd, shell=True, stdout=subprocess.PIPE, stderr=subprocess.PIPE).communicate()
        cmd_output = common.decode(cmd_output[0]).split('\n')[:-1] + common.decode(cmd_output[1]).split('\n')[:-1]
        for line in cmd_output:
            hits = self.version_regex.search(line)
            if hits:
                return hits.group(1)
        return None


    def version_at_least(self, min_version):
        v = self.version()
        if v is None:
            return None
        return LooseVersion(v) >= LooseVersion(min_version) 


    def exe(self):
        '''Returns exectuable that can be used in system call'''
        return self.path

