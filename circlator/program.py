import shutil
import os
import re
import subprocess
from distutils.version import LooseVersion
from circlator import common


class Program:
    def __init__(self, name, version_cmd, version_regex, environment_var=None, debug=False):
        self.name = name
        self.debug = debug

        if environment_var is not None and environment_var in os.environ:
            if self.debug:
                print(self.name, '- getting path from environment variable', environment_var, '=', os.environ[environment_var], flush=True)
            self.path = os.environ[environment_var]
        else:
            if self.debug:
                print(self.name, '- not using environment variable', flush=True)
            self.path = name


        self.version_cmd = version_cmd
        self.version_regex = version_regex

        if self.debug:
            print(self.name, '- checking which(' + self.path + ')', flush=True)

        self.from_which = shutil.which(self.path)

        if self.debug:
            print('   ... got: "', self.from_which, '"', sep='', flush=True)

        self._set_version()


    def in_path(self):
        '''Returns true iff it is in the path'''
        return self.from_which is not None


    def _set_version(self):
        self.version = None

        if self.debug:
            print(self.name, '- checking version ...')

        if not self.in_path():
            if self.debug:
                print(' ... not in path so cannot get version', flush=True)
                return

        cmd = self.exe() + ' ' + self.version_cmd
        if self.debug:
            print('Running this command to get version:', cmd)
        cmd_output = subprocess.Popen(cmd, shell=True, stdout=subprocess.PIPE, stderr=subprocess.PIPE).communicate()
        cmd_output = common.decode(cmd_output[0]).split('\n')[:-1] + common.decode(cmd_output[1]).split('\n')[:-1]
        if self.debug:
            print('__________ (begin output from ', cmd, ')___________', sep='')
            print(*cmd_output, sep='\n')
            print('__________ (end of output from ', cmd, ')___________', sep='')
            print('Looking far a match to the regex "', self.version_regex.pattern, '" in the above output', sep='', flush=True)
        for line in cmd_output:
            hits = self.version_regex.search(line)
            if hits:
                if self.debug:
                    print('Match to this line:', line)
                    print('Got version:', hits.group(1), flush=True)
                self.version = hits.group(1)
                break
        else:
            if self.debug:
                print('No match found to the regex', flush=True)



    def version_at_least(self, min_version):
        v = self.version
        if v is None:
            return None
        return LooseVersion(v) >= LooseVersion(min_version)


    def version_at_most(self, max_version):
        v = self.version
        if v is None:
            return None
        return LooseVersion(v) <= LooseVersion(max_version)


    def exe(self):
        '''Returns exectuable that can be used in system call'''
        return self.path

