import unittest
import os
from circlator import common

modules_dir = os.path.dirname(os.path.abspath(common.__file__))
data_dir = os.path.join(modules_dir, 'tests', 'data')


class TestCommon(unittest.TestCase):
    def test_check_files_exist(self):
        '''test check_files_exist'''
        file_exists = os.path.join(data_dir, 'common_test_file_exists')
        common.check_files_exist([file_exists])
        with self.assertRaises(common.Error):
            common.check_files_exist([file_exists, 'thisisnotafileandshouldcauseanerror'])
