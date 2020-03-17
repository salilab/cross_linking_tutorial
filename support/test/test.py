#!/usr/bin/env python

import unittest
import sys
import os
import subprocess
import tempfile
import shutil

TOPDIR = os.path.abspath(os.path.join(os.path.dirname(__file__), '..', '..'))
os.chdir(TOPDIR)

class Tests(unittest.TestCase):
    def setUp(self):
        self._tmpdir = tempfile.mkdtemp()
        with open(os.path.join(self._tmpdir, 'sitecustomize.py'), 'w') as fh:
            fh.write("""
import matplotlib
matplotlib.use('agg')
""")
        self.pypath = os.environ.get('PYTHONPATH', '')
        os.environ['PYTHONPATH'] = self._tmpdir + os.pathsep + self.pypath

    def tearDown(self):
        os.environ['PYTHONPATH'] = self.pypath
        shutil.rmtree(self._tmpdir, ignore_errors=True)

    def test_xl_ms(self):
        """Test the intro to cross-linking script"""
        p = subprocess.check_call([sys.executable, "cross-link_ms.py"])
        for gen in ('xlinks.csv', 'excluded.None.xl.db',
                    'missing.None.xl.db', 'included.None.xl.db'):
            os.unlink(gen)

    def test_ambiguity(self):
        """Test the cross-linking ambiguity script"""
        p = subprocess.check_call([sys.executable,
                                   "cross-link_ms-ambiguity.py"])
        os.unlink('xlinks.csv')


if __name__ == '__main__':
    unittest.main()
