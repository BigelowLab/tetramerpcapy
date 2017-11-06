
import os
import io
from setuptools import setup


def read(fname):
    return open(os.path.join(os.path.dirname(__file__), fname)).read()

def get_version(relpath):
    '''Read version info from a file without importing it'''
    for line in io.open(os.path.join(os.path.dirname(__file__), relpath), encoding='cp437'):
        if '__version__' in line:
            if '"' in line:
                return line.split('"')[1]
            elif "'" in line:
                return line.split("'")[1]

setup(
    name = "tetramerpca",
    version = get_version('tetramerpca.py'),
    author = "Ben Tupper",
    author_email = "btupper@bigelow.org",
    description = ("Tetramer analysis with PCA"),
    py_modules=[
        'tetramerpca',
        'tabulate',
        'misc',
        'lut',
        'inout',
        'blast',
        'draw'],
    license = "MIT",
    install_requires=[
        'click',
        'pandas',
        'Bio',
        'matplotlib',
        'adjustText'],
    entry_points='''
        [console_scripts]
        tetramerpca=tetramerpca:main
        '''
) 