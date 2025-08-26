from setuptools import setup
from Cython.Build import cythonize

setup(
    name='fastsearch',
    ext_modules=cythonize("cython_modules/fastsearch.pyx"),
    package_dir={'': '.'},
) 