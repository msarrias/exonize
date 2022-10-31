#!/usr/bin/env python3

from setuptools import setup, find_packages

setup(  # See https://stackoverflow.com/a/58534041 for a description of handled keyword arguments.
    name='lib', 
    maintainer='msarrias',
    maintainer_email='msarrias@math.su.se',
    packages=find_packages(),
    python_requires='==3.8.*',
)
