# -*- coding: utf-8 -*-
"""
setup file for creating python package pjInvert
"""
from setuptools import setup
from setuptools import find_packages

setup(name='gpsMSSA',
      version='0.1',
      description='Multi Channel Singular Spectrum Analysis',
      url='https://github.com/USFgeodesy/gpsMSSA.git,
      author='Nick Voss',
      author_email= 'nvoss@mail.usf.edu',
      license='MIT',
      packages=find_packages(),
      zip_safe=False)
