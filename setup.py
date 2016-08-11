# -*- coding: utf-8 -*-
"""
setup file for creating python package pjInvert
"""
from setuptools import setup, find_packages
compile_args = ['-O3']
setup(name='gpsMSSA',
      version='0.1',
      description='Multi Channel Singular Spectrum Analysis',
      url='https://github.com/USFgeodesy/gpsMSSA.git',
      author='Nick Voss',
      author_email= 'nvoss@mail.usf.edu',
      packages = ['gpsMSSA'],
      zip_safe = False,
      )
