#!/usr/bin/env python3
# -*- coding: utf-8 -*-

from setuptools import find_packages, setup

__author__ = 'Yuto Kimura'
__version__ = '0.0.5'

with open('README.md') as f:
    readme = f.read()

with open('LICENSE') as f:
    license = f.read()

setup(
    name='nsdm',
    version=__version__,
    description='NGS Software Development Module',
    long_description=readme,
    author=__author__,
    author_email='biostudy28@gmail.com',
    url='https://github.com/CompBio-TDU-Japan/nsdm',
    license=license,
    packages=find_packages(),
    install_requires=["PyVCF", "biopython"],  # 依存ライブラリ
)
