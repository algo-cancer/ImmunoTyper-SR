from setuptools import setup, find_packages
import os

setup(
    name="immunotyper",
    packages=['immunotyper'],
    package_dir={'':'src'},
    version='0.1.0',
    description='IGHV genotyper tool using PacBio WGS Long Reads',
    author='Michael Ford',
    license='GNU General Public License v3 (GPL-3) modified by the Commons Clause',
)
