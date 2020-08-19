from setuptools import find_packages, setup

setup(
    name='immunotyper',
    packages=find_packages("src"),
    package_dir={"": "src"},
    version='0.1.0',
    description='IGHV genotyper tool using PacBio WGS Long Reads',
    author='Michael Ford',
    license='GNU General Public License v3 (GPL-3) modified by the Commons Clause',
)
