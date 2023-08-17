#!/usr/bin/env python
from setuptools import setup, find_packages
from os import path

def find_api(name):
    return f"{name} @ file://localhost/{here}/../{name}#egg={name}"

here = path.abspath(path.dirname(__file__))

with open(path.join(here, '..', 'version.py')) as f:
    exec(f.read())

with open(path.join(here, 'requirements.txt')) as f:
    requirements = f.read().split()

requirements += [find_api('microns-proximities-api')]

setup(
    name='microns-proximities',
    version=__version__,
    description='EM/Functional proximities for MICrONS',
    author='Andreas Tolias Lab',
    author_email='astolias@bcm.edu',
    packages=find_packages(exclude=[]),
    install_requires=requirements
)