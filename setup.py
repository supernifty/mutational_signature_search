#!/usr/bin/env python

import setuptools
import pathlib

name = "mutational_signature_search"
version = "0.1"
release = "0.1.0"
here = pathlib.Path(__file__).parent.resolve()

setuptools.setup(
    name=name,
    version=version,
    install_requires=[ 'matplotlib', 'pip @ https://github.com/supernifty/mutational_signature'],
    packages=setuptools.find_packages()
)
