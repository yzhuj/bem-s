#!/usr/bin/python
# -*- coding: utf8 -*-
#
#   bem: triangulation and fmm/bem electrostatics tools 
#
#   Copyright (C) 2011-2012 Robert Jordens <jordens@gmail.com>
#
#   This program is free software: you can redistribute it and/or modify
#   it under the terms of the GNU General Public License as published by
#   the Free Software Foundation, either version 3 of the License, or
#   (at your option) any later version.
#
#   This program is distributed in the hope that it will be useful,
#   but WITHOUT ANY WARRANTY; without even the implied warranty of
#   MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
#   GNU General Public License for more details.
#
#   You should have received a copy of the GNU General Public License
#   along with this program.  If not, see <http://www.gnu.org/licenses/>.

# Hack to prevent stupid error on exit of `python setup.py test`. (See
# http://www.eby-sarna.com/pipermail/peak/2010-May/003357.html.)
# https://github.com/erikrose/more-itertools/commit/da7e3c771523711adeaef3c6a67ba99de5e2e81a

# This setup uses cythonize(), which is an updated cython buildin function.
# See https://cython.readthedocs.io/en/latest/src/changes.html#id42.  wwc
try:
    import multiprocessing  # parallel compilation 
except ImportError:
    pass

try:
    from setuptools import setup, Extension, find_packages
except ImportError:
    from distutils.core import setup    # distutils
    from distutils.extension import Extension

from Cython.Build import cythonize
import numpy
import os

setup(
    name="pyfastlap",
    description="BEM FMM Laplace solver",
    long_description= """Python bindings for Fastlap""",
    version="0.0+dev",
    author="Robert Jordens",
    author_email="jordens@gmail.com",
    url="http://launchpad.net/pyfastlap",
    license="multiple",
    install_requires=["numpy", "mayavi", "cython"],
    packages = find_packages(),
    test_suite = "bem.tests",
    # zip_safe = False, # See https://cython.readthedocs.io/en/latest/src/reference/compilation.html#configuring-the-c-build
    ext_modules = cythonize([
        Extension("bem.fastlap",
            define_macros = [],
            # extra_compile_args=["-ffast-math"],
            sources = [
                "bem/fastlap.pyx",
                "bem/fastlap_support.c",
                "fastlap/fastlap.c",
                "fastlap/calcp.c",
                "fastlap/direct.c",
                "fastlap/memtracker.c",
                "fastlap/mulDisplay.c",
                "fastlap/mulDo.c",
                "fastlap/mulMats.c",
                "fastlap/mulGlobal.c",
                "fastlap/mulMulti.c",
                "fastlap/mulLocal.c",
                "fastlap/mulSetup.c",],
            include_dirs = [
                "fastlap",
                numpy.get_include(),],
        ),
        
        Extension("bem.pytriangle",
            define_macros = [
                ("TRILIBRARY", "1"),
                ("NO_TIMER", "1"),
                # ("REDUCED", "1"),
                ("REAL", "double"),
                ("EXTERNAL_TEST", "1"),],
            # extra_compile_args=["-ffast-math"],
            sources = [
                "bem/pytriangle.pyx",
                "triangle/triangle.c",],
            include_dirs = [
                "triangle",
                numpy.get_include(),],
        ),
    ]),
)