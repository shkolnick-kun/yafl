#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
    This file is based on
    https://github.com/FedericoStra/cython-package-example/blob/master/setup.py
"""
"""
    MIT License

    Copyright (c) 2019 Federico Stra

    Permission is hereby granted, free of charge, to any person obtaining a copy
    of this software and associated documentation files (the "Software"), to deal
    in the Software without restriction, including without limitation the rights
    to use, copy, modify, merge, publish, distribute, sublicense, and/or sell
    copies of the Software, and to permit persons to whom the Software is
    furnished to do so, subject to the following conditions:

        The above copyright notice and this permission notice shall be included in all
        copies or substantial portions of the Software.

    THE SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND, EXPRESS OR
    IMPLIED, INCLUDING BUT NOT LIMITED TO THE WARRANTIES OF MERCHANTABILITY,
    FITNESS FOR A PARTICULAR PURPOSE AND NONINFRINGEMENT. IN NO EVENT SHALL THE
    AUTHORS OR COPYRIGHT HOLDERS BE LIABLE FOR ANY CLAIM, DAMAGES OR OTHER
    LIABILITY, WHETHER IN AN ACTION OF CONTRACT, TORT OR OTHERWISE, ARISING FROM,
    OUT OF OR IN CONNECTION WITH THE SOFTWARE OR THE USE OR OTHER DEALINGS IN THE
    SOFTWARE.
"""
"""
    Copyright 2021 anonimous <shkolnick-kun@gmail.com> and contributors.

    Licensed under the Apache License, Version 2.0 (the "License");
    you may not use this file except in compliance with the License.
    You may obtain a copy of the License at

    http://www.apache.org/licenses/LICENSE-2.0

    Unless required by applicable law or agreed to in writing,
    software distributed under the License is distributed on an "AS IS" BASIS,
    WITHOUT WARRANTIES OR CONDITIONS OF ANY KIND, either express or implied.

    See the License for the specific language governing permissions
"""

import numpy
import os
from os.path import dirname, isfile, join, splitext
from setuptools import Extension, setup

try:
    from Cython.Build import cythonize
    from Cython.Compiler.Main import default_options
except ImportError:
    cythonize = None

# https://cython.readthedocs.io/en/latest/src/userguide/source_files_and_compilation.html#distributing-cython-modules
def no_cythonize(extensions, **_ignore):
    for extension in extensions:
        sources = []
        for sfile in extension.sources:
            path, ext = splitext(sfile)
            if ext in (".pyx", ".py"):
                if extension.language == "c++":
                    ext = ".cpp"
                else:
                    ext = ".c"
                sfile = path + ext
            sources.append(sfile)
        extension.sources[:] = sources
    return extensions

EXT_NAME  = "yaflpy"
SETUP_DIR = dirname(__file__)
SRC_DIR   = join(SETUP_DIR, "src")
EXT_DIR   = join(SRC_DIR, EXT_NAME)

DEFINES = [("NPY_NO_DEPRECATED_API", "NPY_1_7_API_VERSION")]
CT_ENV  = {"YAFLPY_USE_64_BIT":False}

if int(os.getenv("YAFL_64", 0)):
    #Override defaults in 64 bit case
    print('YAFL: The extensuion will be compiled in 64 bit variant!!!')
    CT_ENV   = {"YAFLPY_USE_64_BIT":True}
    DEFINES += [("YAFL_USE_64_BIT",1)]

extensions = [
    Extension(EXT_NAME, [join(EXT_DIR, EXT_NAME + ".pyx")],
              language='c',
              include_dirs=[numpy.get_include(), SRC_DIR, EXT_DIR],
              define_macros=DEFINES
              )
    ]

if not isfile(join(EXT_DIR, EXT_NAME + ".c")):
    CYTHONIZE = cythonize is not None
else:
    CYTHONIZE = bool(int(os.getenv("CYTHONIZE", 0))) and cythonize is not None

if CYTHONIZE:
    default_options["compile_time_env"] = CT_ENV
    compiler_directives = {
        "language_level": 3,
        "embedsignature": True
        }
    extensions = cythonize(extensions, compiler_directives=compiler_directives)
else:
    extensions = no_cythonize(extensions)

with open("requirements.txt") as fp:
    install_requires = fp.read().strip().split("\n")

with open("requirements-dev.txt") as fp:
    dev_requires = fp.read().strip().split("\n")

setup(
    ext_modules=extensions,
    install_requires=install_requires,
    extras_require={
        "dev": dev_requires,
        "docs": ["sphinx", "sphinx-rtd-theme"]
    },
)
