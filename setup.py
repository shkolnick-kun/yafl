# -*- coding: utf-8 -*-
import numpy
import os
from os.path import dirname, join, splitext
from setuptools import Extension, setup

try:
    from Cython.Build import cythonize
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

extensions = [
    Extension(EXT_NAME, [join(EXT_DIR, EXT_NAME + ".pyx")],
        include_dirs=[numpy.get_include(), SRC_DIR, EXT_DIR],
        define_macros=[("NPY_NO_DEPRECATED_API", "NPY_1_7_API_VERSION")])
    ]

CYTHONIZE = bool(int(os.getenv("CYTHONIZE", 0))) and cythonize is not None

if CYTHONIZE:
    compiler_directives = {"language_level": 3, "embedsignature": True}
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
