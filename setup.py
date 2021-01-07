# -*- coding: utf-8 -*-
import numpy
from setuptools import Extension, setup
from Cython.Build import cythonize

description = ['Yet Another Filtering Library']

extensions = [
    Extension('yaflpy', ['src/yaflpy.pyx'],
        include_dirs=[numpy.get_include(), 'src', 'src/configpy'],
        define_macros=[("NPY_NO_DEPRECATED_API", "NPY_1_7_API_VERSION")])
    ]

for e in extensions:
    e.cython_directives = {'language_level': "3"} 

setup(
      name="yaflpy",
      version = '0.1.0',
      license = 'Apache License, Version 2.0',
      author = 'anonimous',
      author_email = 'shkolnick-kun@gmail.com',
      url = 'https://github.com/shkolnick-kun/yafl',
      classifiers = [
        'Development Status :: 3 - Alpha',
        'Intended Audience :: Science/Research',
        'License :: OSI Approved :: Apache Software License',
        'Programming Language :: Python :: 3',
        'Programming Language :: Python :: 3.3',
        'Programming Language :: Python :: 3.4',
        'Programming Language :: Python :: 3.5',
        'Programming Language :: Python :: 3.6',
        'Programming Language :: Python :: 3.7',
        'Programming Language :: Python :: 3.8',
        'Programming Language :: Python :: 3.9',
        'Topic :: Software Development :: Embedded Systems',
        ],
      #packages = ['src'],
      ext_modules = cythonize(extensions,
                              compiler_directives={'language_level' : "3"})
)
