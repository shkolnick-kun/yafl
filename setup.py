# -*- coding: utf-8 -*-
import numpy
from os.path import join, dirname
from setuptools import Extension, setup
from Cython.Build import cythonize #must be after setuptools

#------------------------------------------------------------------------------
setup_dir = dirname(__file__)
src_dir   = join(setup_dir, 'src')

extensions = [
    Extension('yaflpy', [join(src_dir, 'yaflpy.pyx')],
        include_dirs=[numpy.get_include(), src_dir, join(src_dir, 'configpy')],
        define_macros=[("NPY_NO_DEPRECATED_API", "NPY_1_7_API_VERSION")])
    ]

#------------------------------------------------------------------------------
deps = ['setuptools', 'Cython', 'numpy']

#------------------------------------------------------------------------------
setup(
      name="yaflpy",
      version = '0.1.0',
      description = ['Yet Another Filtering Library'],
      long_description = open(join(setup_dir, 'Readme.md')).read(),
      long_description_content_type = 'text/markdown',
      license = 'Apache License, Version 2.0',
      license_file = join(setup_dir, 'LICENSE'),
      author = 'anonimous',
      author_email = 'shkolnick-kun@gmail.com',
      url = 'https://github.com/shkolnick-kun/yafl',
      classifiers = [
        'Development Status :: 3 - Alpha',
        'Topic :: Software Development :: Embedded Systems',
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
        ],
      setup_requires=deps,
      install_requires=deps,
      platforms = ['any'],
      ext_modules = cythonize(extensions,
                              compiler_directives={'language_level' : '3'})
)
