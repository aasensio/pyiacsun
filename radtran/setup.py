#! /usr/bin/env python
# -*- coding: utf-8 -*-
"""
Python bindings for MILNE-EDDINGTON f90 Library.
For more information see: https://github.com/aasensio/pyiacsun
::
    Main Changes in 1.0
    ---------------------
    * Working version
:copyright:
    A. Asensio Ramos
:license:
    GNU General Public License (GPL)
"""
from distutils.ccompiler import CCompiler
from distutils.errors import DistutilsExecError, CompileError
from distutils.unixccompiler import UnixCCompiler
from setuptools import find_packages, setup
from setuptools.extension import Extension

import inspect
import os
import platform
from subprocess import Popen, PIPE
import sys


def _compile(self, obj, src, ext, cc_args, extra_postargs, pp_opts):
    compiler_so = self.compiler_so
    arch = platform.architecture()[0].lower()
    if ext == ".f90":
        if sys.platform == 'darwin' or sys.platform == 'linux':
            compiler_so = ["gfortran"]
            cc_args = ["-O", "-fPIC", "-c", "-ffree-form", "-ffree-line-length-none"]
            # Force architecture of shared library.
            if arch == "32bit":
                cc_args.append("-m32")
            elif arch == "64bit":
                cc_args.append("-m64")
            else:
                print("\nPlatform has architecture '%s' which is unknown to "
                      "the setup script. Proceed with caution\n" % arch)
    try:
        self.spawn(compiler_so + cc_args + [src, '-o', obj] +
                   extra_postargs)
    except DistutilsExecError as msg:
        raise CompileError(msg)
UnixCCompiler._compile = _compile


# Hack to prevent build_ext from trying to append "init" to the export symbols.
class finallist(list):
    def append(self, object):
        return


class MyExtension(Extension):
    def __init__(self, *args, **kwargs):
        Extension.__init__(self, *args, **kwargs)
        self.export_symbols = finallist(self.export_symbols)


def get_libgfortran_dir():
    """
    Helper function returning the library directory of libgfortran. Useful
    on OSX where the C compiler oftentimes has no knowledge of the library
    directories of the Fortran compiler. I don't think it can do any harm on
    Linux.
    """
    for ending in [".3.dylib", ".dylib", ".3.so", ".so"]:
        try:
            p = Popen(['gfortran', "-print-file-name=libgfortran" + ending],
                      stdout=PIPE, stderr=PIPE)
            p.stderr.close()
            line = p.stdout.readline().decode().strip()
            p.stdout.close()
            if os.path.exists(line):
                return [os.path.dirname(line)]
        except:
            continue
        return []

# Monkey patch the compilers to treat Fortran files like C files.
CCompiler.language_map['.f90'] = "c"
UnixCCompiler.src_extensions.append(".f90")

# Milne-Eddington
# Import the version string.
path = os.path.join(os.path.abspath(os.path.dirname(inspect.getfile(inspect.currentframe()))), "sourceMilne")
with open(os.path.join(path, "VERSION.txt"), "rt") as fh:
  VERSION = fh.read().strip()

DOCSTRING = __doc__.strip().split("\n")

lib = MyExtension('milne',
                  libraries=["gfortran"],
                  library_dirs=get_libgfortran_dir(),
                  sources=['sourceMilne/pymilne.pyx', 'sourceMilne/vars.f90', 'sourceMilne/maths.f90', 'sourceMilne/atomic.f90', 'sourceMilne/milne.f90'])

setup_configMilne = dict(
    name='milne',
    version=VERSION,
    description=DOCSTRING[0],
    long_description="\n".join(DOCSTRING[2:]),
    author=' A. Asensio Ramos',
    author_email='aasensio@iac.es',
    url='https://github.com/aasensio/pyiacsun',
    license='GNU General Public License, version 3 (GPLv3)',
    platforms='OS Independent',
    install_requires=[
        'numpy',
    ],
    ext_package='',
    ext_modules=[lib],
    classifiers=[
        'Development Status :: 4 - Beta',
        'Environment :: Console',
        'Intended Audience :: Science/Research',
        'Intended Audience :: Developers',
        "Operating System :: Unix",
        "Operating System :: POSIX",
        "Operating System :: MacOS",
        'License :: OSI Approved :: GNU General Public License v3 (GPLv3)',
        'Programming Language :: Python',
        "Programming Language :: Python :: 3.4",
        "Programming Language :: Python :: Implementation :: CPython",
        'Topic :: Scientific/Engineering',
        'Topic :: Scientific/Engineering :: Physics',
    ],
    keywords=['milne', 'radiative transfer'],
    packages=find_packages(),
    zip_safe=False,
    include_package_data=True
)


path = os.path.join(os.path.abspath(os.path.dirname(inspect.getfile(inspect.currentframe()))), "sourceLTE")
with open(os.path.join(path, "VERSION.txt"), "rt") as fh:
  VERSION = fh.read().strip()

DOCSTRING = __doc__.strip().split("\n")

lib = MyExtension('lte',
                  libraries=["gfortran"],
                  library_dirs=get_libgfortran_dir(),
                  sources=['sourceLTE/pylte.pyx', 'sourceLTE/vars.f90', 'sourceLTE/partition.f90', 'sourceLTE/maths.f90', 'sourceLTE/background.f90', 
                  'sourceLTE/hydros.f90', 'sourceLTE/synth.f90', 'sourceLTE/lte.f90'])
setup_configLTE = dict(
    name='lte',
    version=VERSION,
    description=DOCSTRING[0],
    long_description="\n".join(DOCSTRING[2:]),
    author=' A. Asensio Ramos',
    author_email='aasensio@iac.es',
    url='https://github.com/aasensio/pyiacsun',
    license='GNU General Public License, version 3 (GPLv3)',
    platforms='OS Independent',
    install_requires=[
        'numpy',
    ],
    ext_package='',
    ext_modules=[lib],
    classifiers=[
        'Development Status :: 4 - Beta',
        'Environment :: Console',
        'Intended Audience :: Science/Research',
        'Intended Audience :: Developers',
        "Operating System :: Unix",
        "Operating System :: POSIX",
        "Operating System :: MacOS",
        'License :: OSI Approved :: GNU General Public License v3 (GPLv3)',
        'Programming Language :: Python',
        "Programming Language :: Python :: 3.4",
        "Programming Language :: Python :: Implementation :: CPython",
        'Topic :: Scientific/Engineering',
        'Topic :: Scientific/Engineering :: Physics',
    ],
    keywords=['lte', 'radiative transfer'],
    packages=find_packages(),
    zip_safe=False,
    include_package_data=True
)


if __name__ == "__main__":
    setup(**setup_configMilne)
    setup(**setup_configLTE)

    # Attempt to remove the mod files once again.
    for filename in ["vars.mod", "atomic_functions.mod", "math_functions.mod", "math_vars.mod", "milnemod.mod", "atomicpartitionmodule.mod",
        "constants_mod.mod", "globalmodule.mod", "ltemod.mod", "synthmodule.mod", "backgroundopacitymodule.mod", "general_routines_mod.mod", 
        "hydrostaticmodule.mod", "mathsmodule.mod"]:
        try:
            os.remove(filename)
        except:
            pass