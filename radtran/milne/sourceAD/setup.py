from distutils.core import setup
from distutils.extension import Extension
from Cython.Distutils import build_ext
# This line only needed if building with NumPy in Cython file.
from numpy import get_include
from os import system

# compile the fortran modules without linking
fortran_mod_comp = 'make clean ; make'
print fortran_mod_comp
system(fortran_mod_comp)

ext_modules = [Extension(# module name:
                         "pymilneAD",
                         # source file:
                         ['pymilne.pyx'],
                         # other compile args for gcc
                         extra_compile_args=['-fPIC', '-O3'],
                         # other files to link to
                         extra_link_args=['vars.o', 'maths.o', 'atomic.o', 'milne.o', 'ad.o', '-lgfortran'])]

setup(name = 'pymilneAD',
      cmdclass = {'build_ext': build_ext},
      # Needed if building with NumPy.
      # This includes the NumPy headers when compiling.
      include_dirs = [get_include()],
      ext_modules = ext_modules,
      version='0.1',
      py_modules=['milne'],
      description='Milne-Eddington synthesis package',
      author='Andres Asensio Ramos',
      author_email='aasensio@iac.es',      
)

system('cp pymilneAD.so ../')
