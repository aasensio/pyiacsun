import re
import os
import sys

try:
	from setuptools import setup
	setup
except ImportError:
	from distutils.core import setup
	setup

# Handle encoding
major, minor1, minor2, release, serial = sys.version_info
if major >= 3:
	def rd(filename):
		f = open(filename, encoding="utf-8")
		r = f.read()
		f.close()
		return r
else:
	def rd(filename):
		f = open(filename)
		r = f.read()
		f.close()
		return r
	
vre = re.compile("__version__ = \"(.*?)\"")
m = rd(os.path.join(os.path.dirname(os.path.abspath(__file__)),
                    "pyFasta", "__init__.py"))
version = vre.findall(m)[0]

setup(
    name="pyFasta",
    version=version,
    author="A. Asensio Ramos",
    author_email="aasensio@iac.es",
    packages=["pyFasta"],
    license="MIT",
    description="Python version of the FASTA algorithm",
    package_data={"": ["LICENSE"]},
    include_package_data=True,
    install_requires=["numpy","scipy"],
    classifiers=[
        "Development Status :: 5 - Production/Stable",
        "Intended Audience :: Developers",
        "Intended Audience :: Science/Research",
        "License :: OSI Approved :: MIT License",
        "Operating System :: OS Independent",
        "Programming Language :: Python",
    ],
)
