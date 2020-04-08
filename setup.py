import os
#from distutils.core import setup, Extension
from setuptools import setup, Extension

os.environ["CC"] = "gcc"
#compile_args = ["-g -Wall -O2"]
link_args    = ["-lm"]

ext = Extension('phanotate_connect',
                language='gcc',
#               extra_compile_args=compile_args,
                extra_link_args=link_args,
                include_dirs=[
                             '.',
                             '...',
                             os.path.join(os.getcwd(), 'include'),
                ],
                library_dirs = [os.getcwd(),],
                sources = ['src/phanotate_connect.c'])

with open("README.md", "r") as fh:
    long_desc = fh.read()


setup (
    name = 'phanotate_connect',
    version = '0.1',
    author = "Katelyn McNair",
    author_email = "deprekate@gmail.com",
    description = 'A package for connecting nearby orfs together in the directed graph',
    long_description = long_desc,
    long_description_content_type="text/markdown",
    url =  "https://github.com/deprekate/phanotate",
    classifiers=[
        "Programming Language :: Python :: 3",
        "License :: OSI Approved :: GNU General Public License v3 or later (GPLv3+)",
        "Operating System :: OS Independent",
    ],
    python_requires='>3.5.2',
    ext_modules = [ext]
)
