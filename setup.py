'''------------------------ Setup Script for Owlpy--------------------------'''
from __future__ import print_function
from setuptools import setup
from setuptools.command.test import test as TestCommand
import codecs
import os
import sys
import re

HERE = os.path.abspath(os.path.dirname(__file__))


def read(*parts):
    """Return multiple read calls to different readable objects as a single
    string."""
    # intentionally *not* adding an encoding option to open
    return codecs.open(os.path.join(HERE, *parts), 'r').read()

LONG_DESCRIPTION = read('README.rst')


class PyTest(TestCommand):
    def finalize_options(self):
        TestCommand.finalize_options(self)
        self.test_args = [
            '--strict',
            '--verbose',
            '--tb=long',
            'tests']
        self.test_suite = True

    def run_tests(self):
        import pytest
        errno = pytest.main(self.test_args)
        sys.exit(errno)

setup(
    name='owlpy',
    version='1.0.0',
    url='http://github.com/jbeleno/owlpy/',
    license='MIT License',
    author='Juan Sebastián Beleño Díaz',
    tests_require=['pytest', 'pytest-cov'],
    install_requires=[
        'numpy>=1.11.1'
        ],
    cmdclass={'test': PyTest},
    author_email='jsbeleno@gmail.com',
    description='An open source time series library for Python using Matrix Profile',
    long_description=LONG_DESCRIPTION,
    packages=['owlpy'],
    include_package_data=True,
    platforms='any',
    test_suite='tests.test_owlpy',
    classifiers=[
        'Programming Language :: Python',
        'Development Status :: 1 - Beta',
        'Natural Language :: English',
        'Environment :: Desktop Environment',
        'Intended Audience :: Data Scientists',
        'License :: MIT License',
        'Operating System :: OS Independent',
        'Topic :: Software Development :: Libraries :: Python Modules',
        'Topic :: Software Development :: Libraries :: Application Frameworks'
        ],
)