import io
from os.path import dirname, join
from setuptools import setup

def get_version(relpath):
  '''Read version info from a file without importing it'''
  for line in io.open(join(dirname(__file__), relpath), encoding='cp437'):
    if '__version__' in line:
      if '"' in line:
        # __version__ = "0.9"
        return line.split('"')[1]
      elif "'" in line:
        return line.split("'")[1]


setup(name='mentos',
      version=get_version('mentos/__init__.py'),
      description='MENTOS: Maximum entropy approach for estimating metabolite concentrations and fluxes',
      url='http://github.com/djinnome/mentos',
      author='Jeremy Zucker',
      author_email='djinnome@gmail.com',
      install_requires=['pandas','scipy','cvxpy','numpy','cobra','python-libsbml','pyOpt>=1.2.0.1','kitchen'],
      dependency_links=[
        "https://github.com/djinnome/pyopt/tarball/master#egg=pyOpt-1.2.0.1"],
      license='MIT',
      packages=['mentos'],
      zip_safe=False)
