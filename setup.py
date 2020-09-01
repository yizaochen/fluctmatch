from setuptools import setup, find_packages

setup(name='fluctmatch', 
      version='0.1',
      packages=find_packages(),
      url='https://github.com/yizaochen/fluctmatch.git',
      author='Yizao Chen',
      author_email='yizaochen@gmail.com',
      install_requires=[
          'numpy',
          'pandas',
          'scipy',
          'MDAnalysis'
      ]
      )