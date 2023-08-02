from setuptools import setup

setup(name='xenaConvert',
      version='0.1',
      description='Python Utilities',
      url='https://github.com/ucscXena/xenaConvert',
      py_modules = ['convert'],
      install_requires=[
        'pandas',
        'scanpy',
        'umap-learn',
        'leidenalg',
        'louvain'
      ],
     )
