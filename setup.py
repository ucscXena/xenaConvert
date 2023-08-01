from setuptools import setup

setup(name='xenaConvert',
      version='0.1',
      description='Python Utilities',
      author='Jing Zhu',
      author_email='jzhu@soe.ucsc.edu',
      url='https://github.com/ucscXena/xenaConvert',
      install_requires=[
        'pandas',
        'scanpy',
        'umap-learn',
        'leidenalg',
        'louvain'
      ],
     )
