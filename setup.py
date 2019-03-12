from setuptools import setup

config = {
    'name': 'sweepvolume',
    'description': 'Package to calculate sweep-plane volume function for unions of polyhedra.',
    'author': 'Lovis Anderson',
    'author_email': 'lovisanderson@gmail.com',
    'version': '0.1',
    'license': 'GPL',
    'packages': ['sweepvolume'],
    'package_dir': {'sweepvolume': 'sweepvolume'},
    'install_requires': [
        'graphviz',
        'matplotlib',
        'numpy',
        'networkx',
        'pplpy>=0.7',
        'pytest',
        'scipy',
        'sympy',
    ]
}

setup(**config)
