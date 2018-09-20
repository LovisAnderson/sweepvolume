from setuptools import setup, find_packages

config = {
    'name': 'sweepvolume',
    'description': 'Package to calculate sweep-plane volume function for unions of polyhedra.',
    'author': 'Lovis Anderson',
    'author_email': 'lovisanderson@gmail.com',
    'version': '0.1',
    'license': 'GPL',
    'packages': find_packages(),
    'install_requires': [
        'graphviz',
        'logging',
        'matplotlib',
        'numpy',
        'networkx',
        'pandas',
        'pplpy',
        'pygraphviz',
        'pytest',
        'ordered_set',
        'scipy',
        'sympy',
    ]
}

setup(**config)