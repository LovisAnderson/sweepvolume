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
        'graphviz==0.9',
        'logging==0.4.9.6',
        'matplotlib==2.2.3',
        'numpy==1.15.1',
        'networkx>=2.1',
        'pandas==0.23.4',
        'pplpy==0.7',
        'pygraphviz==1.5',
        'pytest==3.8.0',
        'ordered_set==3.0.1',
        'scipy==1.1.0',
        'sympy==1.3',
    ]
}

setup(**config)
