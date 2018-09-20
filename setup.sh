#!/usr/bin/env bash
virtualenv venv
source venv/bin/activate
pip install gmpy2==2.1.0a1 --no-binary ":all:"
python setup.py install