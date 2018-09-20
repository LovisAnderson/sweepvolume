#!/usr/bin/env bash
virtualenv venv
source venv/bin/activate
pip install gmpy2==2.1.0a1 --no-binary ":all:"
pip install -r requirements.txt
python setup.py install