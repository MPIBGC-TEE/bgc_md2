#!/bin/bash
python -m pip install --upgrade pip
pip install  -r requirements.non_src -r requirements.only_pip -r requirements.test 
pip install -r requirements.github

