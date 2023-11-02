#!/bin/bash

g++ -I ../eigen/ projL2_v3.cpp -o projL2
./projL2

g++ -o erroL2 errorL2_v2.cpp
./erroL2

python3 post_proc_plots.py
python3 post_proc_errorL2.py