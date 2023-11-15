#!/bin/bash

g++ -I ../eigen/ mef_app.cpp -o mef_app
./mef_app

g++ -o erroL2 errorL2_v2.cpp
./erroL2

python3 post_proc_plots.py
python3 post_proc_errorL2.py