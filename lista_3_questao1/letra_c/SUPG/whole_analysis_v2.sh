#!/bin/bash

g++ -o mef_app mef_app.cpp
./mef_app

python3 post_proc_plots.py

# python3 post_proc_errorL2.py