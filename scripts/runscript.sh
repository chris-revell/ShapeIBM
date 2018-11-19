#!bash
set -e
make all
./basicIBM
python3 plotter.py
open -a safari output/animated.gif
