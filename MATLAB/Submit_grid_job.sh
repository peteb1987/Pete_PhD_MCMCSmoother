#!/bin/bash
qsub -l lr=1 -m be -t 1-300 ge_script.sh