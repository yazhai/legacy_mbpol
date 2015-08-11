#!/bin/bash

set -o nounset

../optimize $1.xyz -o $1_opt.xyz > $1_opt.log
../normal-modes $1_opt.xyz > $1_frequencies.log
