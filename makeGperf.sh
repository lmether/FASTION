#!/bin/bash

env CPUPROFILE=profg.txt CPUPROFILE_FREQUENCY=1000 ./fastion
pprof --dot ./fastion profg.txt > profg.dot
dot -Tpdf -oprofg.pdf profg.dot

exit