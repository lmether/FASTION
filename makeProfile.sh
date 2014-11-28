#!/bin/bash

gprof -l fastion gmon.out > profl.txt
gprof fastion gmon.out > prof.txt
gprof2dot prof.txt > prof.dot
dot -Tpng -oprof.png prof.dot

exit