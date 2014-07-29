#!/usr/bin/env python

import sys
from pylab import *

print(" ========================================================== ")
print(" Usage: "+sys.argv[0]+" some.ts "+"[some.svg ...]")
print(" ========================================================== ")

if len(sys.argv) == 1:
    print(" ERROR: Need to specify the .ts file!")
    exit()

tsfile = sys.argv[1]

step = []
coor = []

for line in open(tsfile):
    items = line.split()
    if len(items) > 0 and items[0] == "#all":
        step.append(int(items[1]))
        coor.append(float(items[4]))

plot(step, coor, color="green", linewidth=2.5, linestyle="-")
xlabel("step")
ylabel(r'$\AA$')

if len(sys.argv) == 3:
    savefig(sys.argv[2])
else:
    show()
