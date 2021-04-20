#!/usr/bin/env python

import os
import sys
import time

#Date 2021/04
#Author Shuta Takimoto

argvs = sys.argv
start_time = time.time()


if len(argvs) != 2:
    print("\n"
          "         PROGRAM DESCRIPTION\n"
          "\n"
          "         This program is going to convert .dat files to .rho files.\n"
          "         COMMAND --> 'python Percolation.py Arg1 Arg2...'\n"
          "\n"
          "         Arg1 should be .rho file (ex. BVSxmap.rho) related to the structure you calculated.\n"
          "         Arg2 should be .dat files (ex. out.xpath). You can set a couple of files sepalately by space.\n"
          "         Like 'out.xpath out.ypath out.zpath'.\n"
          "\n"
          "         Hope your success. Bye:)\n"
          )

    sys.exit()

with open(argvs[1]) as f:
    maplines = f.readlines()

print("\n")
print("---------------------")

for datfile in argvs[2:]:
    with open(datfile) as f:
        pathlines = [0 - float(i.replace('\n', '')) for i in f.readlines()]

    loop = len(pathlines) // 5
    loopre = len(pathlines) % 5

    with open('xpath.rho', 'w') as f:
        f.write(maplines[0])
        f.write(maplines[1])
        f.write(maplines[2])
        f.write(maplines[3])
        for uu in range(loop):
            for uuu in range(5):
                try:
                    f.write('{:.8e}  '.format(pathlines[(uu + 1) * 5 - 5 + uuu]))

                except:
                    f.write('{:.8e}  '.format(0))

            f.write('\n')

        for uuu in range(loopre):
            try:
                f.write('{:.8e}  '.format(pathlines[loop * 5 + uuu]))

            except:
                f.write('{:.8e}  '.format(0))

    print("output: {}".format(datfile))

end_time = time.time()
elapsed_time = end_time - start_time
print("\n")
print("Elapsed time: {} s\n".format(int(elapsed_time)))

