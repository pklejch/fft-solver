#!/usr/bin/python3

import random
import sys


if len(sys.argv) < 2:
	length=64
else:
	length = int(sys.argv[1])

f= open("input.txt","w")


f.write(str(length)+'\n')
for i in range(0,length):
	num = random.random()
	real = random.randint(1,10)
	imag = random.randint(1,10)
	real = real * num
	imag = imag * num
	f.write(str(real)+" "+str(imag)+'\n')

f.close()
