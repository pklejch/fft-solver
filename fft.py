#!/usr/bin/python3

from scipy import fft
import numpy as np
import time

print("Nacitani vstupniho vektoru...")
f = open("input.txt","r")


input=[]
i=0
for line in f:
	if i==0:
		i=i+1
		continue
	real, imag =line.split(" ")
	imag = imag.strip()
	comp = real+"+"+imag+"j"
	#print(comp)	
	input.append(complex(comp))

f.close()
x = np.array(input)
print("Zahajuji vypocet...")
start = time.clock()
output = fft(x)
end = time.clock()

diff = end - start
print("Vypocet trval %.4f s" % diff)
print("Zapis vysledneho vektoru...")
f = open("PYTHON_FFT.txt","w")
for num in output:
	f.write("(%.4f,%.4f)\n" % (float(num.real), float(num.imag)))
f.close()
