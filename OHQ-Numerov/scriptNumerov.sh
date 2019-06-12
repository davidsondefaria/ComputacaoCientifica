#!/bin/bash
#icc NumerovDerivada.c -o numDeriv.exe
rm Energia.out
gcc NumerovDerivada.c -o numDeriv.exe -lm
for n in 0 1 2 3 41 42
do
	for tam in 1000
	do
		for xf in 10.0
		do
			./numDeriv.exe $n $tam $xf > wave$n-$tam-$xf.dat
		done
	done
done
