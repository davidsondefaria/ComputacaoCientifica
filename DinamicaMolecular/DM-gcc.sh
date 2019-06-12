#!/bin/bash
echo "#################################################################################" >> profile.txt
echo "Início da simulação" >> profile.txt
echo " " >> profile.txt
gcc DM.c -pg -o dm.exe -O0
./dm.exe >> profile.txt
gprof dm.exe >> profile.txt
echo " " >> profile.txt
Rscript script_Gr.R
echo "Fim da simulação" >> profile.txt
echo "#################################################################################" >> profile.txt
echo " " >> profile.txt
mkdir Resultado
rm arqaux.dat
mv *.dat Resultado
mv *.xsf Resultado
rm *.out