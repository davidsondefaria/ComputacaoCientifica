#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <time.h>
//#include "ranf.h"

//tamanho da rede
#define tamX 101
#define tamY 101
//raio do fio curvo
#define raioI 10
#define raioE 30
//valor dos potenciais
#define v1 -50.0
#define v2 100.0
//critério de convergência
#define erro 1e-12

double **pot;

double **alocaLin(int l){
	double **m;
	m = (double **)malloc(l*sizeof(double *));
	/*if(!codif_final){
		puts("Não há espaço para alocar memória");
	}*/
	return m;
}
void alocaCol(double **m, int l, int c){
	int i;
	for(i=0; i<l; i++){
		*(m+i) = (double *)malloc(c*sizeof(double));
		/*if(!*(codif_final+i)){
			prdoublef("Não há espaço para alocar a linha %d", i);
			exit(1);
		}*/
	}
}
void iniPot(){
	int i, j;

	for (j=0; j<tamY; j++){
		for (i=0; i<tamX; i++){
			pot[i][j] = 0.;
		}
	}
}

void fioC(){
	int i, j, a;
	double dx, dy;
	for(a=30; a<=180; a++){
		//fio interno
		dx = raioI*cos(M_PI*a/180.);
		dy = raioI*sin(M_PI*a/180.);
		i = (int)dx + tamX/2;
		j = (int)dy + tamY/2;
		pot[i][j] = v1;
		pot[i][tamY-j-1] = v1;

		//fio externo
		dx = raioE*cos(M_PI*a/180.);
		dy = raioE*sin(M_PI*a/180.);
		i = (int)dx + tamX/2;
		j = (int)dy + tamY/2;
		pot[i][j] = v2;
		pot[i][tamY-j-1] = v2;
	}
}

//métodos de relaxação
void jacobi(){
	int i, j, it=0;
	double dv=0, vn;
	do{
		dv=0;
		it++;
		for (j=1; j<=(int)(tamY/2); j++){
			for(i=1; i<tamX-1; i++){
				vn = pot[i][j];
				if (pot[i][j] != v1 && pot[i][j] != -v1 && pot[i][j] != v2 && pot[i][j] != -v2){
					pot[i][j] = (1./4.)*(pot[i-1][j]+pot[i+1][j]+pot[i][j-1]+pot[i][j+1]);
					//pot[tamX-1-i][j] = pot[i][j];
					pot[i][tamY-1-j] = pot[i][j];
					//pot[tamX-1-i][tamY-1-j] = pot[i][j];
				}
				dv += fabs(vn - pot[i][j]);
			}
		}
		//printf("%lf\n", dv);
	}while(dv>erro);
	printf("Iterações Jacobi: %d\n", it);
}
void gaussSeidel(){
	int i, j, it=0;
	double dv=0, vn;
	double **potNew;

	potNew = alocaLin(tamX);
	alocaCol(potNew, tamX, tamY);

	for (j=0; j<tamY; j++){
		for (i=0; i<tamX; i++){
			potNew[i][j] = 0.0;
		}
	}

	do{
		dv=0;
		it++;
		for (j=1; j<=(int)(tamY/2); j++){
			for(i=1; i<tamX-1; i++){
				vn = pot[i][j];
				if (pot[i][j] != v1 && pot[i][j] != -v1 && pot[i][j] != v2 && pot[i][j] != -v2){
					potNew[i][j] = (1./4.)*(potNew[i-1][j]+pot[i+1][j]+potNew[i][j-1]+pot[i][j+1]);
					//potNew[tamX-1-i][j] = potNew[i][j];
					potNew[i][tamY-1-j] = potNew[i][j];
					//potNew[tamX-1-i][tamY-1-j] = potNew[i][j];

					pot[i][j] = potNew[i][j];
					//pot[tamX-1-i][j] = pot[i][j];
					pot[i][tamY-1-j] = pot[i][j];
					//pot[tamX-1-i][tamY-1-j] = pot[i][j];

					dv += fabs(vn - potNew[i][j]);
				}
			}
		}
	//printf("%lf\n", dv);
	}while(dv>erro);
	printf("Iterações Gauss-Seidel: %d\n", it);
}
void gaussSeidelSOR(double w){
	int i, j, it=0;
	double dv=0, vn, r;
	double **potNew;

	potNew = alocaLin(tamX);
	alocaCol(potNew, tamX, tamY);

	for (j=0; j<tamY; j++){
		for (i=0; i<tamX; i++){
			potNew[i][j] = 0.0;
		}
	}

	do{
		dv=0;
		it++;
		for (j=1; j<=(int)(tamY/2); j++){
			for(i=1; i<tamX-1; i++){
				vn = pot[i][j];
				if (pot[i][j] != v1 && pot[i][j] != -v1 && pot[i][j] != v2 && pot[i][j] != -v2){
					r = (1./4.)*(potNew[i-1][j]+pot[i+1][j]+potNew[i][j-1]+pot[i][j+1]) - pot[i][j];
					potNew[i][j] = pot[i][j] + w*r;
					//potNew[tamX-1-i][j] = potNew[i][j];
					potNew[i][tamY-1-j] = potNew[i][j];
					//potNew[tamX-1-i][tamY-1-j] = potNew[i][j];

					pot[i][j] = potNew[i][j];
					//pot[tamX-1-i][j] = pot[i][j];
					pot[i][tamY-1-j] = pot[i][j];
					//pot[tamX-1-i][tamY-1-j] = pot[i][j];

					dv += fabs(vn - potNew[i][j]);
				}
			}
		}
	//printf("%lf\n", dv);
	}while(dv>erro);
	printf("Iterações Gauss-Seidel SOR: %d\n", it);
}

void main(){
	int i, j;
	char nome[15];
	double Ex, Ey;

	FILE *potencial;
	FILE *campo;

	pot = alocaLin(tamX);
	alocaCol(pot, tamX, tamY);

	printf("Potencial em Dois Fios Circulares com Abertura\n");
	iniPot();
	fioC();

	potencial = fopen("potFC-Jacobi.dat", "w");
	campo = fopen("campoFC-Jacobi.dat", "w");
	jacobi();
	//escreve pot
	for (j=0; j<tamY; j++){
		for (i=0; i<tamX; i++){
			fprintf(potencial, "%lf\t%lf\t%lf\n", i*0.1, j*0.1, pot[i][j]);
			if(i > 0 && i < tamX-1 && j > 0 && j < tamY-1){
				Ex = -(pot[i+1][j] - pot[i-1][j])/2.;
				Ey = -(pot[i][j+1] - pot[i][j-1])/2.;
				fprintf(campo, "%d\t%d\t%lf\t%lf\n", i, j, Ex, Ey);
			}
		}
		fprintf(potencial, "\n");
	}
	fclose(potencial);
	fclose(campo);

	//GaussSeidel	
	iniPot();
	fioC();
	potencial = fopen("potFC-GaussSeidel.dat", "w");
	campo = fopen("campoFC-GaussSeidel.dat", "w");
	gaussSeidel();
	//escreve pot
	for (j=0; j<tamY; j++){
		for (i=0; i<tamX; i++){
			fprintf(potencial, "%lf\t%lf\t%lf\n", i*0.1, j*0.1, pot[i][j]);
			if(i > 0 && i < tamX-1 && j > 0 && j < tamY-1){
				Ex = -(pot[i+1][j] - pot[i-1][j])/2.;
				Ey = -(pot[i][j+1] - pot[i][j-1])/2.;
				fprintf(campo, "%d\t%d\t%lf\t%lf\n", i, j, Ex, Ey);
			}
		}
		fprintf(potencial, "\n");
	}
	fclose(potencial);
	fclose(campo);

	//GaussSeidelSOR
	iniPot();
	fioC();
	potencial = fopen("potFC-GaussSeidelSOR.dat", "w");
	campo = fopen("campoFC-GaussSeidelSOr.dat", "w");
	gaussSeidelSOR(1.5);
	//escreve pot
	for (j=0; j<tamY; j++){
		for (i=0; i<tamX; i++){
			fprintf(potencial, "%lf\t%lf\t%lf\n", i*0.1, j*0.1, pot[i][j]);
			if(i > 0 && i < tamX-1 && j > 0 && j < tamY-1){
				Ex = -(pot[i+1][j] - pot[i-1][j])/2.;
				Ey = -(pot[i][j+1] - pot[i][j-1])/2.;
				fprintf(campo, "%d\t%d\t%lf\t%lf\n", i, j, Ex, Ey);
			}
		}
		fprintf(potencial, "\n");
	}
	fclose(potencial);
	fclose(campo);
}
