#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <time.h>
//#include "ranf.h"

//tamanho da rede
#define tamX 101
#define tamY 101
//localização das barras/fios/caixa
#define potX1 25
#define potY1 10
#define quadX1 40
#define quadY1 40
//raio do fio curvo
#define raioI 10
#define raioE 30
//valor dos potenciais
#define v1 50.0
#define v2 100.0
//critério de convergência
#define erro 1e-6

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
//inicialização
void placasParalelas(){
	int i, j;
	for(j=potY1; j<=(int)(tamY/2.); j++){
		pot[potX1][j] = v1;
		pot[potX1][tamY-1-j] = v1;
		pot[tamX-1-potX1][j] = -v1;
		pot[tamX-1-potX1][tamY-1-j] = -v1;
	}
}
void quadradoCondutor(){
	int i,j;
	for(j=quadY1; j<=(int)(tamY/2.); j++){
		for(i=quadX1; i<=(int)(tamX/2.); i++){
			pot[i][j] = v1;
			pot[tamX-1-i][j] = v1;
			pot[i][tamY-1-j] = v1;
			pot[tamX-1-i][tamY-1-j] = v1;
		}
	}
}
void quadPlacas(){
	int i,j;
	quadradoCondutor();
	for(j=potY1; j<=(int)(tamY/2.); j++){
		pot[potX1][j] = -v1;
		pot[potX1][tamY-1-j] = -v1;
		pot[tamX-1-potX1][j] = -v1;
		pot[tamX-1-potX1][tamY-1-j] = -v1;
	}
}
void fioQuad(){
	int i, j;
	for(j=potY1; j<=(int)(tamY/2.); j++){
		pot[potX1][j] = v1;
		pot[potX1][tamY-1-j] = v1;
		pot[tamX-1-potX1][j] = v1;
		pot[tamX-1-potX1][tamY-1-j] = v1;
	}
	for(i=potX1; i<=(int)(tamY/2.); i++){
		pot[i][potY1] = v1;
		pot[tamX-1-i][potY1] = v1;
		pot[i][tamY-1-potY1] = v1;
		pot[tamX-1-i][tamY-1-potY1] = v1;
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
			for(i=1; i<=(int)(tamX/2); i++){
				vn = pot[i][j];
				if (pot[i][j] != v1 && pot[i][j] != -v1 && pot[i][j] != v2 && pot[i][j] != -v2){
					pot[i][j] = (1./4.)*(pot[i-1][j]+pot[i+1][j]+pot[i][j-1]+pot[i][j+1]);
					pot[tamX-1-i][j] = -pot[i][j];
					pot[i][tamY-1-j] = pot[i][j];
					pot[tamX-1-i][tamY-1-j] = -pot[i][j];
				}
				dv += abs(vn - pot[i][j]);
			}
		}
		printf("%lf\n", dv);
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
			for(i=1; i<=(int)(tamX/2); i++){
				vn = potNew[i][j];
				if (abs(pot[i][j]) == v1 && abs(pot[i][j]) != v2){
					potNew[i][j] = pot[i][j];
					potNew[tamX-1-i][j] = -potNew[i][j];
					potNew[i][tamY-1-j] = potNew[i][j];
					potNew[tamX-1-i][tamY-1-j] = -potNew[i][j];
				}
				else{
					potNew[i][j] = (1./4.)*(potNew[i-1][j]+pot[i+1][j]+potNew[i][j-1]+pot[i][j+1]);
					potNew[tamX-1-i][j] = -potNew[i][j];
					potNew[i][tamY-1-j] = potNew[i][j];
					potNew[tamX-1-i][tamY-1-j] = -potNew[i][j];

					pot[i][j] = potNew[i][j];
					pot[tamX-1-i][j] = -pot[i][j];
					pot[i][tamY-1-j] = pot[i][j];
					pot[tamX-1-i][tamY-1-j] = -pot[i][j];

					dv += abs(vn - potNew[i][j]);
				}
			}
		}
	printf("%lf\n", dv);
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
			for(i=1; i<=(int)(tamX/2); i++){
				vn = potNew[i][j];
				if (abs(pot[i][j]) == v1 && abs(pot[i][j]) != v2){
					potNew[i][j] = pot[i][j];
					potNew[tamX-1-i][j] = -potNew[i][j];
					potNew[i][tamY-1-j] = potNew[i][j];
					potNew[tamX-1-i][tamY-1-j] = -potNew[i][j];
				}
				else{
					r = (1./4.)*(potNew[i-1][j]+pot[i+1][j]+potNew[i][j-1]+pot[i][j+1]) - pot[i][j];
					potNew[i][j] = pot[i][j] + w*r;
					potNew[tamX-1-i][j] = -potNew[i][j];
					potNew[i][tamY-1-j] = potNew[i][j];
					potNew[tamX-1-i][tamY-1-j] = -potNew[i][j];

					pot[i][j] = potNew[i][j];
					pot[tamX-1-i][j] = -pot[i][j];
					pot[i][tamY-1-j] = pot[i][j];
					pot[tamX-1-i][tamY-1-j] = -pot[i][j];

					dv += abs(vn - potNew[i][j]);
				}
			}
		}
	printf("%lf\n", dv);
	}while(dv>erro);
	printf("Iterações Gauss-Seidel SOR: %d\n", it);
}



void main(){
	int i, j;
	char nome[15];

	FILE *potencial;
	FILE *campo;

	pot = alocaLin(tamX);
	alocaCol(pot, tamX, tamY);

//-------------------------------------------------------//
	printf("Potencial Placas Paralelas\n");

	//Jacobi	
	iniPot();
	placasParalelas();
	potencial = fopen("potPP-Jacobi.dat", "w");
	campo = fopen("campoPP-Jacobi.dat", "w");
	jacobi();
	//escreve pot
	for (j=0; j<tamY; j++){
		for (i=0; i<tamX; i++){
			fprintf(potencial, "%.0lf\t", pot[i][j]);
		}
		fprintf(potencial, "\n");
	}
	fclose(potencial);
	fclose(campo);

	//GaussSeidel	
	iniPot();
	placasParalelas();
	potencial = fopen("potPP-GaussSeidel.dat", "w");
	campo = fopen("campoPP-GaussSeidel.dat", "w");
	gaussSeidel();
	//escreve pot
	for (j=0; j<tamY; j++){
		for (i=0; i<tamX; i++){
			fprintf(potencial, "%.0lf\t", pot[i][j]);
		}
		fprintf(potencial, "\n");
	}
	fclose(potencial);
	fclose(campo);

	//GaussSeidelSOR
	iniPot();
	placasParalelas();
	potencial = fopen("potPP-GaussSeidelSOR.dat", "w");
	campo = fopen("campoPP-GaussSeidelSOr.dat", "w");
	gaussSeidelSOR(1.9);
	//escreve pot
	for (j=0; j<tamY; j++){
		for (i=0; i<tamX; i++){
			fprintf(potencial, "%.0lf\t", pot[i][j]);
		}
		fprintf(potencial, "\n");
	}
	fclose(potencial);
	fclose(campo);

/*
//-------------------------------------------------------//
	printf("Potencial Condutor Quadrado\n");
	iniPot();
	quadradoCondutor();

	potencial = fopen("potCQ-Jacobi.dat", "w");
	campo = fopen("campoCQ-Jacobi.dat", "w");
	jacobi();
	//escreve pot
	for (j=0; j<tamY; j++){
		for (i=0; i<tamX; i++){
			fprintf(potencial, "%.0lf\t", pot[i][j]);
		}
		fprintf(potencial, "\n");
	}
	fclose(potencial);
	fclose(campo);

//-------------------------------------------------------//
	printf("Potencial Condutor Quadrado com Placas Paralelas\n");
	iniPot();
	quadPlacas();

	potencial = fopen("potQP-Jacobi.dat", "w");
	campo = fopen("campoQP-Jacobi.dat", "w");
	jacobi();
	//escreve pot
	for (j=0; j<tamY; j++){
		for (i=0; i<tamX; i++){
			fprintf(potencial, "%.0lf\t", pot[i][j]);
		}
		fprintf(potencial, "\n");
	}
	fclose(potencial);
	fclose(campo);

//-------------------------------------------------------//
	printf("Potencial num Fio Quadrado\n");
	iniPot();
	fioQuad();

	potencial = fopen("potQP-Jacobi.dat", "w");
	campo = fopen("campoQP-Jacobi.dat", "w");
	jacobi();
	//escreve pot
	for (j=0; j<tamY; j++){
		for (i=0; i<tamX; i++){
			fprintf(potencial, "%.0lf\t", pot[i][j]);
		}
		fprintf(potencial, "\n");
	}
	fclose(potencial);
	fclose(campo);

//-------------------------------------------------------//
	printf("Potencial em Dois Fios Circulares com Abertura\n");
	iniPot();
	fioC();

	potencial = fopen("potFC-Jacobi.dat", "w");
	campo = fopen("campoFC-Jacobi.dat", "w");
	jacobi();
	//escreve pot
	for (j=0; j<tamY; j++){
		for (i=0; i<tamX; i++){
			fprintf(potencial, "%.0lf\t", pot[i][j]);
		}
		fprintf(potencial, "\n");
	}
	fclose(potencial);
	fclose(campo);
*/

}