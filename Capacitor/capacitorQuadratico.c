#include <stdio.h>
#include <math.h>


#define NUM 103
#define TENSAO 1000


void defPotencial(double V[NUM][NUM],double Vfut[NUM][NUM]);

void mediaVizinhos(double V[NUM][NUM],double Vfut[NUM][NUM]) {
	int i,j;
	for (i=1;i<NUM-1;i++){
		for(j=1;j<NUM-1;j++){
			if( V[i][j] != TENSAO)
				Vfut[i][j] =  (V[i+1][j] + V[i-1][j] +V[i][j+1]+V[i][j-1])/4;	
		}
	}
}


void inicializa(double V[NUM][NUM],double Vfut[NUM][NUM]){
	int i,j;
	for(i=0;i<NUM;i++){
		for(j=0;j<NUM;j++) {
			V[i][j]=0;
			Vfut[i][j]=0;
		}
	}
}

void copia(double V[NUM][NUM],double Vfut[NUM][NUM]){
	int i,j;
	for(i=0;i<NUM;i++){
		for(j=0;j<NUM;j++) {
			V[i][j]=Vfut[i][j];
		}
	}
}

void defPotencial(double V[NUM][NUM],double Vfut[NUM][NUM]){
	int afastx,afasty,i,j;	
	

	for (j=20;j<80;j++){
		V[20][j]=TENSAO;
		V[101-20][j]=TENSAO;
		Vfut[20][j]=TENSAO;
		Vfut[101-20][j]=TENSAO;
	}
	
	for (i=20;i<80;i++){
	
		V[i][20]=TENSAO;
		V[i][101-20]=TENSAO;
		Vfut[i][20]=TENSAO;
		Vfut[i][101-20]=TENSAO;
	}

	
	
	
	
	
	
	
	
	
}

void imprime(double V[NUM][NUM]) {
	int i,j;
	for (i=1;i<NUM-1;i++){
		for(j=1;j<NUM-1;j++) {
			printf ("%f \t", V[i][j]);
		}
		printf("\n");	
	}
}






void main(){
	int i,j;
	double V[NUM][NUM],Vfut[NUM][NUM];
	
	inicializa(V,Vfut);
	
	defPotencial(V,Vfut);

	for (i=0;i<1000;i++){
		mediaVizinhos(V,Vfut);
		copia(V,Vfut);
	}
	imprime(V);
}


