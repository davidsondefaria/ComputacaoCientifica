#include <stdio.h>
#include <math.h>


#define NUM 103
#define TENSAO 20


void defPotencial(double V[NUM][NUM],double Vfut[NUM][NUM]);

void mediaVizinhos(double V[NUM][NUM],double Vfut[NUM][NUM]) {
	int i,j;
	for (i=1;i<NUM-1;i++){
		for(j=1;j<NUM-1;j++){
				if(V[i][j]!=TENSAO) {
					Vfut[i][j] =  (V[i+1][j] + V[i-1][j] +V[i][j+1]+V[i][j-1])/4;
				}
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
	
	for (i=25;i<65;i++){	
		for (j=25;j<65;j++){
			V[i][j]=TENSAO;
			Vfut[i][j]=TENSAO;
		}
	}
	
	
	
	
	
}

void imprime(double V[NUM][NUM]) {
	int i,j;
	for (i=1;i<NUM-1;i++){
		for(j=1;j<NUM-1;j++) {
			printf ("%.1f  \t", V[i][j]);
		}
		printf("\n");	
	}
}







void main(){
	int i,j;
	double V[NUM][NUM],Vfut[NUM][NUM];
	
	inicializa(V,Vfut);
	
	
	
	defPotencial(V,Vfut);

	for (i=0;i<10000;i++){
		mediaVizinhos(V,Vfut);
		copia(V,Vfut);
	}
	
	
	imprime(V);
	
		
	
		

}


