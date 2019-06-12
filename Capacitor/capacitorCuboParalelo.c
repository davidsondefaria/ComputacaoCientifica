#include <stdio.h>
#include <math.h>


#define NUM 103
#define TENSAO 20
#define TENSAO2 -20


void defPotencial(double V[NUM][NUM],double Vfut[NUM][NUM]);

void mediaVizinhos(double V[NUM][NUM],double Vfut[NUM][NUM]) {
	int i,j;
	for (i=1;i<NUM-1;i++){
		for(j=1;j<NUM-1;j++){
				if( (V[i][j]!=TENSAO) &&  V[i][j] != TENSAO2) {
					Vfut[i][j] =  (V[i+1][j] + V[i-1][j] +V[i][j+1]+V[i][j-1])/4;
				}
		}
		
		
	}
	
}

void campoEletrico (double V[NUM][NUM]){
	int i,j;
	for (i=1;i<NUM-1;i++){
		for (j=1;j<NUM-1;j++){
			printf("%d \t %d \t %f \t %f \n ", i,j, -(V[i+1][j]-V[i-1][j])/2,-(V[i][j+1]-V[i][j-1])/2  );
		}	
		printf("\n");
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
	
	for (i=40;i<60;i++){	
		for (j=40;j<60;j++){
			V[i][j]=TENSAO;
			Vfut[i][j]=TENSAO;
		}
	}
	

	
	for(i=25;i<75;i++) {
			j=25;
			V[i][j]=TENSAO2;
			Vfut[i][j]=TENSAO2;
	}

	for(i=25;i<75;i++) {
			j=75;
			V[i][j]=TENSAO2;
			Vfut[i][j]=TENSAO2;
	}
	
	

	
	
}

void imprime(double V[NUM][NUM]) {
	int i,j;
	for (i=1;i<NUM-1;i++){
		for(j=1;j<NUM-1;j++) {
			printf ("%f  \t", V[i][j]);
		}
		printf("\n");	
	}
}







void main(){
	int i,j;
	double V[NUM][NUM],Vfut[NUM][NUM];
	double Ex[NUM][NUM];
	double Ey[NUM][NUM];
	
	inicializa(V,Vfut);

	
	defPotencial(V,Vfut);

	for (i=0;i<10000;i++){
		mediaVizinhos(V,Vfut);
		copia(V,Vfut);
	}
	
	
	imprime(V);
	
		
	//campoEletrico(V);
	
		

}


