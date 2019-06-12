#include <stdio.h>
#include <math.h>


#define NUM 103
#define TENSAO 100


void defPotencial(double V[NUM][NUM]);

void mediaVizinhos(double V[NUM][NUM],double Vfut[NUM][NUM]) {
	int i,j;
	for (i=1;i<NUM-1;i++){
		for(j=1;j<NUM-1;j++){
			Vfut[i][j] =  (V[i+1][j] + V[i-1][j] +V[i][j+1]+V[i][j-1])/4;	
		}
	}
	defPotencial(V);
}


void derivada(double V[NUM][NUM],double Vfut[NUM][NUM]){
	int i,j;

	for (i=2;i<NUM-2;i++){
		for(j=2;j<NUM-2;j++){
			Vfut[i][j] =  (V[i+1][j]+V[i-1][j])/2;	
		}
	}


}


void inicializa(double V[NUM][NUM]){
	int i,j;
	for(i=0;i<NUM;i++){
		for(j=0;j<NUM;j++) {
			V[i][j]=0;
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

void defPotencial(double V[NUM][NUM]){
		
	int dx,dy,t;
	for(t=30;t<330;t++){
		dx = (int) 30*cos(t*M_PI/180);
		dy = (int) 30*sin(t*M_PI/180);
		V[50+dx][50+dy] = 50;
		dx = (int) 10*cos(t*M_PI/180);
		dy = (int) 10*sin(t*M_PI/180);
		V[50+dx][50+dy] = -100;
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
	
	inicializa(V);
	defPotencial(V);

	
	for(i=0;i<100;i++){
		
		mediaVizinhos(V,Vfut);
		copia(V,Vfut);
	}
	//derivada(V,Vfut);
	imprime(V);
	
		
	
		

}


