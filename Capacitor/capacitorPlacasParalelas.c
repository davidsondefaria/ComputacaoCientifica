#include <stdio.h>
#include <math.h>


#define NUM 101
#define TENSAO 20

/*
void derivada (int V[NUM+2][NUM+2])  {
	int i,j;
	double deltaX=1,deltaY=1;
	double devX[NUM][NUM],devY[NUM][NUM];
	for(i=1;i<100;i++;) {	
		for(j=1;i<100;j++){
			devX[i][j]=  (V[i+1][j] - V[i-1][j])/2*(deltaX) ;
			devY[i][j]=  (V[i][j+1] - V[i][j-1])/2*(deltaY) ;
		}
}
*/
void defBarra(double V[NUM+2][NUM+2]);

void mediaVizinhos (double V[NUM+2][NUM+2],double Vfut[NUM+2][NUM+2])  {
	int i,j;
	
	
	
	for(i=1;i<101;i++) {	
		for(j=1;j<101;j++) {
			Vfut[i][j] =  (V[i+1][j] + V[i-1][j] +V[i][j+1]+V[i][j-1])/4 ;
		}
	}

	for(i=0;i<NUM+2;i++) {	
		for(j=0;j<NUM+2;j++) {
			V[i][j] =  Vfut[i][j];
		}
	}
	defBarra(V);
	
}

void inicializa(double V[NUM+2][NUM+2]){
	int i,j;
	for(i=0;i<NUM+2;i++){
		for(j=0;j<NUM+2;j++) {
			V[i][j]=0;
		}
	}

}




void defBarra(double V[NUM+2][NUM+2]) {
	int i,j;
	
	for(i=1;i<NUM;i++) {
			j=26;
			V[i][j]=TENSAO;
	}

	for(i=1;i<NUM;i++) {
			j=75;
			V[i][j]=-TENSAO;
	}

	for(i=1;i<NUM;i++) {
			j=26;
			V[i][j]=TENSAO;
	}

	for(i=1;i<NUM;i++) {
			j=75;
			V[i][j]=-TENSAO;
	}


}

void imprime(double V[NUM+2][NUM+2]) {
	int i,j;
	for (i=1;i<NUM;i++){
		for(j=1;j<NUM;j++) {
			printf ("%f \t", V[i][j]);
			
		}
		printf("\n");	
	}
}



void main() {


double V[NUM+2][NUM+2],Vfut[NUM+2][NUM+2];
int i;

inicializa(V);
defBarra(V);

for (i=0;i<100;i++){
	mediaVizinhos(V,Vfut);
}

imprime(V);


















}


