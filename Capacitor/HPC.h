#define NUM 103


double max;

void defPotencial(double V[NUM][NUM],double Vfut[NUM][NUM]);



void defPotencialCuboParalelo(double V[NUM][NUM],double Vfut[NUM][NUM],int tensao, int tensao2){
	int afastx,afasty,i,j;	
	
	for (i=40;i<60;i++){	
		for (j=40;j<60;j++){
			V[i][j]=tensao;
			Vfut[i][j]=tensao;
		}
	}
	

	
	for(i=25;i<75;i++) {
			j=25;
			V[i][j]=tensao2;
			Vfut[i][j]=tensao2;
	}

	for(i=25;i<75;i++) {
			j=75;
			V[i][j]=tensao2;
			Vfut[i][j]=tensao2;
	}
	
	

	
	
}

void defPotencialPlacasParalelas(double V[NUM][NUM],double Vfut[NUM][NUM],int tensao,int tensao2) {
	int i,j;
	
	for(i=1;i<NUM-1;i++) {
			j=25;
			V[i][j]=tensao;
			Vfut[i][j]=tensao;
	}

	for(i=1;i<NUM-1;i++) {
			j=75;
			V[i][j]=tensao2;
			Vfut[i][j]=tensao2;
	}

	
	

}


void defPotencialQuadrado(double V[NUM][NUM],double Vfut[NUM][NUM],int tensao){
	int afastx,afasty,i,j;	
	

	for (j=20;j<80;j++){
		V[20][j]=tensao;
		V[101-20][j]=tensao;
		Vfut[20][j]=tensao;
		Vfut[101-20][j]=tensao;
	}
	
	for (i=20;i<80;i++){	
		V[i][20]=tensao;
		V[i][101-20]=tensao;
		Vfut[i][20]=tensao;
		Vfut[i][101-20]=tensao;
	}
}





void defPotencialCubo(double V[NUM][NUM],double Vfut[NUM][NUM],int tensao){
	int afastx,afasty,i,j;	
	
	for (i=25;i<65;i++){	
		for (j=25;j<65;j++){
			V[i][j]=tensao;
			Vfut[i][j]=tensao;
		}
	}
}


void defPotencialCircular(double V[NUM][NUM],double Vfut[NUM][NUM],int tensao, int tensao2){
		
	int dx,dy,t;
	for(t=30;t<330;t++){
		dx = (int) 30*cos(t*M_PI/180);
		dy = (int) 30*sin(t*M_PI/180);
		V[50+dx][50+dy] = tensao;
		Vfut[50+dx][50+dy] = tensao;
		dx = (int) 10*cos(t*M_PI/180);
		dy = (int) 10*sin(t*M_PI/180);
		V[50+dx][50+dy] = tensao2;
		Vfut[50+dx][50+dy] = tensao2;
	}	

}


void mediaVizinhosUnicoPotencial (double V[NUM][NUM],double Vfut[NUM][NUM],int tensao)  {
	int i,j;
	max =-1;
	for(i=1;i<NUM-1;i++) {	
		for(j=1;j<NUM-1;j++) {
			if(V[i][j]!=tensao)
				Vfut[i][j] =  (V[i+1][j] + V[i-1][j] +V[i][j+1]+V[i][j-1])/4 ;
			if(max<fabs( Vfut[i][j] - V[i][j] )) 
				max = fabs( Vfut[i][j] - V[i][j] );
		}
	}
	
}


void gaussSeidelUnicoPotencial (double V[NUM][NUM],double Vfut[NUM][NUM],int tensao)  {
	int i,j;
	max =-1;
	for(i=1;i<NUM-1;i++) {	
		for(j=1;j<NUM-1;j++) {
			if(V[i][j]!=tensao)
				Vfut[i][j] =  (V[i+1][j] + Vfut[i-1][j] +V[i][j+1]+Vfut[i][j-1])/4 ;
			if(max<fabs( Vfut[i][j] - V[i][j] )) 
				max = fabs( Vfut[i][j] - V[i][j] );
		}
	}
	
}


void gaussSeidelDuploPotencial (double V[NUM][NUM],double Vfut[NUM][NUM],int tensao,int tensao2)  {
	int i,j;
	max =-1;
	for(i=1;i<NUM-1;i++) {	
		for(j=1;j<NUM-1;j++) {
			if( (V[i][j]!= tensao) && (V[i][j]!= tensao2))
				Vfut[i][j] =  (V[i+1][j] + Vfut[i-1][j] +V[i][j+1]+Vfut[i][j-1])/4 ;
			if(max<fabs( Vfut[i][j] - V[i][j] )) 
				max = fabs( Vfut[i][j] - V[i][j] );
		}
	}
	
}


void sorUnicoPotencial (double V[NUM][NUM],double Vfut[NUM][NUM],int tensao)  {
	int i,j;
	double R[NUM][NUM];
	double omega=1.8;
	max =-1;
	for(i=1;i<NUM-1;i++) {	
		for(j=1;j<NUM-1;j++) {
			if( (V[i][j]!= tensao)){
				R[i][j] =  (V[i+1][j] + Vfut[i-1][j] +V[i][j+1]+Vfut[i][j-1])/4 -V[i][j];
				Vfut[i][j]=V[i][j]+omega*R[i][j];
			}
			if(max<fabs( Vfut[i][j] - V[i][j] )) 
				max = fabs( Vfut[i][j] - V[i][j] );
		}
	}
}





void sorDuploPotencial (double V[NUM][NUM],double Vfut[NUM][NUM],int tensao,int tensao2)  {
	int i,j;
	double R[NUM][NUM];
	double omega=1.8;
	max =-1;
	for(i=1;i<NUM-1;i++) {	
		for(j=1;j<NUM-1;j++) {
			if( (V[i][j]!= tensao) && (V[i][j]!= tensao2)){
				R[i][j] =  (V[i+1][j] + Vfut[i-1][j] +V[i][j+1]+Vfut[i][j-1])/4 -V[i][j];
				Vfut[i][j]=V[i][j]+omega*R[i][j];
			}
			if(max<fabs( Vfut[i][j] - V[i][j] )) 
				max = fabs( Vfut[i][j] - V[i][j] );
		}
	}
}


void mediaVizinhosDuploPotencial (double V[NUM][NUM],double Vfut[NUM][NUM],int tensao,int tensao2)  {
	int i,j;
	max =-1;
	for(i=1;i<NUM-1;i++) {	
		for(j=1;j<NUM-1;j++) {
			if( (V[i][j]!= tensao) && (V[i][j]!= tensao2))
				Vfut[i][j] =  (V[i+1][j] + V[i-1][j] +V[i][j+1]+V[i][j-1])/4 ;
			if(max<fabs( Vfut[i][j] - V[i][j] )) 
				max = fabs( Vfut[i][j] - V[i][j] ); 
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






/*
void imprimenparceirosgnuplot(double V[NUM][NUM]){
int i,j;
FILE * gnuplotPipe = popen ("gnuplot -persistent", "w");
fprintf(gnuplotPipe,"set pm3d '-'");
fprintf(gnuplotPipe,"splot  '-'"); 
for(i=0;i<NUM;i++){
	
	for(j=0;j<NUM;j++){
	 fprintf(gnuplotPipe, "%f \t", V[i][j]);
	}
	fprintf(gnuplotPipe, "\n");
}
printf("\n");
puts("e");
//puts("pause 10");

}
*/



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

void imprime(double V[NUM][NUM]) {
	int i,j;
	FILE *pont_arq;
	pont_arq = fopen("Potencial.dat", "w");
	for (i=1;i<NUM-1;i++){
		for(j=1;j<NUM-1;j++) {
			fprintf (pont_arq,"%d %d %f  \t \n",i,j, V[i][j]);
		}
		fprintf(pont_arq,"\n");	
	}
}

void imprimeGnu(double V[NUM][NUM]) {
	int i,j;
	float deltax,deltay;
	printf("Digite o Dx e Dy, para definir a expessura da malha");
	scanf("%f",&deltax);
	scanf("%f",&deltay);
	FILE *pont_arq;
	pont_arq = fopen("Potencial.gnu", "w");
	
	fprintf(pont_arq,"set pm3d \n");
	fprintf(pont_arq,"splot '-' using 1:2:3 w p pt 0 \n");
	for (i=1;i<NUM-1;i++){
		for(j=1;j<NUM-1;j++) {
			fprintf (pont_arq,"%f %f %f  \t \n",i*deltax,j*deltay, V[i][j]);
		}
		fprintf(pont_arq,"\n");	
	}
	fprintf(pont_arq,"e\n");
	fprintf(pont_arq,"pause -1 \n");
	
}


