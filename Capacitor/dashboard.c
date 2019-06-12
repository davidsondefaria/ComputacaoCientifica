#include <stdio.h>
#include <math.h>
#include "HPC.h"



void main(){
	int i,j,x=-1;
	double V[NUM][NUM],Vfut[NUM][NUM];
	double Ex[NUM][NUM];
	double Ey[NUM][NUM];
	int tensao, tensao2,passos=0,metodo;
	inicializa(V,Vfut);
	
	

	printf("Selecione a opção desejada: \n");
	printf("Digite: 1 para Capacitor de placas paralelas\n   Digite 2 para Capacitor barra quadrada\n Digite 3 para Capacitor fio quadrado\n Digite 4 para capacitor fio circular\n Digite 5 para capacitor com cubo e placas paralelas\n ");

		
		scanf("%d", &x);
		
		if (x==1){ 
			printf("Digite a tensao na primeira placa 1 \n");
			scanf("%d", &tensao);
			printf("Digite a tensao na primeira placa 2 \n");
			scanf("%d", &tensao2);
			printf("Digite 0 para utilizar o método de Jacobi, 1 para Gauss-Seidel ou 2 para SOR \n");
			scanf("%d", &metodo);

			defPotencialPlacasParalelas(V,Vfut,tensao,tensao2);
			if(metodo==0){
				do{		
						//gaussSeidelDuploPotencial(V,Vfut,tensao,tensao2);
						mediaVizinhosDuploPotencial(V,Vfut,tensao,tensao2);
						copia(V,Vfut);
						passos++;
				
				}while(max >0.000001);
			
			} else if(metodo==1){
				do{		
						gaussSeidelDuploPotencial(V,Vfut,tensao,tensao2);
						//mediaVizinhosDuploPotencial(V,Vfut,tensao,tensao2);
						copia(V,Vfut);
						passos++;
				
				}while(max >0.000001);

			} else if(metodo==2){
				do{		
						sorDuploPotencial(V,Vfut,tensao,tensao2);
						//mediaVizinhosDuploPotencial(V,Vfut,tensao,tensao2);
						copia(V,Vfut);
						passos++;
				
				}while(max >0.000001);
			}


			imprimeGnu(V);
		}else if (x==2){
			printf("Digite a tensao na barra quadrada \n");
			scanf("%d", &tensao);
			printf("Digite 0 para utilizar o método de Jacobi, 1 para Gauss-Seidel ou 2 para SOR  \n");
			scanf("%d", &metodo);

			defPotencialCubo(V,Vfut,tensao);
			if(metodo==0){
				do{	
						//gaussSeidelUnicoPotencial(V,Vfut,tensao);
						mediaVizinhosUnicoPotencial(V,Vfut,tensao);
						copia(V,Vfut);
						passos++;
				}while(max >0.00000001);
			} else if(metodo==1){
				do{	
						gaussSeidelUnicoPotencial(V,Vfut,tensao);
						//mediaVizinhosUnicoPotencial(V,Vfut,tensao);
						copia(V,Vfut);
						passos++;
				}while(max >0.00000001);
		
			}else if(metodo==2){
				do{		
						sorUnicoPotencial(V,Vfut,tensao);
						//mediaVizinhosDuploPotencial(V,Vfut,tensao,tensao2);
						copia(V,Vfut);
						passos++;
				
				}while(max >0.000001);
			}
	
		
			imprimeGnu(V);
		}else if (x==3){
			printf("Digite a tensao na barra quadrada \n");
			scanf("%d", &tensao);
			printf("Digite 0 para utilizar o método de Jacobi, 1 para Gauss-Seidel ou 2 para SOR  \n");
			scanf("%d", &metodo);
			defPotencialQuadrado(V,Vfut,tensao);
			if(metodo==0){
				do{	
						//gaussSeidelUnicoPotencial(V,Vfut,tensao);
						mediaVizinhosUnicoPotencial(V,Vfut,tensao);
						copia(V,Vfut);
						passos++;
				}while(max >0.000001);
			} else if(metodo==1){
				do{	
						gaussSeidelUnicoPotencial(V,Vfut,tensao);
						//mediaVizinhosUnicoPotencial(V,Vfut,tensao);
						copia(V,Vfut);
						passos++;
				}while(max >0.000001);

			}else if(metodo==2){
				do{		
						sorUnicoPotencial(V,Vfut,tensao);
						//mediaVizinhosDuploPotencial(V,Vfut,tensao,tensao2);
						copia(V,Vfut);
						passos++;
				
				}while(max >0.000001);
			}


			imprimeGnu(V);


			
		}else if (x==4){
			printf("Digite a tensão no círculo externo \n");
			scanf("%d", &tensao);
			printf("Digite a tensão no círculo interno \n");
			scanf("%d", &tensao2);
			printf("Digite 0 para utilizar o método de Jacobi, 1 para Gauss-Seidel ou 2 para SOR  \n");
			scanf("%d", &metodo);
			
			defPotencialCircular(V,Vfut,tensao,tensao2);
			if(metodo==0){

				do{
						//gaussSeidelDuploPotencial(V,Vfut,tensao,tensao2);
						mediaVizinhosDuploPotencial(V,Vfut,tensao,tensao2);
						copia(V,Vfut);
						passos++;
				}while(max >0.000001);
			}else if(metodo==1){
				do{
						gaussSeidelDuploPotencial(V,Vfut,tensao,tensao2);
						//mediaVizinhosDuploPotencial(V,Vfut,tensao,tensao2);
						copia(V,Vfut);
						passos++;
				}while(max >0.000001);

				
			} else if(metodo==2){
				do{		
						sorDuploPotencial(V,Vfut,tensao,tensao2);
						//mediaVizinhosDuploPotencial(V,Vfut,tensao,tensao2);
						copia(V,Vfut);
						passos++;
				
				}while(max >0.000001);
			}


			
			imprimeGnu(V);
		}else if (x==5){
			printf("Digite a tensão na barra quadrada \n");
			scanf("%d", &tensao);
			printf("Digite a tensão nas placas \n");
			scanf("%d", &tensao2);
			printf("Digite 0 para utilizar o método de Jacobi, 1 para Gauss-Seidel ou 2 para SOR \n");
			scanf("%d", &metodo);


			defPotencialCuboParalelo(V,Vfut,tensao,tensao2);
			if(metodo==0){
				do{
						//gaussSeidelDuploPotencial(V,Vfut,tensao,tensao2);
						mediaVizinhosDuploPotencial(V,Vfut,tensao,tensao2);
						copia(V,Vfut);
						passos++;
				}while(max >0.000001);
			}else if(metodo==1){
				do{
						gaussSeidelDuploPotencial(V,Vfut,tensao,tensao2);
						//mediaVizinhosDuploPotencial(V,Vfut,tensao,tensao2);
						copia(V,Vfut);
						passos++;
				}while(max >0.000001);

			}else if(metodo==2){
				do{		
						sorDuploPotencial(V,Vfut,tensao,tensao2);
						//mediaVizinhosDuploPotencial(V,Vfut,tensao,tensao2);
						copia(V,Vfut);
						passos++;
				
				}while(max >0.000001);
			}

			imprimeGnu(V);
		}

		
	
	printf("Passos %d", passos);
	
	
		
	//campoEletrico(V);
	
		

}


