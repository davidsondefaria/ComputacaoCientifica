#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <time.h>
#include "ranf.h"
#include "estrutura.h"
#include "propriedades.h"
#include "forces.h"
#include "integrationv1.h"
//#include "forcaepotencia.h"
//#include "ingrecao.h"
 
void metIntegracao();
void metIniVelocidade();
double metPotencial(double);

void main(){

	//Inicialização de Variáveis iteração, loop, arquivos e condições
	double t=0, t_tpt, amostragem;
	double auxxx;
	int i, ntemp, j=0;	
	int contador=0;
	int contFilme=1, npassos=500;	//500 passos temporais para o filme
	char nomeAux[100];
	FILE *energias_temp;
	FILE *filme;
	FILE *teste;
	teste = fopen("teste.dat", "w");

	//le dados de entrada
	leInput();
	amostragem = 0.1*nt;
	printf("%g\n", amostragem);
	ntemp = (nt/dt)/npassos;

	//remove o arquivo de displacement
	snprintf(nomeAux, 100, "disp-T%.0lf-I%d-V%d-P%d.dat", temperatura_ideal, inpIntegracao, inpIniVel, inpPot);
	remove(nomeAux);

	//cria o arquivo de energia
	snprintf(nomeAux, 100, "energia-T%.0lf-I%d-V%d-P%d.dat", temperatura_ideal, inpIntegracao, inpIniVel, inpPot);
	energias_temp=fopen(nomeAux,"w");

	//cria o arquivo do filme
	snprintf(nomeAux, 100, "filme-T%.0lf-I%d-V%d-P%d.xsf", temperatura_ideal, inpIntegracao, inpIniVel, inpPot);
	filme = fopen(nomeAux,"w");

	//Inicialização da Rede, das posições, velocidade, aceleração
	redeFCC();
	metIniVelocidade();//iniVelTeoPart();
	iniAc0();
	copia_particula();
	primeira_integracao();

	//Inicialização do parâmetros de G(r), mean square e Normalização da temperatura
	Gr_init();
	diff_init();

	//Impressões para criação do filme
	fprintf(filme, "ANIMSTEPS %d\nCRYSTAL\nPRIMVEC\n", 1+npassos);//1+10*(int)(nt/amostragem));
	j=0;
	for(i=0; i<3; i++){
		fprintf(filme, " %lf\t%lf\t%lf\n", vetorRede[j]*latticeParameter*4e10, vetorRede[j+1]*latticeParameter*4e10, vetorRede[j+2]*latticeParameter*4e10);
		j+=3;
	}
	
	fprintf(filme, "PRIMCOORD %d\n%d 1\n", contFilme, tam);
	contFilme++;
	for(i=0; i<tam; i++){
		fprintf(filme, "  18\t%.8lf\t%.8lf\t%.8lf\t%.8lf\t%.8lf\t%.8lf\n", part[i].px*1e10, part[i].py*1e10, part[i].pz*1e10, part[i].fx*1e12, part[i].fy*1e12, part[i].fz*1e12);
		fprintf(teste, "%g\t%.8lf\t%.8lf\t%.8lf\n", t, part[i].vx/part[i].ax, part[i].vy/part[i].ay, part[i].vz/part[i].az);
	}


	metPotencial(0);//forces(t,1);
	metIntegracao();//veloc_ver(tam);

	//normTemp(tpt());
	//Loop nos primeiros passos de termalização
	for(t=0;t<amostragem;t+=dt){
		//Impressão das partículas durante a termalização
		//if((int)(t/dt)%ntemp == 0){
			fprintf(filme, "PRIMCOORD %d\n%d 1\n", contFilme, tam);
			contFilme++;
			for(i=0; i<tam; i++){
				fprintf(filme, "  18\t%.8lf\t%.8lf\t%.8lf\t%.8lf\t%.8lf\t%.8lf\n", part[i].px*1e10, part[i].py*1e10, part[i].pz*1e10, part[i].fx*1e12, part[i].fy*1e12, part[i].fz*1e12);
				fprintf(teste, "%g\t%.8lf\t%.8lf\t%.8lf\n", t, part[i].vx/part[i].ax, part[i].vy/part[i].ay, part[i].vz/part[i].az);
			}
		//}
		fprintf(energias_temp,"%.20g\t%.20g\t%.20g\t%.20g\n",t,e_K(),metPotencial(t),tpt());			
		metIntegracao();//veloc_ver(tam);
		//normTemp(tpt());

	}

	//Loop de Simulação e cálculo das propriedades
	inicializa_difusao(sqrt(e_K()/mass) );

	for(t=amostragem;t<nt;t+=dt){
		
	//Impressão das partículas após a termalização
		//if((int)(t/dt)%ntemp == 0){
			fprintf(filme, "PRIMCOORD %d\n%d 1\n", contFilme, tam);
			contFilme++;
			for(i=0; i<tam; i++){
				fprintf(filme, "  18\t%.8lf\t%.8lf\t%.8lf\t%.8lf\t%.8lf\t%.8lf\n", part[i].px*1e10, part[i].py*1e10, part[i].pz*1e10, part[i].fx*1e12, part[i].fy*1e12, part[i].fz*1e12);
			}
		//}
		fprintf(energias_temp,"%.20g\t%.20g\t%.20g\t%.20g\n",t,e_K(),metPotencial(t),tpt());			
		metIntegracao();
		//normTemp(tpt());

		//Cálculo da média das velocidades para G(r) e difusão
		if((int)(t/dt)%1==0){
			Gr_sample();
			contador++;
			diff_sample();
			diffusion_sample(sqrt(e_K()/mass),dt);
		}	
	}
	for(i=0; i<tam; i++)
		fprintf(teste, "%g\t%.8lf\t%.8lf\t%.8lf\n", t, part[i].vx, part[i].vy, part[i].vz);
	snprintf(nomeAux, 100, "T%.0lf-I%d-V%d-P%d.dat", temperatura_ideal, inpIntegracao, inpIniVel, inpPot);
	printf("Resultado de %s:\n", nomeAux);
	printf("Tempo de Simulação: %gs\n", nt);
	printf("Temperatura: %.1lfK\n", temperatura_ideal);
	if(inpIntegracao == 1){
		printf("Verlet (p=0.5) \n");
	}else if(inpIntegracao == 2){
		printf("Verlet (p=2.0) \n");
	}else if(inpIntegracao == 3){
		printf("Velocity Verlet (Algoritmo 1)\n");
	}else if(inpIntegracao == 4){
		printf("Velocity Verlet (Algoritmo 2)\n");
	}else if(inpIntegracao == 5){
		printf("Leap\n");
	}else{
		printf("Método de Integração Inválido\n");
	}

	if(inpIniVel == 0){
		printf("Velocidade Zero\n");
	}else if(inpIniVel == 1){
		printf("Velocidade Teo Equipartição\n");
	}else if(inpIniVel == 2){
		printf("Velocidade Boltzman\n");
	}else{
		printf("Inicialização da Velocidade Inválida\n");
	}

	if(inpPot == 1){
		printf("Lennard Jones\n");
	}else if(inpPot == 2){
		printf("Outro\n");
	}else{
		printf("Potencial Inválido\n");
	}

	printf("coeficiente de difusao por distancia media %.20g\n",diff_out(nt,amostragem) );
	printf("coeficiente de difusao por integracao de velocidade %.20g\n",Velocidade_Media );
	Gr_output();
	fclose(energias_temp);
	fclose(filme);
}


void metIntegracao(){
	if(inpIntegracao == 1){
		verlet(tam);
	}else if(inpIntegracao == 2){
		verl(tam);
	}else if(inpIntegracao == 3){
		velocity_verlet(tam);
	}else if(inpIntegracao == 4){
		veloc_ver(tam);
	}else if(inpIntegracao == 5){
		leap(tam);
	}else{
		printf("Método de Integração Inválido\n");
	}
}
void metIniVelocidade(){
	if(inpIniVel == 0){
		//printf("Velocidade Zero\n");
		iniVel0();
	}else if(inpIniVel == 1){
		//printf("Velocidade Teo Equipartição\n");
		iniVelTeoPart();
	}else if(inpIniVel == 2){
		//printf("Velocidade Boltzman\n");
		iniVelTeoPartgauss();
	}else{
		printf("Inicialização da Velocidade Inválida\n");
	}
}
double metPotencial(double t){
	/*if(inpPot == 1){
		//printf("Lennard Jones\n");
	}else if(inpPot == 2){
		//printf("Outro\n");
	}else{*/
	if(inpPot != 1 && inpPot != 2){
		printf("Potencial Inválido\n");
	}
	return forces(t,inpPot);
}



/*CÓDIGOS USADOS EM VERSÕES ANTERIORES
void main(){
	double t,nt=1e-8, t_tpt,amostragem=0.1*nt;
	int i, ntemp, j=0;
	double auxxx;
	ntemp = 10;
	FILE *energias_temp;
	FILE *filme;
	energias_temp=fopen("energias.dat","w");
	filme=fopen("filme.plot","w");
	redeFCC();
	iniVelTeoPart();
	iniAc0();
	copia_particula();
	primeira_integracao();
	//for(i=0; i<tam; i++)
	//	printf("%.20lf\t%.20lf\t%.20lf\n",part[i].px,part[i].py,part[i].pz);

	//iniVelTeoPartgauss();
		
	Gr_init();
	diff_init;
	normTemp(tpt());
	//velTeoPart();
	//puts("set terminal gif animate delay 100");
	//puts("set output 'particulas_particao.gif'");
	//puts("set xrange [35.039999748803964508:35.039975457981157092]");
	//puts("set yrange [7.71806191834435662:7.7180855457524373264]");
	//puts("set zrange [7.7180614817064068234:7.7180428859374323736]");
	//printf("Temperatura inicial: %.20lf\n", tpt());
	//for(i=0; i<tam; i++){
//		printf("%.20g\t%.20g\t%.20g\n",part[i].px,part[i].py,part[i].pz);
//	}
	//forces(0.);
			
	int contador=0;
	for(t=0;t<amostragem;t+=dt){
		verlet(tam);
		forces(t);
	}
	inicializa_difusao(sqrt(e_K()/mass) );
	for(t=amostragem;t<=nt;t+=dt){

		//puts("set xrange[0:53]");
		//puts("set yrange[0:53]");
		//puts("set zrange[0:53]");
		//printf("%g\t%g\t%d\n", t, dt, (int)(t/dt)%10);
		
		verlet(tam);
		
		fprintf(energias_temp,"%.20g\t%.20g\t%.20g\t%.20g\n",t,e_K(),forces(t),tpt());

		if((int)(t/dt)%1==0){

				Gr_sample();
					contador++;
				diff_sample();
				diffusion_sample(sqrt(e_K()/mass),dt);
		/*for(i=0; i<tam; i++){
			if(part[i].px > Lx || part[i].py > Ly || part[i].pz > Lz){
				printf("é > 52.56;\t %d\t%d\t%.2g\n",i, (int)(t/dt), t);
				printf("%.20g\t%.20g\t%.20g\n",part[i].px,part[i].py,part[i].pz);
				//scanf("%lf", &auxxx);
			}
			else if(part[i].px < 0 || part[i].py < 0 || part[i].pz < 0){
				printf("é < 00.00;\t %d\t%d\t%.2g\n",i, (int)(t/dt), t);
				printf("%.20g\t%.20g\t%.20g\n",part[i].px,part[i].py,part[i].pz);
				//scanf("%lf", &auxxx);
			}
		}* /

		//normTemp(tpt());
			//printf("%.20g\n", t_tpt);
			//printf("splot '-' title '%g\t%g'\n", t, t_tpt);
			//for(i=0;i<tam;i++){
			//	printf("%.20g\t%.20g\t%.20g\n",part[i].px,part[i].py,part[i].pz);
				//printf("%.20g\t%.20g\t%.20g\n",part[i].px,part[i].py,part[i].pz);
			//}
			//if(contador==54){
			//i=74;
		//	for(i=0;i<tam;i++){
			//if(i==111)
		//		printf("%.20g\t%.20g\t%.20g\t%d\n",part[i].fx,part[i].fy,part[i].fz,2);
			//else
			//	printf("%.20g\t%.20g\t%.20g\t%d\n",part[i].px,part[i].py,part[i].pz,1);
			
			//}
			//}	
			//puts("e");	
		}	
	}
	printf("coeficiente de difusao por distancia media %.20g\n",diff_out(nt,amostragem) );
	printf("coeficiente de difusao por integracao de velocidade %.20g\n",Velocidade_Media );
	Gr_output();
	//printf("diff %.20g\n",diff_out(dt,nt));
	printf("energias em energias.dat\nposicoes em posicao media.dat\nDisplacamente em displacement.dat (delete este arquivo antes de rodar outra simulação)\n");
	
	//for(i=0;i<tam;i++)	
	//	printf("%.20g\t%.20g\t%.20g\n",part[i].px,part[i].py,part[i].pz);
	//i=111;
				//printf("%.20g\t%.20g\t%.20g\t%d\n",part[i].px,part[i].py,part[i].pz,2);
}



*/
