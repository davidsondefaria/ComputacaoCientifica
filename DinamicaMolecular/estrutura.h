//#include "brunu.h" 	//para uso da função sign
#include "boxmuller.h"
//cubo

#define latticeParameter 5.256e-10
#define nAtoms 4

#define lpX 2.1024e-09	//4*5.256e-10
#define lpY 2.1024e-09
#define lpZ 2.1024e-09

#define Lx 8.4096e-09	//4*4*5.256e-10
#define Ly 8.4096e-09
#define Lz 8.4096e-09

//#define dLx 2.56e-10
#define dLx 0.0
#define dLy 0.0
#define dLz 0.0
#define mass 6.69e-26		//kg
#define ep 1.65e-21			//J 		//1.66E-14		//erg			//1.036091		//eV
#define SIGMA 3.405e-10		//m
#define tam 256
//#define dt 1e-11
#define k_b	1.3806503e-23	//J/K 	//1.3806488E-16	//erg/K			//8.6173324e-5	//eV/K

int inpIntegracao, inpIniVel, inpPot;
double temperatura_ideal;
double nt, dt;
double vetorRede[9];
double posvet[12];

//Funções de Input. Chamar na DM
void leInput(){
	int i, j;
	FILE *input;
	if(!(input = fopen("entrada.in", "r"))){
		printf("Verifique o arquivo de entrada.\n");
	}
	fscanf(input, "Temperatura (K): %lf\n", &temperatura_ideal);
	fscanf(input, "Método de Integração: 1/2. Verlet(p=0.5/p=2), 3/4. Velocity Verlet(Alg. 1/2), 5. Leap => %d\n", &inpIntegracao);
	fscanf(input, "Método de Inicialização de Velocidade: 0. Zero, 1. Teorema Equipartição, 2. Distribuição de Boltzman => %d\n", &inpIniVel);
	fscanf(input, "Potencial: 1. Lennard Jones, 2. Morse => %d\n", &inpPot);
	fscanf(input, "Tempo de Simulação: %lf\n", &nt);
	fscanf(input, "Passo Temporal: %lf\n", &dt);
	fscanf(input, "\nVetor da Rede:\n");
	j=0;
	for(i=0; i<3; i++){
		fscanf(input, "%lf\t%lf\t%lf\n", &vetorRede[j], &vetorRede[j+1], &vetorRede[j+2]);
		j+=3;
	}
	fscanf(input, "\nVetor de Posição:\n");
	j=0;
	for(i=0; i<4; i++){
		fscanf(input, "%lf\t%lf\t%lf\n", &posvet[j], &posvet[j+1], &posvet[j+2]);
		j+=3;
	}
	fclose(input);
}


static struct redeCristalina{
	double px, py, pz;	//posições
}rede[tam];

//Estrutura da partícula	
struct particula{
	double px, py, pz;	//posições
	double fx, fy, fz;	//forças
	double en;			//potencial
	double vx, vy, vz;	//velocidades
	double pm, pcg;		//massa, carga
	double ax, ay, az;	//aceleração
	int natom;			//num atômico	
};


struct particula part[tam];		//vetor de particulas
struct particula part2[tam];

/*//lista encadeada
struct plista{
	int conteudo;
	struct particula *prox;
};*/

//cria a rede e posiciona as particulas em cima
void redeFCC(){
	int i, j, k, w;
	int nx=4, ny=4, nz=4;	//num de translações
	int npcu=4;				//num particulas de celula unitaria
	int np;
	int N, M, mm;

	np = npcu*nx*ny*nz;
	N=npcu;
	mm=0;
	for(i=0; i<nx; i++){
		for(j=0; j<ny; j++){
			for(k=0; k<nz; k++){
				N = N+npcu;
				for(w=1; w<=npcu; w++){
					M = 3*w;
					//mm = N+w;

					//inicializa a rede
					rede[mm].px = ((i+posvet[M-3]))*latticeParameter;
					rede[mm].py = ((j+posvet[M-2]))*latticeParameter;
					rede[mm].pz = ((k+posvet[M-1]))*latticeParameter;
					//printf("%g\t%g\t%g\n",i+posvet[M-2],j+posvet[M-1],k+posvet[M]);
					//inicializa as posições das particulas na rede
					part[mm].px = rede[mm].px;
					part[mm].py = rede[mm].py;
					part[mm].pz = rede[mm].pz;
					part[mm].pm = mass;

					mm++;
				}
			}
		}
	}
}


void iniVel0(){
	int i;
	for(i=0; i<tam; i++){
		//inicializa as velocidades como zero
		part[i].vx = 0;
		part[i].vy = 0;
		part[i].vz = 0;
	}
}

void iniAc0(){
	int i;
	for(i=0; i<tam; i++){
		//inicializa as velocidades como zero
		part[i].ax = 0;
		part[i].ay = 0;
		part[i].az = 0;
	}
}
void primeira_integracao(){
int i;
for(i=0; i<tam; i++){
		//part[i].vx = vvx[i];
		//part[i].vy = vvy[i];
		//part[i].vz = vvz[i];
		part2[i].px += part[i].vx*dt;
		part2[i].py += part[i].vy*dt;
		part2[i].pz += part[i].vz*dt;
		part2[i].px=part2[i].px +(part2[i].px<0)*Lx -(part2[i].px>Lx)*Lx;
		part2[i].py=part2[i].py +(part2[i].py<0)*Ly -(part2[i].py>Ly)*Ly;
		part2[i].pz=part2[i].pz +(part2[i].pz<0)*Lz -(part2[i].pz>Lz)*Lz;
}
}
double *velTeoPart(int semente, double temp){
	int i;
	double sumv=0, sumv2=0;
	double *v;
	v = (double*)malloc(tam*sizeof(double));
	double fs;
	
	rinitialize(semente);

	for(i=0; i<tam; i++){
		v[i] = (ranf() - 0.5);
		sumv = sumv + v[i];
		sumv2 = sumv2 + v[i]*v[i];
	}

	sumv = sumv/tam;
	sumv2 = sumv2/tam;
	fs = sqrt(3*temp*k_b/(sumv2*mass));

	for(i=0; i<tam; i++){
		v[i] = (v[i] - sumv)*fs;
		//v[i] = (v[i] - sumv);
	}
	return v;
}
void iniVelTeoPart(){
	double *vvx, *vvy, *vvz, vv, vm = 0;
	int i;
	double temperatura = temperatura_ideal;

	vvx = velTeoPart(5023424, temperatura);
	vvy = velTeoPart(9576521, temperatura);
	vvz = velTeoPart(1632174, temperatura);

	for(i=0; i<tam; i++){
		part[i].vx = vvx[i];
		part[i].vy = vvy[i];
		part[i].vz = vvz[i];
		//part[i].px += part[i].vx*dt;
		//part[i].py += part[i].vy*dt;
		//part[i].pz += part[i].vz*dt;
		vv = sqrt(part[i].vx*part[i].vx + part[i].vy*part[i].vy + part[i].vz*part[i].vz);
		//printf("%.8g\t%.8g\t%.8g\t%lf\n", part[i].vx, part[i].vy, part[i].vz, vv);
		vm += vv/tam;
	}

	//printf("Velocidade Inicial Média: %.20lf\n", vm);
}

double *velTeoPartgauss(int semente,double ecin,double temp){
	int i;
	double *v;
	v = (double*)malloc(tam*sizeof(double));
	double fs;
	//double temp = 200;	//em Kelvin e deve ser input
	
	rinitialize(semente);

	for(i=0; i<tam; i++){
		v[i] = dist(sqrt(ecin),sqrt(temp) );
	}

	return v;
}

void iniVelTeoPartgauss(){
	double *vvx, *vvy, *vvz, vv, vm = 0;
	int i;

	vvx = velTeoPartgauss(5210284,200*0.5*k_b,200);
	vvy = velTeoPartgauss(97852102,200*0.5*k_b,200);
	vvz = velTeoPartgauss(16987421,200*0.5*k_b,200);

	for(i=0; i<tam; i++){
		part[i].vx = vvx[i];
		part[i].vy = vvy[i];
		part[i].vz = vvz[i];
		vv = sqrt(part[i].vx*part[i].vx + part[i].vy*part[i].vy + part[i].vz*part[i].vz);
		//printf("%lf\t%lf\t%lf\t%lf\n", part[i].vx, part[i].vy, part[i].vz, vv);
		vm += vv/tam;
	}

	//printf("%lf\n", vm);
}

void normTemp(double T){
	double temp = temperatura_ideal;
	double sumx = 0, sumy = 0, sumz = 0, sum2x = 0, sum2y = 0, sum2z = 0;
	int i;
	for(i=0; i<tam; i++){
		sumx = sumx + part[i].vx;
		sumy = sumy + part[i].vy;
		sumz = sumz + part[i].vz;
		sum2x = sum2x + part2[i].vx;
		sum2y = sum2y + part2[i].vy;
		sum2z = sum2z + part2[i].vz;
	}

	sumx = sumx/tam;
	sumy = sumy/tam;
	sumz = sumz/tam;
	sum2x = sum2x/tam;
	sum2y = sum2y/tam;
	sum2z = sum2z/tam;

	for (i=0; i<tam; i++){
		//printf("temp %g\n",sqrt(temp/T));
		part[i].vx = (part[i].vx - sumx)*sqrt(temp/T);
		part[i].vy = (part[i].vy - sumy)*sqrt(temp/T);
		part[i].vz = (part[i].vz - sumz)*sqrt(temp/T);
		part2[i].vx = (part2[i].vx - sum2x)*sqrt(temp/T);
		part2[i].vy = (part2[i].vy - sum2y)*sqrt(temp/T);
		part2[i].vz = (part2[i].vz - sum2z)*sqrt(temp/T);
	}

}

double tpt(){
	int i;
	double sum = 0, nf, T;
	nf = 3.*tam;


	for(i=0; i<tam; i++){
		sum += ((part[i].vx*part[i].vx + part[i].vy*part[i].vy + part[i].vz*part[i].vz));
	}

	T=mass*sum/(k_b*nf);
	//printf("%g\n", T);
	return T;
}




/*
void condContornoPeriodica(){
	double Lx2, Ly2, Lz2;
	Lx2 = Lx/2.;
	Ly2 = Ly/2.;
	Lz2 = Lz/2.;
	int i;

	for(i=0; i<tam; i++){
		part[i].px = part[i].px - sign(Lx2, part[i].px) - sign(Lx2, x-Lx);
		part[i].py = part[i].py - sign(Ly2, part[i].py) - sign(Ly2, y-Ly);
		part[i].pz = part[i].pz - sign(Lz2, part[i].pz) - sign(Lz2, z-Lz);
	}
}*/
