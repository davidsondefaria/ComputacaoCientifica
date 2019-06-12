// calculo da difusão
double Velocidade_Media;
double V_0;
double mean_square_displacamente;
struct particula Gr[tam];
int amostragem=0;

//D=<mean_square_displacamente>
void diff_init(){
	mean_square_displacamente=0;
}

void diff_sample(){
	int i;
	char nomeAux[100];
	double sum=0;
	FILE *arq;
	snprintf(nomeAux, 100, "disp-T%.0lf-I%d-V%d-P%d.dat", temperatura_ideal, inpIntegracao, inpIniVel, inpPot);
	arq=fopen(nomeAux,"a");
	for(i=0;i<tam;i++){
		sum+= pow(part2[i].px-part[i].px,2)+ pow(part2[i].py-part[i].py,2)+pow(part2[i].pz-part[i].pz,2);
	}
	sum=sum/tam;
	mean_square_displacamente+=sum*dt;
	fprintf(arq,"%.20g\n",mean_square_displacamente);
	fclose(arq);
}

double diff_out(double T,double t){
	return mean_square_displacamente/(6*(T-t));
}

void inicializa_difusao(double V){
	Velocidade_Media=0.0;
	V_0=V;
}

void diffusion_sample(double V,double deltat){
	Velocidade_Media+=V_0*V*deltat;
}



void Gr_init(){
	int i;
	for(i=0;i<tam;i++){
		Gr[i].px=0;
		Gr[i].py=0;
		Gr[i].pz=0;
	}
}

void Gr_sample(){
	int i;
	for(i=0;i<tam;i++){
		Gr[i].px+=part2[i].px;
		Gr[i].py+=part2[i].py;
		Gr[i].pz+=part2[i].pz;
	}
	amostragem++;
}

void Gr_output(){
	int i;
	char nomeAux[100];
	snprintf(nomeAux, 100, "PosMed-T%.0lf-I%d-V%d-P%d.dat", temperatura_ideal, inpIntegracao, inpIniVel, inpPot);
	FILE *posicoes;
	posicoes=fopen(nomeAux,"w");
	for(i=0;i<tam;i++){
		Gr[i].px=Gr[i].px/amostragem;
		Gr[i].py=Gr[i].py/amostragem;
		Gr[i].pz=Gr[i].pz/amostragem;
		fprintf(posicoes,"%.20g\t%.20g\t%.20g\n",Gr[i].px,Gr[i].py,Gr[i].pz);
	}
	FILE *aux;
	aux = fopen("arqaux.dat", "w");
	fprintf(aux, "%s", nomeAux);
}

