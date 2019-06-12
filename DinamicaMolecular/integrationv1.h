//struct particula part[tam];


void copia_particula(){
	int i;
	for(i=0;i<tam;i++){
		part2[i].px=part[i].px;
		part2[i].py=part[i].py;
		part2[i].pz=part[i].pz;//posições
		part2[i].fx=part[i].fx;
		part2[i].fy=part[i].fy;
		part2[i].fz=part[i].fz; //forcas
		part2[i].vx=part[i].vx;
		part2[i].vy=part[i].vy;
		part2[i].vz=part[i].vz;
		part2[i].en=part[i].en;
		part2[i].pm=part[i].pm;
		part2[i].pcg=part[i].pcg;
		part2[i].ax=part[i].ax;
		part2[i].ay=part[i].ay;
		part2[i].az=part[i].az;
		part2[i].natom=part[i].natom;
	}
}

double integra_velocidade(double v,double a){
	return v+a*dt;

}

double integra_posicao(double v,double S,double a){
	return v*dt+0.5*a*dt*dt;
}

void copia(struct particula *p,struct particula *q){
	//struct particula aux;
	p->vx=q->vx;
	p->vy=q->vy;
	p->vz=q->vz;
	p->px=q->px;
	p->py=q->py;
	p->pz=q->pz;
}
double e_K(){
int i;
double sum=0;
for(i=0; i<tam; i++){
		sum += (part[i].pm*(part[i].vx*part[i].vx + part[i].vy*part[i].vy + part[i].vz*part[i].vz));
	}
return 0.5*sum;

}

struct particula verlet(int n){
	int i;
	struct particula *q;
	q=(struct particula *)malloc(sizeof(struct particula));
	q=part2;
	for(i=0;i<n;i++){
		q->vx=0.5*integra_velocidade(part2[i].vx,part2[i].ax)+0.5*integra_velocidade(part[i].vx,part[i].ax);
		q->vy=0.5*integra_velocidade(part2[i].vy,part2[i].ay)+0.5*integra_velocidade(part[i].vy,part[i].ay);
		q->vz=0.5*integra_velocidade(part2[i].vz,part2[i].az)+0.5*integra_velocidade(part[i].vz,part[i].az);

		q->px+=0.5*integra_posicao(part2[i].vx,part2[i].px,part2[i].ax)+0.5*integra_posicao(part[i].vx,part[i].px,part[i].ax);
		q->py+=0.5*integra_posicao(part2[i].vy,part2[i].py,part2[i].ay)+0.5*integra_posicao(part[i].vy,part[i].py,part[i].ay);
		q->pz+=0.5*integra_posicao(part2[i].vz,part2[i].pz,part2[i].az)+0.5*integra_posicao(part[i].vz,part[i].pz,part[i].az);


		/*q->px=q->px +(q->px<0)*Lx -(q->px>Lx)*Lx;
		q->py=q->py +(q->py<0)*Ly -(q->py>Ly)*Ly;
		q->pz=q->pz +(q->pz<0)*Lz -(q->pz>Lz)*Lz;*/

		if(q->px >= Lx){
			q->px = q->px - Lx-dLx;
		}else if(q->px < 0){
			q->px = q->px + Lx+dLx;
		}

		if(q->py >= Ly){
			q->py = q->py - Ly-dLy;
		}else if(q->py < 0){
			q->py = q->py + Ly+dLy;
		}

		if(q->pz >= Lz){
			q->pz = q->pz - Lz-dLz;
		}else if(q->pz < 0){
			q->pz = q->pz + Lz+dLz;
		}
		//copias

		part[i].vx=part2[i].vx;
		part[i].vy=part2[i].vy;
		part[i].vz=part2[i].vz;

		part[i].px=part2[i].px;
		part[i].py=part2[i].py;
		part[i].pz=part2[i].pz;

		part2[i].vx=q->vx;
		part2[i].vy=q->vy;
		part2[i].vz=q->vz;

		part2[i].px=q->px;
		part2[i].py=q->py;
		part2[i].pz=q->pz;
	}
}



struct particula velocity_verlet(int n){
	int i;
	struct particula *q;
	q=(struct particula *)malloc(sizeof(struct particula));
	for(i=0;i<n;i++){
		q->vx=0.5*integra_velocidade(part2[i].vx,part2[i].ax)+0.5*integra_velocidade(part[i].vx,part[i].ax);
		q->vy=0.5*integra_velocidade(part2[i].vy,part2[i].ay)+0.5*integra_velocidade(part[i].vy,part[i].ay);
		q->vz=0.5*integra_velocidade(part2[i].vz,part2[i].az)+0.5*integra_velocidade(part[i].vz,part[i].az);

		q->px=0.5*integra_posicao(part2[i].vx,part2[i].px,part2[i].ax)+0.5*integra_posicao(q->vx,part[i].px,part[i].ax);
		q->py=0.5*integra_posicao(part2[i].vy,part2[i].py,part2[i].ay)+0.5*integra_posicao(q->vy,part[i].py,part[i].ay);
		q->pz=0.5*integra_posicao(part2[i].vz,part2[i].pz,part2[i].az)+0.5*integra_posicao(q->vz,part[i].pz,part[i].az);
		
		/*q->px=q->px +(q->px<0)*Lx -(q->px>Lx)*Lx;
		q->py=q->py +(q->py<0)*Ly -(q->py>Ly)*Ly;
		q->pz=q->pz +(q->pz<0)*Lz -(q->pz>Lz)*Lz;*/

		if(q->px >= Lx){
			q->px = q->px - Lx-dLx;
		}else if(q->px < 0){
			q->px = q->px + Lx+dLx;
		}

		if(q->py >= Ly){
			q->py = q->py - Ly-dLy;
		}else if(q->py < 0){
			q->py = q->py + Ly+dLy;
		}

		if(q->pz >= Lz){
			q->pz = q->pz - Lz-dLz;
		}else if(q->pz < 0){
			q->pz = q->pz + Lz+dLz;
		}
		//copias

		part[i].vx=part2[i].vx;
		part[i].vy=part2[i].vy;
		part[i].vz=part2[i].vz;

		part[i].px=part2[i].px;
		part[i].py=part2[i].py;
		part[i].pz=part2[i].pz;

		part2[i].vx=q->vx;
		part2[i].vy=q->vy;
		part2[i].vz=q->vz;

		part2[i].px=q->px;
		part2[i].py=q->py;
		part2[i].pz=q->pz;
	}
}

void verl(int n){
	int i;
	struct particula *q;
	q=(struct particula *)malloc(sizeof(struct particula));
	for(i=0;i<n;i++){
		q->px=2*part2[i].px-part[i].px +part2[i].ax*dt*dt;
		q->py=2*part2[i].py-part[i].py +part2[i].ay*dt*dt;
		q->pz=2*part2[i].pz-part[i].pz +part2[i].az*dt*dt;

		if(q->px >= Lx){
			q->px = q->px - Lx;
		}else if(q->px < 0){
			q->px = q->px + Lx;
		}

		if(q->py >= Ly){
			q->py = q->py - Ly;
		}else if(q->py < 0){
			q->py = q->py + Ly;
		}

		if(q->pz >= Lz){
			q->pz = q->pz - Lz;
		}else if(q->pz < 0){
			q->pz = q->pz + Lz;
		}

		//q->px=q->px +(q->px<0)*Lx -(q->px>Lx)*Lx;
		//q->py=q->py +(q->py<0)*Ly -(q->py>Ly)*Ly;
		//q->pz=q->pz +(q->pz<0)*Lz -(q->pz>Lz)*Lz;

		part[i].vx=part2[i].vx;
		part[i].vy=part2[i].vy;
		part[i].vz=part2[i].vz;

		part[i].px=part2[i].px;
		part[i].py=part2[i].py;
		part[i].pz=part2[i].pz;

		//part2[i].vx=q->vx;
		//part2[i].vy=q->vy;
		//part2[i].vz=q->vz;

		part2[i].px=q->px;
		part2[i].py=q->py;
		part2[i].pz=q->pz;
	}
}

void veloc_ver(int n){
	int i;
	struct particula *q;
	q=(struct particula *)malloc(sizeof(struct particula));
	for(i=0;i<n;i++){
		q->px=part2[i].px+part2[i].vx*dt +0.5*part2[i].ax*dt*dt;
		q->py=part2[i].py+part2[i].vy*dt +0.5*part2[i].ay*dt*dt;
		q->pz=part2[i].pz+part2[i].vz*dt +0.5*part2[i].az*dt*dt;

		/*q->px=q->px +(q->px<0)*Lx -(q->px>Lx)*Lx;
		q->py=q->py +(q->py<0)*Ly -(q->py>Ly)*Ly;
		q->pz=q->pz +(q->pz<0)*Lz -(q->pz>Lz)*Lz;
*/
		if(q->px > Lx){
			q->px = q->px - Lx;
		}else if(q->px < 0){
			q->px = q->px + Lx;
		}

		if(q->py > Ly){
			q->py = q->py - Ly;
		}else if(q->py < 0){
			q->py = q->py + Ly;
		}

		if(q->pz > Lz){
			q->pz = q->pz - Lz;
		}else if(q->pz < 0){
			q->pz = q->pz + Lz;
		}

		q->vx=part2[i].vx+(0.5*part2[i].ax +0.5*part[i].ax)*dt;
		q->vy=part2[i].vy+(0.5*part2[i].ay +0.5*part[i].ay)*dt;		
		q->vz=part2[i].vz+(0.5*part2[i].az +0.5*part[i].az)*dt;
		part[i].ax=part2[i].ax;
		part[i].ay=part2[i].ay;
		part[i].az=part2[i].az;

		part[i].vx=part2[i].vx;
		part[i].vy=part2[i].vy;
		part[i].vz=part2[i].vz;

		part[i].px=part2[i].px;
		part[i].py=part2[i].py;
		part[i].pz=part2[i].pz;

		part2[i].vx=q->vx;
		part2[i].vy=q->vy;
		part2[i].vz=q->vz;

		part2[i].px=q->px;
		part2[i].py=q->py;
		part2[i].pz=q->pz;
	}
}

//leap frog

void leap(int n){
	struct particula *q;
	int i;
	q=(struct particula *)malloc(sizeof(struct particula));
	for(i=0;i<n;i++){
		q->px=part2[i].px+ dt*part2[i].vx;
		q->vx=(part2[i].px-part[i].px)/dt +dt*part2[i].ax;
		q->py=part2[i].py+ dt*part2[i].vy;
		q->vy=(part2[i].py-part[i].py)/dt +dt*part2[i].ay;
		q->pz=part2[i].pz+ dt*part2[i].vz;
		q->vz=(part2[i].pz-part[i].pz)/dt +dt*part2[i].az; 

		//q->px=q->px +(q->px<0)*Lx -(q->px>Lx)*Lx;
		//q->py=q->py +(q->py<0)*Ly -(q->py>Ly)*Ly;
		//q->pz=q->pz +(q->pz<0)*Lz -(q->pz>Lz)*Lz;

		if(q->px > Lx){
			q->px = q->px - Lx;
		}else if(q->px < 0){
			q->px = q->px + Lx;
		}

		if(q->py > Ly){
			q->py = q->py - Ly;
		}else if(q->py < 0){
			q->py = q->py + Ly;
		}

		if(q->pz > Lz){
			q->pz = q->pz - Lz;
		}else if(q->pz < 0){
			q->pz = q->pz + Lz;
		}

		part[i].vx=part2[i].vx;
		part[i].vy=part2[i].vy;
		part[i].vz=part2[i].vz;

		part[i].px=part2[i].px;
		part[i].py=part2[i].py;
		part[i].pz=part2[i].pz;

		part2[i].vx=q->vx;
		part2[i].vy=q->vy;
		part2[i].vz=q->vz;

		part2[i].px=q->px;
		part2[i].py=q->py;
		part2[i].pz=q->pz;

	}
}
/*


double E_k(){
int i;
double e=0;
for(i=0;i<tam;i++)
	e+=0.5*part2[i].m* (pow(part2[i].vx,2)+pow(part2[i].vy,2),pow(part2[i].vz,2) );

return e/n;


}

double temperatura(int graus){

return 2*E_k()/(graus*k); //k=constante de boltzman

}*/

/*
struct particula *velocity_verlet(struct particula *p1,struct particula *p2,int n){
//p1 se torna a particula no tempo presente enquanto p2 no tempo futuro n sera o numero de particulas da rede
struct particula q;
int i;
for(i=0;i<n;i++){
q->vx=p2[i]->vx+0->5*integra_v(p2[i]->vx,p2[i]->ax)+0->5*integra_v(p1[i]->vx,p1[i]->ax);
q->vy=p2[i]->vy+0->5*integra_v(p2[i]->vy,p2[i]->ay)+0->5*integra_v(p1[i]->vy,p1[i]->ay);
q->vz=p2[i]->vz+0->5*integra_v(p2[i]->vz,p2[i]->az)+0->5*integra_v(p1[i]->vz,p1[i]->az);

q->px=p2[i]->px+0->5*integra_p(p2[i]->px,p2[i]->px,p2[i]->ax)+0->5*integra_p(q->vx,p1[i]->px,p1[i]->ax);
q->py=p2[i]->py+0->5*integra_p(p2[i]->py,p2[i]->py,p2[i]->ay)+0->5*integra_p(q->vy,p1[i]->py,p1[i]->ay);
q->pz=p2[i]->pz+0->5*integra_p(p2[i]->pz,p2[i]->pz,p2[i]->az)+0->5*integra_p(q->vz,p1[i]->pz,p1[i]->az);
//inicio da copia
copia(p1,p2);
copia(p2,&q);
}
return &q;
}

struct particula *leap_frog(struct particula *p1,struct particula *p2,int n){
//p1 se torna a particula no tempo presente enquanto p2 no tempo futuro n sera o numero de particulas da rede
struct particula q;
int i;
for(i=0;i<n;i++){
q->vx=p2[i]->vx+0->25*integra_v(p2[i]->vx,p2[i]->ax)+0->25*integra_v(p1[i]->vx,p1[i]->ax);
q->vy=p2[i]->vy+0->25*integra_v(p2[i]->vy,p2[i]->ay)+0->25*integra_v(p1[i]->vy,p1[i]->ay);
q->vz=p2[i]->vz+0->25*integra_v(p2[i]->vz,p2[i]->az)+0->25*integra_v(p1[i]->vz,p1[i]->az);

q->px=p2[i]->px+0->5*integra_p(p2[i]->px,p2[i]->px,p2[i]->ax)+0->5*integra_p(q->vx,p1[i]->px,p1[i]->ax);
q->py=p2[i]->py+0->5*integra_p(p2[i]->py,p2[i]->py,p2[i]->ay)+0->5*integra_p(q->vy,p1[i]->py,p1[i]->ay);
q->pz=p2[i]->pz+0->5*integra_p(p2[i]->pz,p2[i]->pz,p2[i]->az)+0->5*integra_p(q->vz,p1[i]->pz,p1[i]->az);
//inicio da copia
copia(p1,p2);
copia(p2,&q);
}
return &q;
}
*/

