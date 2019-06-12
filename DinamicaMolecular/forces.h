
#include  <math.h>


double forces (double t,int type) {
	int i,j,cont=0;
	double Fx[256],Fy[256],Fz[256],F[256],fx=0.0,fy=0.0,fz=0.0,x=0.0,y=0.0,z=0.0,xFut=0.0,yFut=0.0,zFut=0.0,r2=0.0,Fij=0.0,rc=0.0,r2i=0.0,r6=0.0,en=0.0;
	rc=Lx/2;
	double rneg=0.5*Lx*0.2;
	iniAc0();
	//rc = 2.37343*SIGMA;
	Fij=1.0;
	for(i=0;i<256;i++){
		F[i]=0.0;
		Fx[i]=0.0;
		Fy[i]=0.0;
		Fz[i]=0.0;
		part[i].fx = 0.0;
		part[i].fy = 0.0;
		part[i].fz = 0.0;
		part2[i].fx = 0.0;
		part2[i].fy = 0.0;
		part2[i].fz = 0.0;
		xFut=0.0;
		yFut=0.0;
		zFut=0.0;
	}
	
	double ecut=4*ep*(pow(SIGMA/rc,12)-pow(SIGMA/rc,6) );

	for (i=0;i<256;i++){
		x=part[i].px;
		y=part[i].py;
		z=part[i].pz;
	
		 for(j=i+1;j<256;j++){

			xFut=-x+part[j].px; 
		 	yFut=-y+part[j].py;
		 	zFut=-z+part[j].pz;
		 	
		 	//minima imagem
			if(xFut<-0.5*Lx+dLx){
				xFut=(part[j].px-x+Lx+dLx);
				//printf("1\n");
			}
			else if(xFut>0.5*Lx+dLx){
				xFut=(part[j].px-x-Lx-dLx);
				//printf("2\n");
			}
			if(yFut<-0.5*Ly+dLy){
				yFut=(part[j].py-y+Ly+dLy);
				//printf("3\n");
			}
			else if(yFut>0.5*Ly+dLy){
				yFut=(part[j].py-y-Ly-dLy);
				//printf("4\n");
			}
			if(zFut<-0.5*Lz+dLz){
				zFut=(part[j].pz-z+Lz+dLz);
				//printf("5\n");
			}
			else if(zFut>0.5*Lz+dLz){
				zFut=(part[j].pz-z-Lz-dLz);
				//printf("6\n");
			}
			
			r2=xFut*xFut+yFut*yFut+zFut*zFut;

			if(r2<rc*rc&&r2>(rneg*rneg)){
				if(type==1){  //1=lennard Jones
					Fij=(48.0*ep/r2)*(pow(SIGMA*SIGMA/r2,6) - 0.5*pow(SIGMA*SIGMA/r2,3)); 
					 // Fij=(48.0*SIGMA*SIGMA/r2)*(pow(SIGMA*SIGMA/r2,6) - 0.5*pow(SIGMA*SIGMA/r2,3));
				} else if(type==2){  //2=Potencial Morse
					//Fij= D*((c/v)*exp(-1.0*(c/v)*sqrt(r2))/sqrt(r2) - (c/v)*exp(-2.0*(c/v))*sqrt(r2)/sqrt(r2));
				} else {  // Colombiano
					//Fij = k*q*q/r2;
				}
				
				
				fx=Fij*xFut;
				fy=Fij*yFut;
				fz=Fij*zFut;
				part[j].fx=part[j].fx-fx;
				part[j].fy=part[j].fy-fy;
				part[j].fz=part[j].fz-fz;
				part[j].ax=((part[j].fx/part[j].pm));
				part[j].ay=((part[j].fy/part[j].pm));
				part[j].az=((part[j].fz/part[j].pm));
				
									
					
				//en=en-4*r6*(r6-1.0)-ecut;

				if(type==1){  //1=lennard Jones
					en=en+4.0*ep*(pow(SIGMA*SIGMA/r2,3))*(pow(SIGMA*SIGMA/r2,3) -1.0) - ecut;	 //potencial lennard jones
				} else if(type==2){
					//en = en+ D*(exp(-2.0*(c/v))*sqrt(r2)) -2.0*exp(-1.0*(c/v)*sqrt(r2) ) -ecut;   //potencial morse	
				} else {  // Colombiano
					//en = en+ k*q*q/sqrt(r2)   -ecut; // mas temos cargas?
				}
			}
			//if(i==127 || j==127){
		
				//printf("%d,%d \t %.2g \t %.20g \t %.20g \t \t %.20g \n",i,j, t, part[j].fx,part[j].fy,part[j].fz);
			//}

		}
		
		part[i].fx=part[i].fx+fx;
		part[i].fy=part[i].fy+fy;
		part[i].fz=part[i].fz+fz;
		//printf("%g %g %g \n",part[i].fx,part[i].fy,part[i].fz);
		part[i].ax=((part[i].fx/part[i].pm));
		part[i].ay=((part[i].fy/part[i].pm));
		part[i].az=((part[i].fz/part[i].pm));
		
	}
	
	


 	for(i=0;i<256;i++){
		
		
		printf("%g \t %g \t %g \t %g \t \n ", part[i].ax, part[i].ay,part[i].az, sqrt(part[i].fx*part[i].fx+part[i].fy*part[i].fy+part[i].fz*part[i].fz ));
	}

	return en;
}



