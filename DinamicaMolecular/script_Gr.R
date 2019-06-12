x = readChar("arqaux.dat", nchars=25)
a=read.table(x);
Dcell=52.56e-10
b=1
b[2]=1
Lx=52.56e-10
i=1
box=50
u=( (a[i,1]-a[,1]>0.5*Lx)*(Lx-a[i,1]+a[,1])^2+ (a[i,1]-a[,1]<=0.5*Lx)*(a[i,1]-a[,1])^2 + (a[i,2]-a[,2]>0.5*Lx)*(Lx-a[i,2]+a[,2])^2+ (a[i,2]-a[,2]<=0.5*Lx)*(a[i,2]-a[,2])^2 +(a[i,3]-a[,3]>0.5*Lx)*(Lx-a[i,3]+a[,3])^2+ (a[i,3]-a[,3]<=0.5*Lx)*(a[i,3]-a[,3])^2 )^0.5  #distancia
#u=hist(u,box,plot=0)
for(i in 2:length(a[,1]) ){
u[ ((i-1)*(256+1)):((i)*(256))]=( (a[i,1]-a[-i,1]>0.5*Lx)*(Lx-a[i,1]+a[-i,1])^2+ (a[i,1]-a[-i,1]<=0.5*Lx)*(a[i,1]-a[-i,1])^2 + (a[i,2]-a[-i,2]>0.5*Lx)*(Lx-a[i,2]+a[-i,2])^2+ (a[i,2]-a[-i,2]<=0.5*Lx)*(a[i,2]-a[-i,2])^2 +(a[i,3]-a[-i,3]>0.5*Lx)*(Lx-a[i,3]+a[-i,3])^2+ (a[i,3]-a[-i,3]<=0.5*Lx)*(a[i,3]-a[-i,3])^2 )^0.5  #distancia
}
#c=u
c=hist(u,box,plot=0)
c$breaks=c$breaks[1:length(c$breaks)-1]
delta=0.5*Dcell/box
#fact=(Dcell^3)/(2.0*pi*256^2*(delta)^3)
#for( i in 1:length(c$counts)){
#c$counts=c$counts/(i-0.5)**2
#}
#for (i=0;i<nB;i++) {
  #  r=dr*(i+0.5);
   # vb=((i+1)*(i+1)*(i+1)-i*i*i)*dr*dr*dr;
   # nid=(4./3.)*M_PI*vb*rho;
c$counts[1]=0
c$breaks=2*c$breaks*1e9/3.405
#for (i in 2:length(c$breaks)-1 ){
#	c$counts[i]=c$counts[i]/(256*0.1*4*3.1415*c$breaks[i]^2)

#}
v=data.frame(c$breaks,c$counts);
write.table(file=paste("G(r)-",x),v,row.names=F,col.names=F)
