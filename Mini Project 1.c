#include<stdio.h>
#include<stdlib.h>
#include<math.h>

double Ufun(double u, double v)
{
    return(1/pow(u,4)-(1/u)-(1.5*v*v)/u);
}
double check(double u,double r)
{
    return((2.0/pow(u,3))*log(u/(1+r))-(2.0/3.0)*(1-pow((1+r)/u,3)));
}
void main()
{
    FILE *fp;
    fp=fopen("data.csv","w");
    double dt=0.0001,ndt,desdel,del,phi;
    int tMAX=50;
    double r=0.2,t,i,sum;
    double ustar,ustar2,ustar3,nu,u,su;
    double vstar,vstar2,vstar3,nv,v,sv;
    v=0;
    u=1+r;
    for(t=0;t<=tMAX;t+=dt)
    {
        nv=v;
        nu=u;
        for(i=0;i<2;i++)
        {
        vstar=nv+0.25*dt*Ufun(nu,nv);
        ustar=nu+0.25*dt*nv;

        vstar2=nv+0.25*dt*Ufun(ustar,vstar);
        ustar2=nu+0.25*dt*vstar;

        vstar3=nv+0.5*dt*Ufun(ustar2,vstar2);
        ustar3=nu+0.5*dt*vstar2;

        nv=nv+(dt/12)*(Ufun(nu,nv)+2*Ufun(ustar,vstar)+2*Ufun(ustar2,vstar2)+Ufun(ustar3,vstar3));
        nu=nu+(dt/12)*(nv+2*vstar+2*vstar2+vstar3);
        }
        vstar=v+0.5*dt*Ufun(u,v);
        ustar=u+0.5*dt*v;

        vstar2=v+0.5*dt*Ufun(ustar,vstar);
        ustar2=u+0.5*dt*vstar;

        vstar3=v+dt*Ufun(ustar2,vstar2);
        ustar3=u+dt*vstar2;

        sv=v+(dt/6)*(Ufun(u,v)+2*Ufun(ustar,vstar)+2*Ufun(ustar2,vstar2)+Ufun(ustar3,vstar3));
        su=u+(dt/6)*(v+2*vstar+2*vstar2+vstar3);

        fprintf(fp,"%lf,%lf\n",t,u);
        desdel=0.0001*(fabs(u)+dt*fabs(v));
        del=fabs(su-nu);
        phi=desdel/del;
        if (phi>1)
            dt*=pow(phi,0.2);
        if (phi<1)
            dt*=pow(phi,0.25);
        printf("%.5lf\t%.5lf\t%.5lf\t%.5lf\t%.5lf\t%.5lf\n",u,del,phi,desdel,nu,dt);
        v=sv;
        u=su;
    }
    fclose(fp);
}
