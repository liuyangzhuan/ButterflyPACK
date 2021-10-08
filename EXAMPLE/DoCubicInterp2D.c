#include "math.h"
#include "stdlib.h"
#include "stdio.h" 
#include "DoCubicInterp2D.h"


double kernelu(double s)
{
  double val=0.0;
  double ss=fabs(s);  
  if(ss<=1.0){
    val=1.5*ss*ss*ss-2.5*ss*ss+1.0;
  }
  else if(ss>1.0 && ss<2.0){
    val=-0.5*ss*ss*ss+2.5*ss*ss-4.0*ss+2.0;
  } 
  else{
    val=0.0;
  }

  /*
  if(ss<1.0e-10){
    val=1;
  }
  else if(ss>1.0e-10 && ss<1-1.0e-10){
    val=1.5*ss*ss*ss-2.5*ss*ss+1;
  }
  else if(ss>1+1.0e-10 && ss<2-1.0e-10){
    val=-0.5*ss*ss*ss+2.5*ss*ss-4*ss+2;
  }
  else if(ss>2+1.0e-10){
    val=0;
  }
  else{
    val=0;
  }
  //*/

  return val;
}

void  CubicInterp2D(double *x, double *y, double *u, int Nx, int Ny, double *xx1, double *yy1, double *v, int No)
{
  //Given domain
  double xmin=x[0], xmax=x[Nx-1];
  double ymin=y[0], ymax=y[Ny-1];
  //printf("\t Domain = [%f %f %f %f %f %f]\n",xmin,xmax,ymin,ymax,zmin,zmax);
  //Given meshsize
  double hx=(xmax-xmin)/(Nx-1);
  double hy=(ymax-ymin)/(Ny-1);
  //printf("\t Meshsize = [%f %f %f]\n",hx,hy,hz);
  //dummy variables
  int ix, iy;
  double x0, y0;
  int iix, iiy; 

  //printf("Nx Ny Nz = [%d %d %d]\n",Nx,Ny,Nz); 
  //printf("N1x N1y N1z = [%d %d %d]\n",N1x,N1y,N1z); 

  //map to unit length h.
  //added by one to fit into the later extension
  //
  double *x1=(double *)malloc(No*sizeof(double));
  double *y1=(double *)malloc(No*sizeof(double));
//  printf("\t Map to unit length ******\n");
  for(ix=0;ix<No;ix++){ 
     x1[ix]=xx1[ix];  
     x1[ix]=(x1[ix]-xmin)/hx;
     if(x1[ix]<=0){ x1[ix]=0;}
     if(x1[ix]-(Nx-1)>=0.0){x1[ix]=Nx-1;}
     x1[ix]=x1[ix]+1;
     //printf("%.20lf\n",x1[ix]);
  }
  //return;
  for(iy=0;iy<No;iy++){ 
     y1[iy]=yy1[iy];  
     y1[iy]=(y1[iy]-ymin)/hy;
     if(y1[iy]<=0){ y1[iy]=0;}
     if(y1[iy]-(Ny-1)>=0.0){y1[iy]=Ny-1;}
     y1[iy]=y1[iy]+1;
     //printf("%.20lf\n",y1[iy]);
  }
  //return; 

  //extension to boundary by 1 layer on each end
//  printf("\t Extension to one more layer ******\n");
  int Nx2=Nx+2;
  int Ny2=Ny+2;
  double *uu=(double *)malloc(Nx2*Ny2*sizeof(double)); 
  for(ix=0;ix<Nx2*Ny2;ix++){uu[ix]=0;} 
  int id=0;
  int id2=0;
  //data in the interior
   for(iy=1;iy<Ny2-1;iy++){
      for(ix=1;ix<Nx2-1;ix++){ 
         id=ix+Nx2*iy; id2=ix-1+Nx*(iy-1);
         if(id>=Nx2*Ny2 || id2>=Nx*Ny){ printf("%d %d\n",id,id2);}
         uu[id]=u[id2];  
      }
   }
  //*
  //data on the extended edges
  //To edges
//  printf("\t Data on edges ******\n");
  for(iy=0;iy<Ny2;iy++){     
        ix=0;
        id=ix+Nx2*iy; 
        uu[id]=uu[ix+3+Nx2*iy]-3*uu[ix+2+Nx2*iy]+3*uu[ix+1+Nx2*iy];

        ix=Nx2-1;
        id=ix+Nx2*iy; 
        uu[id]=3*uu[ix-1+Nx2*iy]-3*uu[ix-2+Nx2*iy]+uu[ix-3+Nx2*iy];
  } 
   
  for(ix=0;ix<Nx2;ix++){ 
        iy=0; 
        id=ix+Nx2*iy; 
        uu[id]=uu[ix+Nx2*(iy+3)]-3*uu[ix+Nx2*(iy+2)]+3*uu[ix+Nx2*(iy+1)];

        iy=Ny2-1; 
        id=ix+Nx2*iy; 
        uu[id]=3*uu[ix+Nx2*(iy-1)]-3*uu[ix+Nx2*(iy-2)]+uu[ix+Nx2*(iy-3)];
  }
  

  
  //To corners
//  printf("\t Data on corners ******\n");
  iy=0; ix=0;
  id=ix+Nx2*iy; 
  uu[id]=uu[ix+3+Nx2*iy]-3.0*uu[ix+2+Nx2*iy]+3.0*uu[ix+1+Nx2*iy];

  iy=0; ix=Nx2-1;
  id=ix+Nx2*iy; 
  uu[id]=3*uu[ix-1+Nx2*iy]-3*uu[ix-2+Nx2*iy]+uu[ix-3+Nx2*iy];

  iy=Ny2-1; ix=0;
  id=ix+Nx2*iy; 
  uu[id]=uu[ix+3+Nx2*iy]-3*uu[ix+2+Nx2*iy]+3*uu[ix+1+Nx2*iy];
  
  iy=Ny2-1; ix=Nx2-1;
  id=ix+Nx2*iy; 
  uu[id]=3*uu[ix-1+Nx2*iy]-3*uu[ix-2+Nx2*iy]+uu[ix-3+Nx2*iy];
  //*/
   
  //convolution 
//  printf("\t Do convolution ******\n");
  double ux[4], uy[4];
  int ii, jj; 
  int i, j;
  double iiix, iiiy; 
  int iii, jjj;
  for(id=0;id<No;id++){
     //
      iix=(int)(floor(x1[id])); if(iix>=Nx){iix=Nx-1;}
      iiix=x1[id]-iix;  
      ux[0]=kernelu(iiix+1);
      ux[1]=kernelu(iiix);
      ux[2]=kernelu(iiix-1);
      ux[3]=kernelu(iiix-2);

      iiy=(int)(floor(y1[id])); if(iiy>=Ny){iiy=Ny-1;}
      iiiy=y1[id]-iiy;
      uy[0]=kernelu(iiiy+1);
      uy[1]=kernelu(iiiy);
      uy[2]=kernelu(iiiy-1);
      uy[3]=kernelu(iiiy-2);

      //convolution
      v[id]=0; 
      for(ii=-1;ii<=2;ii++){
        i=iix+ii; 
        iii=ii+1;
        for(jj=-1;jj<=2;jj++){ 
            j=iiy+jj; 
            jjj=jj+1;;  
            v[id]=v[id]+uu[i+Nx2*j]*uy[jjj]*ux[iii]; 
        }
      }
      //       
  } 

  free(uu); free(x1); free(y1);
 
  return;
}
