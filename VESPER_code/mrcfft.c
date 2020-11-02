#include <fftw3.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <stdbool.h>
#include <math.h>
#include <omp.h>
#include <float.h>
#include "struct.h"
#include "mrc.h"
#include "mrcfft.h"

#define ZERO 0.000000

extern CMD cmd;

//ZXY
void Euler2mtxZYX(double a,double b, double c,double m[3][3]){

 m[0][0]= cos(b)*cos(c);
 m[0][1]= sin(a)*sin(b)*cos(c)-cos(a)*sin(c);
 m[0][2]= cos(a)*sin(b)*cos(c)+sin(a)*sin(c);

 m[1][0]= cos(b)*sin(c);
 m[1][1]= cos(a)*cos(c)+sin(a)*sin(b)*sin(c);
 m[1][2]= cos(a)*sin(b)*sin(c)-sin(a)*cos(c);

 m[2][0]=-sin(b);
 m[2][1]= sin(a)*cos(b);
 m[2][2]= cos(a)*cos(b);
/*
 printf("%f %f %f\n",m[0][0],m[0][1],m[0][2]);
 printf("%f %f %f\n",m[1][0],m[1][1],m[1][2]);
 printf("%f %f %f\n",m[2][0],m[2][1],m[2][2]);
 puts("//");
*/
}

void Euler2mtxXYZ(double a,double b, double c,double m[3][3]){

 m[0][0]=  cos(b)*cos(c);
 m[0][1]=  cos(a)*sin(c)+sin(a)*sin(b)*cos(c);
 m[0][2]=  sin(a)*sin(c)-cos(a)*sin(b)*cos(c);

 m[1][0]= -cos(b)*sin(c);
 m[1][1]= cos(a)*cos(c)-sin(a)*sin(b)*sin(c);
 m[1][2]= sin(a)*cos(c)+cos(a)*sin(b)*sin(c);

 m[2][0]= sin(b);
 m[2][1]= -sin(a)*cos(b);
 m[2][2]= cos(a)*cos(b);

}

void mtxmtx(double a[3][3],double b[3][3],double c[3][3]){
 //puts("Product");
 for(int i=0;i<3;i++){
  for(int j=0;j<3;j++){
   c[i][j]=a[i][0]*b[0][j]+a[i][1]*b[1][j]+a[i][2]*b[2][j];
  //printf("%f ",c[i][j]);
  }
  //puts("");
 }
}


void rot_cd(double in[3],double out[3],double m[3][3]){

 out[0]=m[0][0]*in[0]+m[0][1]*in[1]+m[0][2]*in[2];
 out[1]=m[1][0]*in[0]+m[1][1]*in[1]+m[1][2]*in[2];
 out[2]=m[2][0]*in[0]+m[2][1]*in[1]+m[2][2]*in[2];

}

void inverse_mtx(double a[3][3],double b[3][3]){
 double det=0;
 int i,j;
 for(i=0;i<3;i++){
  for(j=0;j<3;j++){
   b[i][j]=a[j][i];
  }
 }
}


void RotMRC(MRC *in,MRC *out,double rx, double ry, double rz){

 int xydim=in->xdim*in->ydim;
 int xdim=in->ydim;
 int ind;
 double mtx[3][3],iv_mtx[3][3],tmp[3][3];
 double cent[3];

 out->std=in->std;
 out->ave=in->ave;
 out->std_norm_ave=in->std_norm_ave;

 cent[0]=0.5*(double)(in->xdim);
 cent[1]=0.5*(double)(in->ydim);
 cent[2]=0.5*(double)(in->zdim);

 Euler2mtxZYX(RAD(rx),RAD(ry),RAD(rz),mtx);
 //Euler2mtxZYX(RAD(360-rx),RAD(ry),RAD(rz),iv_mtx);

 inverse_mtx(mtx,iv_mtx);

 mtxmtx(mtx,iv_mtx,tmp);

 //printf("rote %f %f %f\n",rx,ry,rz);
 //Euler2mtxZYX(0,0,0,mtx);
 
 //#pragma omp parallel for schedule(dynamic,5)
 for(int x=0;x<in->xdim;x++){
  double pos[3],pos2[3],pos3[3];
  double vec[3],vec2[3];
  int nx,ny,nz,ind2;
  pos[0]=(double)x-cent[0];
 for(int y=0;y<in->ydim;y++){
  pos[1]=(double)y-cent[1];
 for(int z=0;z<in->zdim;z++){
  pos[2]=(double)z-cent[2];
  ind=xydim*z+xdim*y+x;
  //convert to relative position from center
  //rotate
  rot_cd(pos,pos2,iv_mtx); 
  //rot_cd(pos2,pos3,mtx); 

  //printf("pos1= %.1f %.1f %.1f\n",pos[0],pos[1],pos[2]);
  //printf("pos2= %.1f %.1f %.1f\n",pos2[0],pos2[1],pos2[2]);
  //printf("pos3= %.1f %.1f %.1f\n",pos3[0],pos3[1],pos3[2]);

  pos2[0]+=cent[0];
  pos2[1]+=cent[1];
  pos2[2]+=cent[2];


  ind2=xydim*(int)(pos2[2])+xdim*(int)pos2[1]+(int)pos2[0];
  //outside
  if(pos2[0]<0||pos2[1]<0||pos2[2]<0||
     pos2[0]>=in->xdim||pos2[1]>=in->ydim||pos2[2]>=in->zdim){

   out->dens[ind]=0.00;
   out->vec[ind][0]=0.00;
   out->vec[ind][1]=0.00;
   out->vec[ind][2]=0.00;
   continue;
  }
  if(in->dens[ind2]==0.00){

   out->dens[ind]=0.00;
   out->vec[ind][0]=0.00;
   out->vec[ind][1]=0.00;
   out->vec[ind][2]=0.00;
   continue;
  }
  

  //printf("*pos2= %.1f %.1f %.1f %d\n",pos2[0],pos2[1],pos2[2],ind2);


  vec[0]=in->vec[ind2][0]; 
  vec[1]=in->vec[ind2][1]; 
  vec[2]=in->vec[ind2][2]; 

  rot_cd(vec,vec2,mtx);

  //printf("%d %d %f %f\n",ind,ind2,vec[0],vec2[0]);

  out->dens[ind]=in->dens[ind2];
  out->vec[ind][0]=vec2[0];
  out->vec[ind][1]=vec2[1];
  out->vec[ind][2]=vec2[2];

 }}}


}


void Complex_Complex(fftwf_complex *a, fftwf_complex *b, fftwf_complex *rt,int n) {
 //#pragma omp parallel for schedule(dynamic,5)
 for(int i=0;i<n;i++){
    rt[i][0] = a[i][0] * b[i][0] - a[i][1] * b[i][1];
    rt[i][1] = a[i][0] * b[i][1] + a[i][1] * b[i][0];
 }
    return;
}


float FindBestTrans(float *x,float *y,float *z,int n,int trans[3]){

 float best=0.000;
 double s;
 int ind;//,bpos[3];
 int n2=n*n;

 trans[0]=trans[1]=trans[2]=0;
 for(int i=0;i<n;i++){
 for(int j=0;j<n;j++){
 for(int k=0;k<n;k++){
  ind=n2*k+n*j+i;
  s=x[ind]+y[ind]+z[ind];
  if(best<s){
   best=s;
   trans[0]=i;
   trans[1]=j;
   trans[2]=k;
  }
 }}}
 return best;
}

float FindBestTrans1D(float *x,int n,int trans[3]){

 double best=0;
 double s;
 int ind;//,bpos[3];
 int n2=n*n;

 //bpos[0]=bpos[1]=bpos[2]=0;
 for(int i=0;i<n;i++){
 for(int j=0;j<n;j++){
 for(int k=0;k<n;k++){
  ind=n2*k+n*j+i;
  //s=x[ind]+y[ind]+z[ind];
  s=x[ind];
  if(best<s){
   best=s;
   trans[0]=i;
   trans[1]=j;
   trans[2]=k;
  }
 }}}
 return best;
}


int cmp_tbl(const void *c1, const void *c2){
 TBL a=*(TBL *)c1;
 TBL b=*(TBL *)c2;

 if(a.sco>b.sco) return -1;
 if(a.sco<b.sco) return 1;
 return 0;
}

int cmp_tbl_code(const void *c1, const void *c2){
 TBL a=*(TBL *)c1;
 TBL b=*(TBL *)c2;

 if(a.code>b.code) return -1;
 if(a.code<b.code) return 1;
 return 0;
}


bool SearchMAPfftMT(MRC *m1,MRC *m2,double ang,bool Emode){

 MRC mtmp,*MT_mtmp;
 int n=m2->xdim;
 int xydim=m2->xdim*m2->ydim;
 int xyzdim=m2->xdim*m2->ydim*m2->zdim;
 int Nth=omp_get_max_threads();
 int TopN=cmd.TopN;

 if(Emode){//Evaluation mode
  TopN=0;
  if(m1->cent[0]!=m2->cent[0]
   ||m1->cent[1]!=m2->cent[1]
   ||m1->cent[2]!=m2->cent[2]
   ||m1->orgxyz[0]!=m2->orgxyz[0]
   ||m1->orgxyz[1]!=m2->orgxyz[1]
   ||m1->orgxyz[2]!=m2->orgxyz[2]
   ){
   puts("#For Evaluation Mode");
   puts("#Please use the same coordinate system and map size for map1 and map2.");
   puts("#Example:");
   puts("#In Chimera command line: open map1 and map2 as #0 and #1, then type");
   puts("#> open map1.mrc");
   puts("#> open map2.mrc");
   puts("#> vop #1 resample onGrid #0");
   puts("#> volume #2 save new.mrc");
   puts("#Chimera will generate the resampled map2.mrc as new.mrc");
   return true;
  }
 }

 fftwf_complex *X1,*Y1,*Z1;
 fftwf_complex *X2,*Y2,*Z2;
 fftwf_complex *X12,*Y12,*Z12;
 float *x1,*y1,*z1;
 float *x2,*y2,*z2;
 float *x12,*y12,*z12;
 fftwf_plan px1,py1,pz1;
 fftwf_plan px2,py2,pz2;
 fftwf_plan px12,py12,pz12;
 int trans[3];

 //Multi-threading
 fftwf_complex **MT_X1,**MT_Y1,**MT_Z1;
 fftwf_complex **MT_X2,**MT_Y2,**MT_Z2;
 fftwf_complex **MT_X12,**MT_Y12,**MT_Z12;
 float **MT_x1,**MT_y1,**MT_z1;
 float **MT_x2,**MT_y2,**MT_z2;
 float **MT_x12,**MT_y12,**MT_z12;
  
 if((x1=(float *)fftwf_malloc(sizeof(float)*xyzdim))==NULL)
  return true;
 if((y1=(float *)fftwf_malloc(sizeof(float)*xyzdim))==NULL)
  return true;
 if((z1=(float *)fftwf_malloc(sizeof(float)*xyzdim))==NULL)
  return true;
 if((X1=(fftwf_complex *)fftwf_malloc(sizeof(fftwf_complex)*xyzdim))==NULL)
  return true;
 if((Y1=(fftwf_complex *)fftwf_malloc(sizeof(fftwf_complex)*xyzdim))==NULL)
  return true;
 if((Z1=(fftwf_complex *)fftwf_malloc(sizeof(fftwf_complex)*xyzdim))==NULL)
  return true;



 if((MT_x2=(float **)fftwf_malloc(sizeof(float *)*Nth))==NULL) return true;
 if((MT_y2=(float **)fftwf_malloc(sizeof(float *)*Nth))==NULL) return true;
 if((MT_z2=(float **)fftwf_malloc(sizeof(float *)*Nth))==NULL) return true;

 if((MT_x12=(float **)fftwf_malloc(sizeof(float *)*Nth))==NULL) return true;
 if((MT_y12=(float **)fftwf_malloc(sizeof(float *)*Nth))==NULL) return true;
 if((MT_z12=(float **)fftwf_malloc(sizeof(float *)*Nth))==NULL) return true;

 if((MT_X2=(fftwf_complex **)fftwf_malloc(sizeof(fftwf_complex *)*Nth))==NULL) return true;
 if((MT_Y2=(fftwf_complex **)fftwf_malloc(sizeof(fftwf_complex *)*Nth))==NULL) return true;
 if((MT_Z2=(fftwf_complex **)fftwf_malloc(sizeof(fftwf_complex *)*Nth))==NULL) return true;
 if((MT_X12=(fftwf_complex **)fftwf_malloc(sizeof(fftwf_complex *)*Nth))==NULL) return true;
 if((MT_Y12=(fftwf_complex **)fftwf_malloc(sizeof(fftwf_complex *)*Nth))==NULL) return true;
 if((MT_Z12=(fftwf_complex **)fftwf_malloc(sizeof(fftwf_complex *)*Nth))==NULL) return true;

 if((MT_mtmp=(MRC*)malloc(sizeof(MRC)*Nth))==NULL) return true;

 //For each thread
 for(int i=0;i<Nth;i++){


  //MRC
  MT_mtmp[i].xdim=m2->xdim;
  MT_mtmp[i].ydim=m2->ydim;
  MT_mtmp[i].zdim=m2->zdim;
  MT_mtmp[i].widthx=m2->widthx;

  MT_mtmp[i].cent[0]=m2->cent[0];
  MT_mtmp[i].cent[1]=m2->cent[1];
  MT_mtmp[i].cent[2]=m2->cent[2];

  MT_mtmp[i].orgxyz[0]=m2->orgxyz[0];
  MT_mtmp[i].orgxyz[1]=m2->orgxyz[1];
  MT_mtmp[i].orgxyz[2]=m2->orgxyz[2];

  if((MT_mtmp[i].dens=(float *)malloc(sizeof(float)*xyzdim))==NULL)
   return true;
  if((MT_mtmp[i].sco=(float *)calloc(sizeof(float),xyzdim))==NULL)
   return true;
  if((MT_mtmp[i].vec=(double **)malloc(sizeof(double *)*xyzdim))==NULL)
   return true;
  for(int j=0;j<xyzdim;j++)
   if((MT_mtmp[i].vec[j]=(double *)malloc(sizeof(double)*3))==NULL)
    return true;

 //Malloc
 if((MT_x2[i]=(float *)fftwf_malloc(sizeof(float)*xyzdim))==NULL)
  return true;
 if((MT_y2[i]=(float *)fftwf_malloc(sizeof(float)*xyzdim))==NULL)
  return true;
 if((MT_z2[i]=(float *)fftwf_malloc(sizeof(float)*xyzdim))==NULL)
  return true;
 if((MT_X2[i]=(fftwf_complex *)fftwf_malloc(sizeof(fftwf_complex)*xyzdim))==NULL)
  return true;
 if((MT_Y2[i]=(fftwf_complex *)fftwf_malloc(sizeof(fftwf_complex)*xyzdim))==NULL)
  return true;
 if((MT_Z2[i]=(fftwf_complex *)fftwf_malloc(sizeof(fftwf_complex)*xyzdim))==NULL)
  return true;

 if((MT_x12[i]=(float *)fftwf_malloc(sizeof(float)*xyzdim))==NULL)
  return true;
 if((MT_y12[i]=(float *)fftwf_malloc(sizeof(float)*xyzdim))==NULL)
  return true;
 if((MT_z12[i]=(float *)fftwf_malloc(sizeof(float)*xyzdim))==NULL)
  return true;
 if((MT_X12[i]=(fftwf_complex *)fftwf_malloc(sizeof(fftwf_complex)*xyzdim))==NULL)
  return true;
 if((MT_Y12[i]=(fftwf_complex *)fftwf_malloc(sizeof(fftwf_complex)*xyzdim))==NULL)
  return true;
 if((MT_Z12[i]=(fftwf_complex *)fftwf_malloc(sizeof(fftwf_complex)*xyzdim))==NULL)
  return true;

 }

 puts("#END malloc");

 //hust copy
 x2=MT_x2[0];
 X2=MT_X2[0];
 x12=MT_x12[0];
 X12=MT_X12[0];

 //input real data in  x1,y1,z1
 for(int i=0;i<xyzdim;i++) x1[i]=m1->vec[i][0];
 for(int i=0;i<xyzdim;i++) y1[i]=m1->vec[i][1];
 for(int i=0;i<xyzdim;i++) z1[i]=m1->vec[i][2];

 for(int i=0;i<xyzdim;i++) x2[i]=m2->vec[i][0];
 //for(int i=0;i<xyzdim;i++) y2[i]=m2->vec[i][1];
 //for(int i=0;i<xyzdim;i++) z2[i]=m2->vec[i][2];

 
 puts("#Make PLANs");
 printf("#Array size= %d * %d * %d = %d\n",n,n,n,n*n*n);
 //FFT PLAN for m1
 int d3=m1->xdim*m1->xdim*m1->xdim;
 double rd3=1.000/(double)d3;

 px1=fftwf_plan_dft_r2c_3d(n,n,n,x1,X1,FFTW_MEASURE); 
 for(int i=0;i<xyzdim;i++) x1[i]=m1->vec[i][0];

 #pragma omp parallel
 #pragma omp sections
 {
  	#pragma omp section
 	{
 	 fftwf_execute_dft_r2c(px1,x1,X1);
 	 for(int i=0;i<xyzdim;i++) X1[i][1]*=-1.000;
	}
	#pragma omp section
 	{
 	 fftwf_execute_dft_r2c(px1,y1,Y1);
 	 for(int i=0;i<xyzdim;i++) Y1[i][1]*=-1.000;
	}
	#pragma omp section
 	{
 	 fftwf_execute_dft_r2c(px1,z1,Z1);
 	 for(int i=0;i<xyzdim;i++) Z1[i][1]*=-1.000;
	}
	#pragma omp section
 	{
	 fftwf_execute_dft_r2c(px1,x2,X2);
	}
 }

 Complex_Complex(X1,X2,X12,xyzdim);
 px12=fftwf_plan_dft_c2r_3d(n,n,n,X12,x12,FFTW_MEASURE); 
 puts("#FIN PLANs");

 double sco;
 TBL *t,*tbl;
 int Ntbl=0;

 if((t=(TBL *)malloc(sizeof(TBL)*TopN))==NULL)
  return true;
 
 puts("#Start Rot");

 //Make job table
 int Njobs=0;
 for(double rx=0;rx<360; rx+=ang){
 for(double ry=0;ry<360; ry+=ang){
 for(double rz=0;rz<=180;rz+=ang){
  Njobs++;
 }}}

 double *jobtbl;
 printf("# NumOfJobs= %d\n",Njobs);
 if((tbl=(TBL *)malloc(sizeof(TBL)*Njobs))==NULL)
  return true;
 int cnt=0;
 for(double rx=0;rx<360; rx+=ang){
 for(double ry=0;ry<360; ry+=ang){
 for(double rz=0;rz<=180;rz+=ang){
  tbl[cnt].r[0]=rx;
  tbl[cnt].r[1]=ry;
  tbl[cnt].r[2]=rz;
  tbl[cnt].sco=0.00;
  cnt++;
 }}}

 Njobs=cnt;
 printf("# NumOfJobs= %d\n",Njobs);
 #pragma omp parallel for schedule(dynamic,5)
 for(int job=0;job<Njobs;job++){
  int th=omp_get_thread_num();
  double rx,ry,rz;
  fftwf_complex *X2,*Y2,*Z2;
  fftwf_complex *X12,*Y12,*Z12;
  float *x2,*y2,*z2;
  float *x12,*y12,*z12;
  double sco;
  int trans[3];
  MRC *mtmp;

  mtmp=&MT_mtmp[th];

  rx=tbl[job].r[0];
  ry=tbl[job].r[1];
  rz=tbl[job].r[2];


  x2=MT_x2[th];
  y2=MT_y2[th];
  z2=MT_z2[th];

  X2=MT_X2[th];
  Y2=MT_Y2[th];
  Z2=MT_Z2[th];

  X12=MT_X12[th];
  Y12=MT_Y12[th];
  Z12=MT_Z12[th];

  x12=MT_x12[th];
  y12=MT_y12[th];
  z12=MT_z12[th];

  RotMRC(m2,mtmp,rx,ry,rz);

  //input
  for(int i=0;i<xyzdim;i++) x2[i]=mtmp->vec[i][0];
  for(int i=0;i<xyzdim;i++) y2[i]=mtmp->vec[i][1];
  for(int i=0;i<xyzdim;i++) z2[i]=mtmp->vec[i][2];


  	 fftwf_execute_dft_r2c(px1,x2,X2);
  	 Complex_Complex(X1,X2,X12,xyzdim);
  	 fftwf_execute_dft_c2r(px12,X12,x12);

  	 fftwf_execute_dft_r2c(px1,y2,Y2);
  	 Complex_Complex(Y1,Y2,Y12,xyzdim);
  	 fftwf_execute_dft_c2r(px12,Y12,y12);

  	 fftwf_execute_dft_r2c(px1,z2,Z2);
  	 Complex_Complex(Z1,Z2,Z12,xyzdim);
  	 fftwf_execute_dft_c2r(px12,Z12,z12);

  //find best
  sco=FindBestTrans(x12,y12,z12,m1->xdim,trans);
  tbl[job].t[0]=trans[0];
  tbl[job].t[1]=trans[1];
  tbl[job].t[2]=trans[2];
  tbl[job].sco=sco*rd3;
  printf("R %.1f %.1f %.1f Best= %d %d %d %f\n",rx,ry,rz,trans[0],trans[1],trans[2],tbl[job].sco);
  

 }

 //STD and Ave
 double Ave=0;
 double Std=0;
 double Sum=0;
 for(int i=0;i<Njobs;i++)
  Sum+=tbl[i].sco;
 Ave=Sum/(double)Njobs;
 Sum=0;
 for(int i=0;i<Njobs;i++)
  Sum+=(tbl[i].sco-Ave)*(tbl[i].sco-Ave);
 Std=sqrt(Sum/(double)Njobs);


 qsort(tbl,Njobs,sizeof(TBL),cmp_tbl);
 //refine topN
 double Bestscore;
 TBL *TopTbl;
 if((TopTbl=(TBL *)malloc(sizeof(TBL)*TopN*27))==NULL)
  return true;

 printf("#refine top %d\n",TopN);

 cnt=0;
 if(ang>5.0){
 for(int i=0;i<TopN;i++){
  double r[3];
  r[0]=tbl[i].r[0];
  r[1]=tbl[i].r[1];
  r[2]=tbl[i].r[2];
 	for(double rx=r[0]-5.0;rx<=r[0]+5.0 && rx<360; rx+=5.0){
 	for(double ry=r[1]-5.0;ry<=r[1]+5.0 && ry<360; ry+=5.0){
 	for(double rz=r[2]-5.0;rz<=r[2]+5.0 && rz<180; rz+=5.0){
 	 TopTbl[cnt].r[0]=rx;
	 if(rx<0)
	  TopTbl[cnt].r[0]+=360.0;
 	 TopTbl[cnt].r[1]=ry;
	 if(ry<0)
	  TopTbl[cnt].r[1]+=360.0;
 	 TopTbl[cnt].r[2]=rz;
	 if(rz<0)
	  TopTbl[cnt].r[2]+=360.0;
 	 TopTbl[cnt].code=360*360*rx+360*ry+rz;
 	 TopTbl[cnt].sco=0.00;
 	 cnt++;
 	}}}
 }
 qsort(TopTbl,cnt,sizeof(TBL),cmp_tbl_code);
 Njobs=0;
 for(int i=0;i<cnt;i++){
  if(i==0){
   Njobs++;
   continue;
  }
  if(TopTbl[i].code!=TopTbl[i-1].code){
   TopTbl[Njobs]=TopTbl[i];
   Njobs++;
  }
 }
 }else{
 	for(int i=0;i<TopN;i++){
		TopTbl[i].r[0]=tbl[i].r[0];
		TopTbl[i].r[1]=tbl[i].r[1];
		TopTbl[i].r[2]=tbl[i].r[2];

		TopTbl[i].t[0]=tbl[i].t[0];
		TopTbl[i].t[1]=tbl[i].t[1];
		TopTbl[i].t[2]=tbl[i].t[2];
		TopTbl[i].sco=tbl[i].sco;
	}
	Njobs=0;//Do not perform refinement
 }
 printf("#REFINE NumOfJobs= %d\n",Njobs);

 #pragma omp parallel for schedule(dynamic,5)
 for(int job=0;job<Njobs;job++){
  int th=omp_get_thread_num();
  double rx,ry,rz;
  fftwf_complex *X2,*Y2,*Z2;
  fftwf_complex *X12,*Y12,*Z12;
  float *x2,*y2,*z2;
  float *x12,*y12,*z12;
  double sco;
  int trans[3];
  MRC *mtmp;

  mtmp=&MT_mtmp[th];

  rx=TopTbl[job].r[0];
  ry=TopTbl[job].r[1];
  rz=TopTbl[job].r[2];


  x2=MT_x2[th];
  y2=MT_y2[th];
  z2=MT_z2[th];

  X2=MT_X2[th];
  Y2=MT_Y2[th];
  Z2=MT_Z2[th];

  X12=MT_X12[th];
  Y12=MT_Y12[th];
  Z12=MT_Z12[th];

  x12=MT_x12[th];
  y12=MT_y12[th];
  z12=MT_z12[th];

  RotMRC(m2,mtmp,rx,ry,rz);

  //input
  for(int i=0;i<xyzdim;i++) x2[i]=mtmp->vec[i][0];
  for(int i=0;i<xyzdim;i++) y2[i]=mtmp->vec[i][1];
  for(int i=0;i<xyzdim;i++) z2[i]=mtmp->vec[i][2];

  	 fftwf_execute_dft_r2c(px1,x2,X2);
  	 Complex_Complex(X1,X2,X12,xyzdim);
  	 fftwf_execute_dft_c2r(px12,X12,x12);

  	 fftwf_execute_dft_r2c(px1,y2,Y2);
  	 Complex_Complex(Y1,Y2,Y12,xyzdim);
  	 fftwf_execute_dft_c2r(px12,Y12,y12);

  	 fftwf_execute_dft_r2c(px1,z2,Z2);
  	 Complex_Complex(Z1,Z2,Z12,xyzdim);
  	 fftwf_execute_dft_c2r(px12,Z12,z12);

  //find best
  sco=FindBestTrans(x12,y12,z12,m1->xdim,trans);
  TopTbl[job].t[0]=trans[0];
  TopTbl[job].t[1]=trans[1];
  TopTbl[job].t[2]=trans[2];
  TopTbl[job].sco=sco*rd3;
  printf("REF %.1f %.1f %.1f Best= %d %d %d %.1f\n",rx,ry,rz,trans[0],trans[1],trans[2],TopTbl[job].sco);
  

 }

 if(ang<=5.0)
	 Njobs=TopN;

 //Show topN
 qsort(TopTbl,Njobs,sizeof(TBL),cmp_tbl);
 for(int i=0;i<TopN;i++){
  PrintTbl(TopTbl,i,m1,m2,Ave,Std);
  RotMRC(m2,&MT_mtmp[0],TopTbl[i].r[0],TopTbl[i].r[1],TopTbl[i].r[2]);
  Bestscore=GetScore(m1,&MT_mtmp[0],TopTbl[i].t);
  if(cmd.ShowGrid)
   ShowVec3(m1,&MT_mtmp[0],TopTbl[i].t);
 }
  //Evaluation
 if(Emode==true){
	 printf("##EVALUATION OF INITIAL POSITION\n");
  int i=0;
  double sco[5];
  TopTbl[i].r[0]=0.00;
  TopTbl[i].r[1]=0.00;
  TopTbl[i].r[2]=0.00;

  TopTbl[i].t[0]=0;
  TopTbl[i].t[1]=0;
  TopTbl[i].t[2]=0;

  RotMRC(m2,&MT_mtmp[0],TopTbl[i].r[0],TopTbl[i].r[1],TopTbl[i].r[2]);
  //Bestscore=GetScore(m1,&MT_mtmp[0],TopTbl[i].t);
  Bestscore=GetScore2(m1,&MT_mtmp[0],TopTbl[i].t,sco);
  TopTbl[i].sco=sco[0];
  PrintTbl(TopTbl,i,m1,m2,Ave,Std);
  if(cmd.ShowGrid)
   ShowVec3(m1,&MT_mtmp[0],TopTbl[i].t);

 }

/*
 fftw_cleanup_threads();

 fftw_destroy_plan(px2);
 fftw_destroy_plan(py2);
 fftw_destroy_plan(pz2);

 fftw_destroy_plan(px12);
 fftw_destroy_plan(py12);
 fftw_destroy_plan(pz12);
*/
 return false;
}

bool SearchMAPfftMT_OVCC(MRC *m1,MRC *m2,double ang, int mode, bool Emode){

 MRC mtmp,*MT_mtmp;
 int n=m2->xdim;
 int xdim=m2->xdim;
 int xydim=m2->xdim*m2->ydim;
 int xyzdim=m2->xdim*m2->ydim*m2->zdim;
 int Nth=omp_get_max_threads();
 int TopN=cmd.TopN;
 if(Emode==true){
	 TopN=0;
  if(m1->cent[0]!=m2->cent[0]
   ||m1->cent[1]!=m2->cent[1]
   ||m1->cent[2]!=m2->cent[2]
   ||m1->orgxyz[0]!=m2->orgxyz[0]
   ||m1->orgxyz[1]!=m2->orgxyz[1]
   ||m1->orgxyz[2]!=m2->orgxyz[2]
   ){
   puts("#For Evaluation Mode");
   puts("#Please use the same coordinate system and map size for map1 and map2.");
   puts("#Example:");
   puts("#In Chimera command line: open map1 and map2 as #0 and #1, then type");
   puts("#> open map1.mrc");
   puts("#> open map2.mrc");
   puts("#> vop #1 resample onGrid #0");
   puts("#> volume #2 save new.mrc");
   puts("#Chimera will generate the resampled map2.mrc as new.mrc");
   return true;
  }
 }

 fftwf_complex *X1,*Y1,*Z1;
 fftwf_complex *X2,*Y2,*Z2;
 fftwf_complex *X12,*Y12,*Z12;
 float *x1,*y1,*z1;
 float *x2,*y2,*z2;
 float *x12,*y12,*z12;
 fftwf_plan px1,py1,pz1;
 fftwf_plan px2,py2,pz2;
 fftwf_plan px12,py12,pz12;
 int trans[3];

 //Multi-threading
 fftwf_complex **MT_X1,**MT_Y1,**MT_Z1;
 fftwf_complex **MT_X2,**MT_Y2,**MT_Z2;
 fftwf_complex **MT_X12,**MT_Y12,**MT_Z12;
 float **MT_x1,**MT_y1,**MT_z1;
 float **MT_x2,**MT_y2,**MT_z2;
 float **MT_x12,**MT_y12,**MT_z12;


 if(mode==0)
  puts("##OVERLAP MODE##");
 if(mode==1)
  puts("##CCC MODE##");
 if(mode==2)
  puts("##PCC MODE##");
 if(mode==3)
  puts("##LAPLACIAN FILTER MODE##");

 double rstd2=1.000/((m1->std*m2->std));
 double rstd3=1.000/((m1->std_norm_ave*m2->std_norm_ave));

 if((x1=(float *)fftwf_malloc(sizeof(float)*xyzdim))==NULL)
  return true;
/*
 if((y1=(float *)fftwf_malloc(sizeof(float)*xyzdim))==NULL)
  return true;
 if((z1=(float *)fftwf_malloc(sizeof(float)*xyzdim))==NULL)
  return true;
*/
 if((X1=(fftwf_complex *)fftwf_malloc(sizeof(fftwf_complex)*xyzdim))==NULL)
  return true;
/*
 if((Y1=(fftwf_complex *)fftwf_malloc(sizeof(fftwf_complex)*xyzdim))==NULL)
  return true;
 if((Z1=(fftwf_complex *)fftwf_malloc(sizeof(fftwf_complex)*xyzdim))==NULL)
  return true;
*/


 if((MT_x2=(float **)fftwf_malloc(sizeof(float *)*Nth))==NULL) return true;
/*
 if((MT_y2=(float **)fftwf_malloc(sizeof(float *)*Nth))==NULL) return true;
 if((MT_z2=(float **)fftwf_malloc(sizeof(float *)*Nth))==NULL) return true;
*/
 if((MT_x12=(float **)fftwf_malloc(sizeof(float *)*Nth))==NULL) return true;

/*
 if((MT_y12=(float **)fftwf_malloc(sizeof(float *)*Nth))==NULL) return true;
 if((MT_z12=(float **)fftwf_malloc(sizeof(float *)*Nth))==NULL) return true;
*/

 if((MT_X2=(fftwf_complex **)fftwf_malloc(sizeof(fftwf_complex *)*Nth))==NULL) return true;
/*
 if((MT_Y2=(fftwf_complex **)fftwf_malloc(sizeof(fftwf_complex *)*Nth))==NULL) return true;
 if((MT_Z2=(fftwf_complex **)fftwf_malloc(sizeof(fftwf_complex *)*Nth))==NULL) return true;
*/
 if((MT_X12=(fftwf_complex **)fftwf_malloc(sizeof(fftwf_complex *)*Nth))==NULL) return true;
/*
 if((MT_Y12=(fftwf_complex **)fftwf_malloc(sizeof(fftwf_complex *)*Nth))==NULL) return true;
 if((MT_Z12=(fftwf_complex **)fftwf_malloc(sizeof(fftwf_complex *)*Nth))==NULL) return true;
*/
 if((MT_mtmp=(MRC*)malloc(sizeof(MRC)*Nth))==NULL) return true;

 //For each thread
 for(int i=0;i<Nth;i++){


  //MRC
  MT_mtmp[i].xdim=m2->xdim;
  MT_mtmp[i].ydim=m2->ydim;
  MT_mtmp[i].zdim=m2->zdim;
  MT_mtmp[i].widthx=m2->widthx;

  MT_mtmp[i].cent[0]=m2->cent[0];
  MT_mtmp[i].cent[1]=m2->cent[1];
  MT_mtmp[i].cent[2]=m2->cent[2];

  MT_mtmp[i].orgxyz[0]=m2->orgxyz[0];
  MT_mtmp[i].orgxyz[1]=m2->orgxyz[1];
  MT_mtmp[i].orgxyz[2]=m2->orgxyz[2];

  if((MT_mtmp[i].dens=(float *)malloc(sizeof(float)*xyzdim))==NULL)
   return true;
  if((MT_mtmp[i].sco=(float *)calloc(sizeof(float),xyzdim))==NULL)
   return true;
  if((MT_mtmp[i].vec=(double **)malloc(sizeof(double *)*xyzdim))==NULL)
   return true;
  for(int j=0;j<xyzdim;j++)
   if((MT_mtmp[i].vec[j]=(double *)malloc(sizeof(double)*3))==NULL)
    return true;

 //Malloc
 if((MT_x2[i]=(float *)fftwf_malloc(sizeof(float)*xyzdim))==NULL)
  return true;
/*
 if((MT_y2[i]=(float *)fftwf_malloc(sizeof(float)*xyzdim))==NULL)
  return true;
 if((MT_z2[i]=(float *)fftwf_malloc(sizeof(float)*xyzdim))==NULL)
  return true;
*/
 if((MT_X2[i]=(fftwf_complex *)fftwf_malloc(sizeof(fftwf_complex)*xyzdim))==NULL)
  return true;
/*
 if((MT_Y2[i]=(fftwf_complex *)fftwf_malloc(sizeof(fftwf_complex)*xyzdim))==NULL)
  return true;
 if((MT_Z2[i]=(fftwf_complex *)fftwf_malloc(sizeof(fftwf_complex)*xyzdim))==NULL)
  return true;
*/

 if((MT_x12[i]=(float *)fftwf_malloc(sizeof(float)*xyzdim))==NULL)
  return true;
/*
 if((MT_y12[i]=(float *)fftwf_malloc(sizeof(float)*xyzdim))==NULL)
  return true;
 if((MT_z12[i]=(float *)fftwf_malloc(sizeof(float)*xyzdim))==NULL)
  return true;
*/
 if((MT_X12[i]=(fftwf_complex *)fftwf_malloc(sizeof(fftwf_complex)*xyzdim))==NULL)
  return true;
/*
 if((MT_Y12[i]=(fftwf_complex *)fftwf_malloc(sizeof(fftwf_complex)*xyzdim))==NULL)
  return true;
 if((MT_Z12[i]=(fftwf_complex *)fftwf_malloc(sizeof(fftwf_complex)*xyzdim))==NULL)
  return true;
*/
 }

 puts("#END malloc");

 //hust copy
 x2=MT_x2[0];
 X2=MT_X2[0];
 x12=MT_x12[0];
 X12=MT_X12[0];

 //input real data in  x1

 if(mode==0){
 for(int i=0;i<xyzdim;i++)
  if(m1->dens[i]>ZERO)
   x1[i]=1.00;
  else
   x1[i]=ZERO;

 for(int i=0;i<xyzdim;i++)
  if(m2->dens[i]>ZERO)
   x2[i]=1.00;
  else
   x2[i]=ZERO;
 }
 if(mode==1){
  for(int i=0;i<xyzdim;i++)
  if(m1->dens[i]>ZERO)
   x1[i]=m1->dens[i];
  else
   x1[i]=ZERO;

 for(int i=0;i<xyzdim;i++)
  if(m2->dens[i]>ZERO)
   x2[i]=m2->dens[i];
  else
   x2[i]=ZERO;
 }
 if(mode==2){ //PCC
  for(int i=0;i<xyzdim;i++)
  if(m1->dens[i]>ZERO)
   x1[i]=m1->dens[i]-m1->ave;
  else
   x1[i]=ZERO;

 for(int i=0;i<xyzdim;i++)
  if(m2->dens[i]>ZERO)
   x2[i]=m2->dens[i]-m2->ave;
  else
   x2[i]=ZERO;
 }
 if(mode==3){ //Laplacian Filter
	 for(int x=0;x<xdim;x++){
	 for(int y=0;y<xdim;y++){
	 for(int z=0;z<xdim;z++){
		 int ind=z*xydim+y*xdim+x;
		 x1[ind]=ZERO;
		 x2[ind]=ZERO;
		 if(m1->dens[ind]>ZERO){
			 //6*d
			 x1[ind]=-6.00*m1->dens[ind];
			 //x+1,y  ,z
			 if(x+1<xdim) x1[ind]+=m1->dens[(z  )*xydim+(y  )*xdim+(x+1)];
			 //x-1,y  ,z
			 if(x-1>=0)   x1[ind]+=m1->dens[(z  )*xydim+(y  )*xdim+(x-1)];
			 //x  ,y+1,z
			 if(y+1<xdim) x1[ind]+=m1->dens[(z  )*xydim+(y+1)*xdim+(x  )];
			 //x  ,y-1,z
			 if(y-1>=0)   x1[ind]+=m1->dens[(z  )*xydim+(y-1)*xdim+(x  )];
			 //x  ,y  ,z+1
			 if(z+1<xdim) x1[ind]+=m1->dens[(z+1)*xydim+(y  )*xdim+(x  )];
			 //x  ,y  ,z-1
			 if(z-1>=0)   x1[ind]+=m1->dens[(z-1)*xydim+(y  )*xdim+(x  )];
		 }
		 if(m2->dens[ind]>ZERO){
			 //6*d
			 x2[ind]=-6.00*m2->dens[ind];
			 //x+1,y  ,z
			 if(x+1<xdim) x2[ind]+=m2->dens[(z  )*xydim+(y  )*xdim+(x+1)];
			 //x-1,y  ,z
			 if(x-1>=0)   x2[ind]+=m2->dens[(z  )*xydim+(y  )*xdim+(x-1)];
			 //x  ,y+1,z
			 if(y+1<xdim) x2[ind]+=m2->dens[(z  )*xydim+(y+1)*xdim+(x  )];
			 //x  ,y-1,z
			 if(y-1>=0)   x2[ind]+=m2->dens[(z  )*xydim+(y-1)*xdim+(x  )];
			 //x  ,y  ,z+1
			 if(z+1<xdim) x2[ind]+=m2->dens[(z+1)*xydim+(y  )*xdim+(x  )];
			 //x  ,y  ,z-1
			 if(z-1>=0)   x2[ind]+=m2->dens[(z-1)*xydim+(y  )*xdim+(x  )];
		 }
	 }}}
 }

 
 puts("#Make PLANs");
 printf("#Array size= %d * %d * %d = %d\n",n,n,n,n*n*n);
 //FFT PLAN for m1
 int d3=m1->xdim*m1->xdim*m1->xdim;
 double rd3=1.000/(double)d3;

 px1=fftwf_plan_dft_r2c_3d(n,n,n,x1,X1,FFTW_MEASURE); 
 //for(int i=0;i<xyzdim;i++) x1[i]=m1->vec[i][0];
 //Overlap
 if(mode==0){
  for(int i=0;i<xyzdim;i++)
   if(m1->dens[i]>ZERO)
    x1[i]=1.00;
   else
    x1[i]=ZERO;
 }
 //CC
 if(mode==1){
  for(int i=0;i<xyzdim;i++)
   if(m1->dens[i]>ZERO)
    x1[i]=m1->dens[i];
   else
    x1[i]=ZERO;
 }
 //PCC
 if(mode==2){
  for(int i=0;i<xyzdim;i++)
   if(m1->dens[i]>ZERO)
    x1[i]=m1->dens[i]-m1->ave;
   else
    x1[i]=ZERO;
 }
 //Laplacian
 if(mode==3){ //Laplacian Filter
	 for(int x=0;x<xdim;x++){
	 for(int y=0;y<xdim;y++){
	 for(int z=0;z<xdim;z++){
		 int ind=z*xydim+y*xdim+x;
		 x1[ind]=ZERO;
		 if(m1->dens[ind]>ZERO){
			 //6*d
			 x1[ind]=-6.00*m1->dens[ind];
			 //x+1,y  ,z
			 if(x+1<xdim) x1[ind]+=m1->dens[(z  )*xydim+(y  )*xdim+(x+1)];
			 //x-1,y  ,z
			 if(x-1>=0)   x1[ind]+=m1->dens[(z  )*xydim+(y  )*xdim+(x-1)];
			 //x  ,y+1,z
			 if(y+1<xdim) x1[ind]+=m1->dens[(z  )*xydim+(y+1)*xdim+(x  )];
			 //x  ,y-1,z
			 if(y-1>=0)   x1[ind]+=m1->dens[(z  )*xydim+(y-1)*xdim+(x  )];
			 //x  ,y  ,z+1
			 if(z+1<xdim) x1[ind]+=m1->dens[(z+1)*xydim+(y  )*xdim+(x  )];
			 //x  ,y  ,z-1
			 if(z-1>=0)   x1[ind]+=m1->dens[(z-1)*xydim+(y  )*xdim+(x  )];
		 }
	 }}}
 }

 #pragma omp parallel
 #pragma omp sections
 {
  	#pragma omp section
 	{
 	 fftwf_execute_dft_r2c(px1,x1,X1);
 	 for(int i=0;i<xyzdim;i++) X1[i][1]*=-1.000;
	}
	#pragma omp section
 	{
	 fftwf_execute_dft_r2c(px1,x2,X2);
	}
 }

 Complex_Complex(X1,X2,X12,xyzdim);
 px12=fftwf_plan_dft_c2r_3d(n,n,n,X12,x12,FFTW_MEASURE); 
 puts("#FIN PLANs");

 //double sco;
 TBL *t,*tbl;
 int Ntbl=0;

 if((t=(TBL *)malloc(sizeof(TBL)*TopN))==NULL)
  return true;
 
 puts("#Start Rot");

 //Make job table
 int Njobs=(360/ang)*(360/ang)*(360/ang);
 double *jobtbl;
 printf("# NumOfJobs= %d\n",Njobs);
 if((tbl=(TBL *)malloc(sizeof(TBL)*Njobs))==NULL)
  return true;
 int cnt=0;
 for(double rx=0;rx<360; rx+=ang){
 for(double ry=0;ry<360; ry+=ang){
 for(double rz=0;rz<=180;rz+=ang){
  tbl[cnt].r[0]=rx;
  tbl[cnt].r[1]=ry;
  tbl[cnt].r[2]=rz;
  tbl[cnt].sco=ZERO;
  cnt++;
 }}}

 Njobs=cnt;
 printf("# NumOfJobs= %d\n",Njobs);
 #pragma omp parallel for schedule(dynamic,5)
 for(int job=0;job<Njobs;job++){
  int th=omp_get_thread_num();
  double rx,ry,rz;
  fftwf_complex *X2,*Y2,*Z2;
  fftwf_complex *X12,*Y12,*Z12;
  float *x2,*y2,*z2;
  float *x12,*y12,*z12;
  double sco;
  int trans[3];
  MRC *mtmp;

  mtmp=&MT_mtmp[th];

  rx=tbl[job].r[0];
  ry=tbl[job].r[1];
  rz=tbl[job].r[2];


  x2=MT_x2[th];
  //y2=MT_y2[th];
  //z2=MT_z2[th];

  X2=MT_X2[th];
  //Y2=MT_Y2[th];
  //Z2=MT_Z2[th];

  X12=MT_X12[th];
  //Y12=MT_Y12[th];
  //Z12=MT_Z12[th];

  x12=MT_x12[th];
  //y12=MT_y12[th];
  //z12=MT_z12[th];

  RotMRC(m2,mtmp,rx,ry,rz);

  //input
  if(mode==0){
   for(int i=0;i<xyzdim;i++) 
    if(mtmp->dens[i]>ZERO)
     x2[i]=1.00;
    else
     x2[i]=ZERO;
  }else if(mode==1){
   for(int i=0;i<xyzdim;i++) 
    if(mtmp->dens[i]>ZERO)
     x2[i]=mtmp->dens[i];
    else
     x2[i]=ZERO;
  }else if(mode==2){
   for(int i=0;i<xyzdim;i++)
    if(mtmp->dens[i]>ZERO)
     x2[i]=mtmp->dens[i]-mtmp->ave;
    else
     x2[i]=ZERO;
  }
  if(mode==3){ //Laplacian Filter
	 for(int x=0;x<xdim;x++){
	 for(int y=0;y<xdim;y++){
	 for(int z=0;z<xdim;z++){
		 int ind=z*xydim+y*xdim+x;
		 x2[ind]=ZERO;
		 if(mtmp->dens[ind]>ZERO){
			 //6*d
			 x2[ind]=-6.00*mtmp->dens[ind];
			 //x+1,y  ,z
			 if(x+1<xdim) x2[ind]+=mtmp->dens[(z  )*xydim+(y  )*xdim+(x+1)];
			 //x-1,y  ,z
			 if(x-1>=0)   x2[ind]+=mtmp->dens[(z  )*xydim+(y  )*xdim+(x-1)];
			 //x  ,y+1,z
			 if(y+1<xdim) x2[ind]+=mtmp->dens[(z  )*xydim+(y+1)*xdim+(x  )];
			 //x  ,y-1,z
			 if(y-1>=0)   x2[ind]+=mtmp->dens[(z  )*xydim+(y-1)*xdim+(x  )];
			 //x  ,y  ,z+1
			 if(z+1<xdim) x2[ind]+=mtmp->dens[(z+1)*xydim+(y  )*xdim+(x  )];
			 //x  ,y  ,z-1
			 if(z-1>=0)   x2[ind]+=mtmp->dens[(z-1)*xydim+(y  )*xdim+(x  )];
		 }
	 }}}
 }
  //for(int i=0;i<xyzdim;i++) x2[i]=mtmp->vec[i][0];
  //for(int i=0;i<xyzdim;i++) y2[i]=mtmp->vec[i][1];
  //for(int i=0;i<xyzdim;i++) z2[i]=mtmp->vec[i][2];

  	 fftwf_execute_dft_r2c(px1,x2,X2);
  	 Complex_Complex(X1,X2,X12,xyzdim);
  	 fftwf_execute_dft_c2r(px12,X12,x12);

/*
  	 fftwf_execute_dft_r2c(px1,y2,Y2);
  	 Complex_Complex(Y1,Y2,Y12,xyzdim);
  	 fftwf_execute_dft_c2r(px12,Y12,y12);

  	 fftwf_execute_dft_r2c(px1,z2,Z2);
  	 Complex_Complex(Z1,Z2,Z12,xyzdim);
  	 fftwf_execute_dft_c2r(px12,Z12,z12);
*/
  //find best
  sco=FindBestTrans1D(x12,m1->xdim,trans);
  if(mode==1)
   sco=sco*rstd2;
  if(mode==2)
   sco=sco*rstd3;
  tbl[job].t[0]=trans[0];
  tbl[job].t[1]=trans[1];
  tbl[job].t[2]=trans[2];
  tbl[job].sco=sco*rd3;
  printf("R %.1f %.1f %.1f Best= %d %d %d %f\n",rx,ry,rz,trans[0],trans[1],trans[2],tbl[job].sco);
 }

 //STD and Ave
 double Ave=0;
 double Std=0;
 double Sum=0;
 for(int i=0;i<Njobs;i++)
  Sum+=tbl[i].sco;
 Ave=Sum/(double)Njobs;
 Sum=0;
 for(int i=0;i<Njobs;i++)
  Sum+=(tbl[i].sco-Ave)*(tbl[i].sco-Ave);
 Std=sqrt(Sum/(double)Njobs);


 qsort(tbl,Njobs,sizeof(TBL),cmp_tbl);
 //refine topN
 double Bestscore;
 TBL *TopTbl;
 if((TopTbl=(TBL *)malloc(sizeof(TBL)*TopN*27))==NULL)
  return true;

 printf("#refine top %d\n",TopN);
 if(ang>5.0){
 cnt=0;
 for(int i=0;i<TopN;i++){
  double r[3];
  r[0]=tbl[i].r[0];
  r[1]=tbl[i].r[1];
  r[2]=tbl[i].r[2];
 	for(double rx=r[0]-5.0;rx<=r[0]+5.0 && rx<360; rx+=5.0){
 	for(double ry=r[1]-5.0;ry<=r[1]+5.0 && ry<360; ry+=5.0){
 	for(double rz=r[2]-5.0;rz<=r[2]+5.0 && rz<180; rz+=5.0){
 	 TopTbl[cnt].r[0]=rx;
	 if(rx<0)
	  TopTbl[cnt].r[0]+=360.0;
 	 TopTbl[cnt].r[1]=ry;
	 if(ry<0)
	  TopTbl[cnt].r[1]+=360.0;
 	 TopTbl[cnt].r[2]=rz;
	 if(rz<0)
	  TopTbl[cnt].r[2]+=360.0;
 	 TopTbl[cnt].code=360*360*rx+360*ry+rz;
 	 TopTbl[cnt].sco=0.00;
 	 cnt++;
 	}}}
 }
 qsort(TopTbl,cnt,sizeof(TBL),cmp_tbl_code);

 Njobs=0;
 for(int i=0;i<cnt;i++){
  if(i==0){
   Njobs++;
   continue;
  }
  if(TopTbl[i].code!=TopTbl[i-1].code){
   TopTbl[Njobs]=TopTbl[i];
   Njobs++;
  }
 }
 }else{
        for(int i=0;i<TopN;i++){
                TopTbl[i].r[0]=tbl[i].r[0];
                TopTbl[i].r[1]=tbl[i].r[1];
                TopTbl[i].r[2]=tbl[i].r[2];

                TopTbl[i].t[0]=tbl[i].t[0];
                TopTbl[i].t[1]=tbl[i].t[1];
                TopTbl[i].t[2]=tbl[i].t[2];
                TopTbl[i].sco=tbl[i].sco;
        }
        Njobs=0;//Do not perform refinement
 }

 printf("#REFINE NumOfJobs= %d\n",Njobs);

 #pragma omp parallel for schedule(dynamic,5)
 for(int job=0;job<Njobs;job++){
  int th=omp_get_thread_num();
  double rx,ry,rz;
  fftwf_complex *X2,*Y2,*Z2;
  fftwf_complex *X12,*Y12,*Z12;
  float *x2,*y2,*z2;
  float *x12,*y12,*z12;
  double sco;
  int trans[3];
  MRC *mtmp;

  mtmp=&MT_mtmp[th];

  rx=TopTbl[job].r[0];
  ry=TopTbl[job].r[1];
  rz=TopTbl[job].r[2];


  x2=MT_x2[th];
  y2=MT_y2[th];
  z2=MT_z2[th];

  X2=MT_X2[th];
  Y2=MT_Y2[th];
  Z2=MT_Z2[th];

  X12=MT_X12[th];
  Y12=MT_Y12[th];
  Z12=MT_Z12[th];

  x12=MT_x12[th];
  y12=MT_y12[th];
  z12=MT_z12[th];

  RotMRC(m2,mtmp,rx,ry,rz);

  //input
  if(mode==0){
   for(int i=0;i<xyzdim;i++) 
    if(mtmp->dens[i]>0.00)
     x2[i]=1.00;
    else
     x2[i]=0.00;
  }
  if(mode==1){
   for(int i=0;i<xyzdim;i++) 
    if(mtmp->dens[i]>0.00)
     x2[i]=mtmp->dens[i];
    else
     x2[i]=0.00;
  }
  if(mode==2){
   for(int i=0;i<xyzdim;i++)
    if(mtmp->dens[i]>0.00)
     x2[i]=mtmp->dens[i]-mtmp->ave;
    else
     x2[i]=0.00;
  }
  if(mode==3){ //Laplacian Filter
	 for(int x=0;x<xdim;x++){
	 for(int y=0;y<xdim;y++){
	 for(int z=0;z<xdim;z++){
		 int ind=z*xydim+y*xdim+x;
		 x2[ind]=ZERO;
		 if(mtmp->dens[ind]>ZERO){
			 //6*d
			 x2[ind]=-6.00*mtmp->dens[ind];
			 //x+1,y  ,z
			 if(x+1<xdim) x2[ind]+=mtmp->dens[(z  )*xydim+(y  )*xdim+(x+1)];
			 //x-1,y  ,z
			 if(x-1>=0)   x2[ind]+=mtmp->dens[(z  )*xydim+(y  )*xdim+(x-1)];
			 //x  ,y+1,z
			 if(y+1<xdim) x2[ind]+=mtmp->dens[(z  )*xydim+(y+1)*xdim+(x  )];
			 //x  ,y-1,z
			 if(y-1>=0)   x2[ind]+=mtmp->dens[(z  )*xydim+(y-1)*xdim+(x  )];
			 //x  ,y  ,z+1
			 if(z+1<xdim) x2[ind]+=mtmp->dens[(z+1)*xydim+(y  )*xdim+(x  )];
			 //x  ,y  ,z-1
			 if(z-1>=0)   x2[ind]+=mtmp->dens[(z-1)*xydim+(y  )*xdim+(x  )];
		 }
	 }}}
 }
  //for(int i=0;i<xyzdim;i++) x2[i]=mtmp->vec[i][0];
  //for(int i=0;i<xyzdim;i++) y2[i]=mtmp->vec[i][1];
  //for(int i=0;i<xyzdim;i++) z2[i]=mtmp->vec[i][2];

  	 fftwf_execute_dft_r2c(px1,x2,X2);
  	 Complex_Complex(X1,X2,X12,xyzdim);
  	 fftwf_execute_dft_c2r(px12,X12,x12);

/*
  	 fftwf_execute_dft_r2c(px1,y2,Y2);
  	 Complex_Complex(Y1,Y2,Y12,xyzdim);
  	 fftwf_execute_dft_c2r(px12,Y12,y12);

  	 fftwf_execute_dft_r2c(px1,z2,Z2);
  	 Complex_Complex(Z1,Z2,Z12,xyzdim);
  	 fftwf_execute_dft_c2r(px12,Z12,z12);
*/
  //find best
  sco=FindBestTrans1D(x12,m1->xdim,trans);
  if(mode==1)
   sco=sco*rstd2;
  if(mode==2)
   sco=sco*rstd3;
  TopTbl[job].t[0]=trans[0];
  TopTbl[job].t[1]=trans[1];
  TopTbl[job].t[2]=trans[2];
  TopTbl[job].sco=sco*rd3;
  printf("REF %.1f %.1f %.1f Best= %d %d %d %f\n",rx,ry,rz,trans[0],trans[1],trans[2],TopTbl[job].sco);
  

 }

 if(ang<=5.0)
	 Njobs=TopN;

 //Show topN
 qsort(TopTbl,Njobs,sizeof(TBL),cmp_tbl);
 for(int i=0;i<TopN;i++){
	 /*
  printf("#%d R={%.1f %.1f %.1f} T={%d %d %d} sco= %f zsco= %f\n",i,
  TopTbl[i].r[0],TopTbl[i].r[1],TopTbl[i].r[2],TopTbl[i].t[0],TopTbl[i].t[1],TopTbl[i].t[2],TopTbl[i].sco,(TopTbl[i].sco-Ave)/Std);
	*/
/*
  int f[3];
	  double t[3];
  f[0]=TopTbl[i].t[0];
  f[1]=TopTbl[i].t[1];
  f[2]=TopTbl[i].t[2];
	 if(TopTbl[i].t[0]>0.5*m2->zdim) f[0]=TopTbl[i].t[0]-m2->xdim;
	 if(TopTbl[i].t[1]>0.5*m2->zdim) f[1]=TopTbl[i].t[1]-m2->xdim;
	 if(TopTbl[i].t[2]>0.5*m2->zdim) f[2]=TopTbl[i].t[2]-m2->xdim;

	 t[0]=m1->orgxyz[0]-m2->orgxyz[0] + f[0]*m2->widthx;
	 t[1]=m1->orgxyz[1]-m2->orgxyz[1] + f[1]*m2->widthx;
	 t[2]=m1->orgxyz[2]-m2->orgxyz[2] + f[2]*m2->widthx;

	 printf("#%d R={%.1f %.1f %.1f} T={%.3f %.3f %.3f} sco= %.1f zsco= %f\n",i,
  	TopTbl[i].r[0],TopTbl[i].r[1],TopTbl[i].r[2],
	t[0],t[1],t[2],
	TopTbl[i].sco,(TopTbl[i].sco-Ave)/Std);
*/
	 PrintTbl(TopTbl,i,m1,m2,Ave,Std);

  RotMRC(m2,&MT_mtmp[0],TopTbl[i].r[0],TopTbl[i].r[1],TopTbl[i].r[2]);
  Bestscore=GetScore(m1,&MT_mtmp[0],TopTbl[i].t);
  if(cmd.ShowGrid)
   ShowVec3(m1,&MT_mtmp[0],TopTbl[i].t);
 }


 //Evaluation
 if(Emode==true){
	 printf("##EVALUATION OF INITIAL POSITION\n");
  int i=0;
  double sco[5];
  TopTbl[i].r[0]=0.00;
  TopTbl[i].r[1]=0.00;
  TopTbl[i].r[2]=0.00;

  TopTbl[i].t[0]=0;
  TopTbl[i].t[1]=0;
  TopTbl[i].t[2]=0;

  RotMRC(m2,&MT_mtmp[0],TopTbl[i].r[0],TopTbl[i].r[1],TopTbl[i].r[2]);
  //Bestscore=GetScore(m1,&MT_mtmp[0],TopTbl[i].t);
  Bestscore=GetScore2(m1,&MT_mtmp[0],TopTbl[i].t,sco);
  TopTbl[i].sco=sco[mode+1];
  PrintTbl(TopTbl,i,m1,m2,Ave,Std);
  if(cmd.ShowGrid)
   ShowVec3(m1,&MT_mtmp[0],TopTbl[i].t);

 }



 return false;
}

void PrintTbl(TBL *TopTbl,int i,MRC *m1,MRC *m2,double Ave,double Std){
 int f[3];
 double t[3],cent[3];
 double mtx[3][3];

 printf("Std= %f Ave= %f\n",Std,Ave);
	  f[0]=TopTbl[i].t[0];
	  f[1]=TopTbl[i].t[1];
	  f[2]=TopTbl[i].t[2];

	  if(TopTbl[i].t[0]>0.5*m2->zdim) f[0]=TopTbl[i].t[0]-m2->xdim;
	  if(TopTbl[i].t[1]>0.5*m2->zdim) f[1]=TopTbl[i].t[1]-m2->xdim;
	  if(TopTbl[i].t[2]>0.5*m2->zdim) f[2]=TopTbl[i].t[2]-m2->xdim;
/*
	 t[0]=m1->orgxyz[0]-m2->orgxyz[0] + f[0]*m2->widthx;
	 t[1]=m1->orgxyz[1]-m2->orgxyz[1] + f[1]*m2->widthx;
	 t[2]=m1->orgxyz[2]-m2->orgxyz[2] + f[2]*m2->widthx;
*/

	 Euler2mtxZYX(RAD(TopTbl[i].r[0]),RAD(TopTbl[i].r[1]),RAD(TopTbl[i].r[2]),mtx);


	 //rotate center pos
	 cent[0]=m2->cent[0]*mtx[0][0]+m2->cent[1]*mtx[0][1]+m2->cent[2]*mtx[0][2];
	 cent[1]=m2->cent[0]*mtx[1][0]+m2->cent[1]*mtx[1][1]+m2->cent[2]*mtx[1][2];
	 cent[2]=m2->cent[0]*mtx[2][0]+m2->cent[1]*mtx[2][1]+m2->cent[2]*mtx[2][2];

	 t[0]=m1->cent[0]-(cent[0]+f[0]*m2->widthx);
	 t[1]=m1->cent[1]-(cent[1]+f[1]*m2->widthx);
	 t[2]=m1->cent[2]-(cent[2]+f[2]*m2->widthx);



	 printf("#%d R={%.1f %.1f %.1f} MTX={%.9f %.9f %.9f %.9f %.9f %.9f %.9f %.9f %.9f} T={%.3f %.3f %.3f} sco= %.3f zsco= %f\n",i,
		TopTbl[i].r[0],TopTbl[i].r[1],TopTbl[i].r[2],
		mtx[0][0],mtx[0][1],mtx[0][2],
		mtx[1][0],mtx[1][1],mtx[1][2],
		mtx[2][0],mtx[2][1],mtx[2][2],
	 t[0],t[1],t[2],
	 TopTbl[i].sco,(TopTbl[i].sco-Ave)/Std);

}




/*
double GetScore(MRC *m1, MRC *m2,int T[3]){
 
 int px,py,pz,t[3];
 int xdim=m1->xdim;
 int xydim=m1->xdim*m1->ydim;
 int xyzdim=m1->xdim*m1->ydim*m1->zdim;
 int ind1,ind2;
 int tot=0;
 int Nm=0;
 double s,sco=0;
 //int hist[20];

 t[0]=T[0];
 t[1]=T[1];
 t[2]=T[2];

 if(t[0]>0.5*xdim) t[0]-=xdim;
 if(t[1]>0.5*xdim) t[1]-=xdim;
 if(t[2]>0.5*xdim) t[2]-=xdim;

 for(int x=0;x<m1->xdim;x++){
  px=x+t[0];
 for(int y=0;y<m1->xdim;y++){
  py=y+t[1];
 for(int z=0;z<m1->xdim;z++){
  pz=z+t[2];
  ind1=xydim*z+xdim*y+x;
  if(m1->dens[ind1]>0.00)
   tot++;

  if(px<0||px>=xdim)
   continue;
  if(py<0||py>=xdim)
   continue;
  if(pz<0||pz>=xdim)
   continue;

  ind2=xydim*pz+xdim*py+px;
  if(m2->dens[ind2]==0.00)
   continue;
  if(m1->dens[ind1]==0.00 && m2->dens[ind2]>0.00){
   tot++;
   continue;
  }

     s=m1->vec[ind1][0]*m2->vec[ind2][0]
     +m1->vec[ind1][1]*m2->vec[ind2][1]
     +m1->vec[ind1][2]*m2->vec[ind2][2];

  m2->sco[ind2]=s;
  sco+=s;
  Nm++;
 }}}
 printf("Overlap= %.4f %d/%d ",(double)Nm/(double)tot,Nm,tot);
 printf("Score= %.1f\n",sco);
 return sco;
}
*/

double GetScore(MRC *m1, MRC *m2,int T[3]){

 int px,py,pz,t[3];
 int xdim=m1->xdim;
 int xydim=m1->xdim*m1->ydim;
 int xyzdim=m1->xdim*m1->ydim*m1->zdim;
 int ind1,ind2;
 int tot=0;
 int Nm=0;
 double s,sco=0;
 double cc_sum=0;
 double pcc_sum=0;
 double std1,std2,d1,d2;
 double pstd1,pstd2,pd1,pd2;
 std1=std2=0;
 int Ncc=0;
 //int hist[20];

 std1=m1->std; std2=m2->std;

 pstd1=m1->std_norm_ave; pstd2=m2->std_norm_ave;

 t[0]=T[0];
 t[1]=T[1];
 t[2]=T[2];

 if(t[0]>0.5*xdim) t[0]-=xdim;
 if(t[1]>0.5*xdim) t[1]-=xdim;
 if(t[2]>0.5*xdim) t[2]-=xdim;

 for(int x=0;x<m1->xdim;x++){
  px=x+t[0];
 for(int y=0;y<m1->xdim;y++){
  py=y+t[1];
 for(int z=0;z<m1->xdim;z++){
  pz=z+t[2];
  ind1=xydim*z+xdim*y+x;
  if(m1->dens[ind1]>0.00)
   tot++;

  if(px<0||px>=xdim)
   continue;
  if(py<0||py>=xdim)
   continue;
  if(pz<0||pz>=xdim)
   continue;

  ind2=xydim*pz+xdim*py+px;


     if(m1->dens[ind1]>0){
      d1=m1->dens[ind1];
      pd1=m1->dens[ind1]-m1->ave;
     }else{
      d1=pd1=0;
     }
     if(m2->dens[ind2]>0){
      d2=m2->dens[ind2];
      pd2=m2->dens[ind2]-m2->ave;
     } else{
      d2=pd2=0;
     }
      cc_sum+=d1*d2;
      pcc_sum+=pd1*pd2;
      //std1+=d1*d1;
      //std2+=d2*d2;

  if(m2->dens[ind2]==0.00)
   continue;
  if(m1->dens[ind1]==0.00 && m2->dens[ind2]>0.00){
   tot++;
   continue;
  }

     s=m1->vec[ind1][0]*m2->vec[ind2][0]
     +m1->vec[ind1][1]*m2->vec[ind2][1]
     +m1->vec[ind1][2]*m2->vec[ind2][2];

  m2->sco[ind2]=s;
  sco+=s;
  Nm++;
  }}}
 //std1=sqrt(std1);
 //std2=sqrt(std2);
 printf("Overlap= %.4f %d/%d CC= %f PCC= %f ",(double)Nm/(double)tot,Nm,tot,cc_sum/(std1*std2),
		 pcc_sum/(pstd1*pstd2));
 printf("Score= %.1f\n",sco);
 return sco;
}

double GetScore2(MRC *m1, MRC *m2,int T[3],double results[5]){

 int px,py,pz,t[3];
 int xdim=m1->xdim;
 int xydim=m1->xdim*m1->ydim;
 int xyzdim=m1->xdim*m1->ydim*m1->zdim;
 int ind1,ind2;
 int tot=0;
 int Nm=0;
 double s,sco=0;
 double cc_sum=0;
 double pcc_sum=0;
 double lap_sum=0;
 double std1,std2,d1,d2;
 double pstd1,pstd2,pd1,pd2;
 double lap1,lap2;
 std1=std2=0;
 int Ncc=0;
 //int hist[20];

 std1=m1->std; std2=m2->std;

 pstd1=m1->std_norm_ave; pstd2=m2->std_norm_ave;

 t[0]=T[0];
 t[1]=T[1];
 t[2]=T[2];

 if(t[0]>0.5*xdim) t[0]-=xdim;
 if(t[1]>0.5*xdim) t[1]-=xdim;
 if(t[2]>0.5*xdim) t[2]-=xdim;

 for(int x=0;x<m1->xdim;x++){
  px=x+t[0];
 for(int y=0;y<m1->xdim;y++){
  py=y+t[1];
 for(int z=0;z<m1->xdim;z++){
  pz=z+t[2];
  ind1=xydim*z+xdim*y+x;
  if(m1->dens[ind1]>0.00)
   tot++;

  if(px<0||px>=xdim)
   continue;
  if(py<0||py>=xdim)
   continue;
  if(pz<0||pz>=xdim)
   continue;

  ind2=xydim*pz+xdim*py+px;


     if(m1->dens[ind1]>0){
      d1=m1->dens[ind1];
      pd1=m1->dens[ind1]-m1->ave;
     }else{
      d1=pd1=0;
     }
     if(m2->dens[ind2]>0){
      d2=m2->dens[ind2];
      pd2=m2->dens[ind2]-m2->ave;
     } else{
      d2=pd2=0;
     }
      cc_sum+=d1*d2;
      pcc_sum+=pd1*pd2;

      	//Laplacian filter
	if(d1>0 && d2> 0){
	lap1=-6.00*d1;
	lap2=-6.00*d2;
	//x+1,y,z
	if(x+1<xdim) 	lap1+=m1->dens[xydim*z +xdim*y +x+1];
	if(px+1<xdim) 	lap2+=m2->dens[xydim*pz+xdim*py+px+1];
	//x-1,y,z
	if(x-1>=0) 	lap1+=m1->dens[xydim*z +xdim*y +x-1];
	if(px-1>=0) 	lap2+=m2->dens[xydim*pz+xdim*py+px-1];
	//x,y+1,z
	if(y+1<xdim) 	lap1+=m1->dens[xydim*z +xdim*(y+1) +x];
	if(py+1<xdim) 	lap2+=m2->dens[xydim*pz+xdim*(py+1)+px];
	//x,y-1,z
	if(y-1>=0) 	lap1+=m1->dens[xydim*z +xdim*(y-1) +x];
	if(py-1>=0) 	lap2+=m2->dens[xydim*pz+xdim*(py-1)+px];
	//x,y,z+1
	if(z+1<xdim) 	lap1+=m1->dens[xydim*(z+1) +xdim*y +x];
	if(pz+1<xdim) 	lap2+=m2->dens[xydim*(pz+1)+xdim*py+px];
	//x,y,z-1
	if(z-1>=0) 	lap1+=m1->dens[xydim*(z-1) +xdim*y +x];
	if(pz-1>=0) 	lap2+=m2->dens[xydim*(pz-1)+xdim*py+px];

	lap_sum+=lap1*lap2;
	}

  //DOT score
  if(m2->dens[ind2]==0.00)
   continue;
  if(m1->dens[ind1]==0.00 && m2->dens[ind2]>0.00){
   tot++;
   continue;
  }

     //dot
     s=m1->vec[ind1][0]*m2->vec[ind2][0]
     +m1->vec[ind1][1]*m2->vec[ind2][1]
     +m1->vec[ind1][2]*m2->vec[ind2][2];

  m2->sco[ind2]=s;
  sco+=s;
  Nm++;
  }}}
 //std1=sqrt(std1);
 //std2=sqrt(std2);
 printf("Overlap= %.4f %d/%d CC= %f PCC= %f LAP= %f ",(double)Nm/(double)tot,Nm,tot,cc_sum/(std1*std2),
		 pcc_sum/(pstd1*pstd2),lap_sum);
 printf("Score= %.1f\n",sco);
 results[0]=sco;//dot score
 results[1]=(double)Nm;//Overlap
 results[2]=cc_sum/(std1*std2);//CC
 results[3]=pcc_sum/(pstd1*pstd2);//PCC
 results[4]=lap_sum;//LAPLACIAN score
 return sco;
}
