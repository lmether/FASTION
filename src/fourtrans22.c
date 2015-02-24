#include "fftw3.h"

void fourtrans2(double *d,int nn[],int dir, int id)
{
  static int first1[2]={1,1},first2[2]={1,1};
  static fftw_plan p[2],p2[2];
  fftw_complex *in;

  if (dir==1){
    if (first1[id]){
      p[id]=NULL;
      in=(fftw_complex*)d;
      p[id] = fftw_plan_dft_2d(nn[1],nn[0],in,in,-1,FFTW_ESTIMATE);
      first1[id]=0;
    }
    fftw_execute(p[id]);
  }
  else{
    if (first2[id]){
      p2[id]=NULL;
      in=(fftw_complex*)d;
      p2[id] = fftw_plan_dft_2d(nn[1],nn[0],in,in,+1,FFTW_ESTIMATE);
      first2[id]=0;
    }	
    fftw_execute(p2[id]);
  }
  return;
}

void fourtrans3(double *d,int nn[],int dir, int id)
{
  static int first1[2]={1,1},first2[2]={1,1};
  static fftw_plan p_a[2],p_b[2],p2_a[2],p2_b[2];
  fftw_complex *in;

  if (dir==1){
    if (first1[id]){
      in=(fftw_complex*)d;
      p_a[id] = fftw_plan_many_dft(1,&nn[0],nn[1],in,NULL,1,nn[0],in,NULL,1,nn[0], \
			       FFTW_FORWARD,FFTW_ESTIMATE);
      p_b[id] = fftw_plan_many_dft(1,&nn[1],nn[0]/2,in,NULL,nn[0],1,in,NULL,nn[0],1, \
			       FFTW_FORWARD,FFTW_ESTIMATE);
      first1[id]=0;
    }
    fftw_execute(p_b[id]);
    fftw_execute(p_a[id]);
  }
  else{
    if (first2[id]){
      in=(fftw_complex*)d;
      p2_a[id] = fftw_plan_many_dft(1,&nn[0],nn[1],in,NULL,1,nn[0],in,NULL,1,nn[0],	\
				FFTW_BACKWARD,FFTW_ESTIMATE);
      p2_b[id] = fftw_plan_many_dft(1,&nn[1],nn[0]/2,in,NULL,nn[0],1,in,NULL,nn[0],1, \
				FFTW_BACKWARD,FFTW_ESTIMATE);
      first2[id]=0;
    }
    fftw_execute(p2_a[id]);
    fftw_execute(p2_b[id]);     
  }
  return;
}

