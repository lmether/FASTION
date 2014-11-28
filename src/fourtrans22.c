/*

© Copyright 2014 CERN. 

This software is distributed under the terms of the GNU General Public 
Licence version 3 (GPL Version 3), copied verbatim in the file LICENCE.
In applying this licence, CERN does not waive the privileges and immunities 
granted to it by virtue of its status as an Intergovernmental Organization 
or submit itself to any jurisdiction.

Authors: Giovanni Rumolo, Lotta Mether

*/

#include "fftw/fftw.h"

void fourtrans2(double *d,int nn[],int dir)
{
    static first1=1,first2=1;
    static fftwnd_plan p,p2;
    fftw_complex *out=NULL,*in=NULL;

    if (dir==1) {
        if (first1) {
/*
            FILE *f;
            f=fopen("wisdom2.test","r");
            fftw_import_wisdom_from_file(f);
            fclose(f);
*/
            in=(fftw_complex*)alloca(sizeof(fftw_complex)*nn[0]*nn[1]);
            out=(fftw_complex*)alloca(sizeof(fftw_complex)*nn[0]*nn[1]);
/*
            p = fftw2d_create_plan_specific(nn[1],nn[0],FFTW_FORWARD,
                                            FFTW_MEASURE | FFTW_IN_PLACE |
FFTW_USE_WISDOM,
                                            in,1,out,1);
*/
            p = fftw2d_create_plan_specific(nn[1],nn[0],FFTW_FORWARD,
                                            FFTW_ESTIMATE | FFTW_IN_PLACE,
                                            in,1,out,1);

            first1=0;
        }
        fftwnd_one(p,(fftw_complex*)d,out);
     }
    else{
        if (first2) {
            in=(fftw_complex*)alloca(sizeof(fftw_complex)*nn[0]*nn[1]);
            out=(fftw_complex*)alloca(sizeof(fftw_complex)*nn[0]*nn[1]);
/*
            p2 = fftw2d_create_plan_specific(nn[1],nn[0],FFTW_BACKWARD,
                                             FFTW_MEASURE | FFTW_IN_PLACE |
FFTW_USE_WISDOM,
                                             in,1,out,1);
*/

            p2 = fftw2d_create_plan_specific(nn[1],nn[0],FFTW_BACKWARD,
                                             FFTW_ESTIMATE | FFTW_IN_PLACE,
                                             in,1,out,1);

            first2=0;
        }
        fftwnd_one(p2,(fftw_complex*)d,out);
    }
    return;
}

void fourtrans3(double *d,int nn[],int dir)
{
    static first1=1,first2=1;
    static fftw_plan p_a,p_b,p2_a,p2_b;
    fftw_complex *out=NULL,*in=NULL;

    if (dir==1) {
        if (first1) {
/*
            FILE *f;
            f=fopen("wisdom2.test","r");
            fftw_import_wisdom_from_file(f);
            fclose(f);
*/
            in=(fftw_complex*)alloca(sizeof(fftw_complex)*nn[0]*nn[1]);
            out=(fftw_complex*)alloca(sizeof(fftw_complex)*nn[0]*nn[1]);
/*
            p = fftw2d_create_plan_specific(nn[1],nn[0],FFTW_FORWARD,
                                            FFTW_MEASURE | FFTW_IN_PLACE |
FFTW_USE_WISDOM,
                                            in,1,out,1);

            p = fftw2d_create_plan_specific(nn[1],nn[0],FFTW_FORWARD,
                                            FFTW_ESTIMATE | FFTW_IN_PLACE,
                                            in,1,out,1);

*/

            p_a = fftw_create_plan_specific(nn[0],FFTW_FORWARD,
                                            FFTW_ESTIMATE | FFTW_IN_PLACE,
                                            in,1,out,1);
            p_b = fftw_create_plan_specific(nn[1],FFTW_FORWARD,
                                            FFTW_ESTIMATE | FFTW_IN_PLACE,
                                            in,nn[0],out,1);
            first1=0;
        }
//      fftwnd_one(p,(fftw_complex*)d,out);
        fftw(p_b,nn[0]/2,(fftw_complex*)d,nn[0],1,NULL,1,1);
        fftw(p_a,nn[1],(fftw_complex*)d,1,nn[0],NULL,1,1);
     }
    else{
        if (first2) {
            in=(fftw_complex*)alloca(sizeof(fftw_complex)*nn[0]*nn[1]);
            out=(fftw_complex*)alloca(sizeof(fftw_complex)*nn[0]*nn[1]);
/*
            p2 = fftw2d_create_plan_specific(nn[1],nn[0],FFTW_BACKWARD,
                                             FFTW_MEASURE | FFTW_IN_PLACE |
FFTW_USE_WISDOM,
                                             in,1,out,1);

            p2 = fftw2d_create_plan_specific(nn[1],nn[0],FFTW_BACKWARD,
                                             FFTW_ESTIMATE | FFTW_IN_PLACE,
                                             in,1,out,1);

*/

            p2_a = fftw_create_plan_specific(nn[0],FFTW_BACKWARD,
                                             FFTW_ESTIMATE | FFTW_IN_PLACE,
                                             in,1,out,1);
            p2_b = fftw_create_plan_specific(nn[1],FFTW_BACKWARD,
                                             FFTW_ESTIMATE | FFTW_IN_PLACE,
                                             in,nn[0],out,1);
            first2=0;
        }
//      fftwnd_one(p2,(fftw_complex*)d,out);
      fftw(p2_a,nn[1],(fftw_complex*)d,1,nn[0],NULL,1,1);
      fftw(p2_b,nn[0]/2,(fftw_complex*)d,nn[0],1,NULL,1,1);     
    }
    return;
}
