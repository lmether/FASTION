/*

© Copyright 2014 CERN. 

This software is distributed under the terms of the GNU General Public 
Licence version 3 (GPL Version 3), copied verbatim in the file LICENCE.
In applying this licence, CERN does not waive the privileges and immunities 
granted to it by virtue of its status as an Intergovernmental Organization 
or submit itself to any jurisdiction.

Authors: Giovanni Rumolo, Lotta Mether

This program is used to study the fast ion instability across a
train of electron bunches which propagates through a lattice.

*/

/* Include all facilities   */

//#include  <gsl_sf.h>
//#include  <gsl_sf_erf.h>
#include  <string.h>
#include  <stdio.h>
//#include  <iostream.h>  
#include  <math.h>
//#include  <libio.h>
#include  <time.h>
#include  <stdlib.h>
#include  <complex.h>  /* with this package one can define complex variables,
                          if needed; for instance, complex<double> z         */

/* Definition of constants  */
/* Physical constants....   */

#define  PI      3.14159265358979
#define  C       2.99792458e8       /* velocity of light          [m/s]      */
#define  E       1.6021892e-19      /* electron charge            [Cb]       */
#define  MU      1.66055e-27        /* proton mass                [kg]       */
#define  ME      9.109534e-31       /* proton-positron mass       [kg]       */
#define  EPS0    8.8542e-12         /* dielectric constant        [C^2/Jm]   */
#define  RP      1.5458775e-18      /* classical proton radius    [m]        */
#define  RE      2.817938e-15       /* classical electron radius  [m]        */
#define  I0      31.356e6           /* ... for the tune shift..   [A]        */
#define  ELMS    9.109534e-31       /* electron mass              [Kg]       */
#define  Z0      376.99112          /* Vacuum impedance           [Ohm]      */

/* Parameters for the computation... */
				

#define  NION    500               /* Number of macro-ions produced per
                                     bunch passage                      */
//#define  NION    1                  /* Number of macro-ions produced per
//                                       bunch passage                      */
#define  NELB    10000              /* Number of macroelectrons
                                       in one bunch                       */
#define  NBUN    312                /* Number of bunches in the train     */
#define  NDC     128                /* Number of cells to produce a
                                       distribution function              */
                                      

#define  NR_END    1
#define  FREE_ARG  char*

/* type definition  */

typedef struct /* GRID */
{
  int n_cell_x,n_cell_y;
  double delta_x,delta_y;
  double offset_x,offset_y;
  double min_x,max_x,min_y,max_y;
  double cut_x,cut_y;
  double *rho1,*rho2,*phi1,*phi2,*dist,*temp,*pot;
  int integration_method,distribution,interpolation;
} GRID;

/* Definition of the external variables   */

FILE  *ions_pr               ,
      *centr_pr              ,
      *trainhdtl_pr          ,
      *ele_phase             ,
      *ion_inphase           ,
      *ion_inix[15]          ,
      *ion_iniy[15]          ,
      *samp_ion[15]          ,
      *samp2_ion[15]         ,
      *sigma_pr              ,
      *wakefield_pr          ,
      *wakefield2_pr         ,
      *grid_pr               ;

char  filename[80]           ; 

long  mmain                 ,   /* down to mmain -> indices used in main for loops  */
      imain                 ,
      jmain                 ,
      kmain                 ,
      lmain                 ,
      mion[NION*NBUN]       ,   /* contains the mass index of the ion               */

      it                    ,   /* index used for the main cycle over the total 
                                   number of interactions: 2*number of fodos        */
      nt                    ,   /* index used for cycle over the number of turns    */
      iion                  ,   /* switch for fast ion instability                  */
      ifi                   ,   /* switch for field ionization                      */
      i_lattice             ,   /* switch for use of external Twiss file            */
      nspe                  ,   /* number of ion species                            */
      nstep                 ,   /* number of interactions                           */
      nturn                 ,   /* number of turns through lattice                  */
      nfodo                 ,   /* number of fodo cells through the line            */
      ili                   ,   /* flag for the function ffrank: it has to be =! 1  */
      i_cool                ,   /* flag for the distribution of the emittance
                                   cooling due to acceleration                      */
      iwake_flag            ,   /* flag for the wake field: input parameter         */
      i_pipe                ,   /* flag for gemoetry in effective wake calculation
                                   0 -> circular (a=b)  1 -> flat  (a>>b)                        */
      i_dip                 ,
      ionoutgrid            ,   /* counter of ions out of the grid                               */
      ntemp                 ,   /* temporary number of ions while bunches are passing            */
      i_kick                ,   /* switch on the type of initial kick given to the train
                                   0 -> no kick
                                   1 -> constant kick along the train (=nkick0*sigma0)
                                   2 -> sinusoidal kick along the train (amplitude nkick0*sigma0
                                        and and wave number h_pert)
				   3 -> random kick with maximum amplitude nkick0*sigma0         */
      h_pertx               ,   /* Harmonic number of the horizontal kick in t=0                 */
      h_perty               ,   /* Harmonic number of the vertical kick in t=0                   */
      n_diag                ,   /* every how many turns the bunch shape is stored                */
      n_diag2               ;   /* every how many turns the ion distributions are stored         */

double xel[NBUN*NELB]       ,   /* horizontal coordinate for electrons              */
       yel[NBUN*NELB]       ,   /* vertical coordinate for electrons                */
       xpel[NBUN*NELB]      ,   /* dx/ds  for electrons                             */
       ypel[NBUN*NELB]      ,   /* dy/ds  for electrons                             */
       xion[NBUN*NION]      ,   /* horizontal coordinate for ions                   */
       yion[NBUN*NION]      ,   /* vertical coordinate for ions                     */
       xpion[NBUN*NION]     ,   /* dxe/dt  for ions                                 */
       ypion[NBUN*NION]     ,   /* dye/dt  for ions                                 */
       qionfi[NBUN]         ,   /* number of ions in one macroion with field ionization     */
       xs[NBUN]             ,
       ys[NBUN]             ,
       sxs[NBUN]            ,
       sys[NBUN]            ,
       sxp2s[NBUN]          ,
       syp2s[NBUN]          ,
       xppr[NBUN]           ,
       yppr[NBUN]           ,
       xel_max[NBUN]        ,   /* [m] maximum electron x position                  */
       yel_max[NBUN]        ,   /* [m] minimum electron y position                  */
       xsion                ,
       ysion                ,
       sxsion               ,
       sysion               ,
       sx0fion              ,   /* horizontal semi-axis of section fully ionized by field ionization         */
       sy0fion              ,   /* verical semi-axis of section fully ionized by field ionization            */
       drrefill             ,   /* radius re-filled between bunch passages                                   */
       vprelx[NELB]         ,
       vprely[NELB]         ,
       vprelxp[NELB]        ,
       vprelyp[NELB]        ,
       vprionx[NION]        ,
       vpriony[NION]        ,
       vprionr[NION]        ,
  //indx[NDC]            ,   /* grid points for the horizontal distribution of electrons        */
  //hde[NDC]             ,   /* horizontal distribution of electrons after having been kicked   */
  //indy[NDC]            ,   /* grid points for the vertical distribution of electrons          */
  //vde[NDC]             ,   /* vertical distribution of electrons after having been kicked     */

       indxp[NDC]           ,   /* grid points for the bunch long. distribution     */
       hdp[NDC]             ,   /* Bunch longitudinal distribution                  */
       indxp2[NDC]          ,   /* grid points for the bunch long. distribution     */
       hdp2[NDC]            ,   /* Bunch longitudinal distribution                  */

  //indyp[NDC]           ,   /* grid points for the vertical distribution of the bunch particles  */
  //vdp[NDC]             ,   /* vertical distribution of bunch particles before interacting with the particles */
  //indxpp[NDC]          ,   /* grid points for the distribution of horizontal momenta of the bunch particles  */
  //hpdp[NDC]            ,   /* distribution of horizontal momenta of the bunch particles before interacting with the electrons  */
  //indypp[NDC]          ,   /* grid points for the distribution of vertical momenta of the bunch particles    */
  //vpdp[NDC]            ,   /* distribution of vertical momenta of the bunch particles before interacting with the electrons    */

  /* now the declaration of the input parameters starts...                          */

       *pss                 ,   /* [nTorr] Partial pressures of gas components 
                                   in the vacuum chamber                            */
       *am                  ,   /* Mass numbers of rest gas molecules               */
       *crsec               ,   /* [MBarn] Ionization cross section for the different
                                   molecule types of rest gas                       */
       *prob                ,   /* Probability distribution for different types of ions (gas ionization)             */
       *prob2               ,   /* Probability distribution for different types of ions (field ionization)           */
       factor1              ,   /* auxiliary parameter for the force calculation    */
       *fscale1             ,   /* auxiliary parameter for the force calculation    */
       nele                 ,   /* number of electrons per bunch                    */
       nionpb               ,   /* number of ions per bunch                         */
       nionpb2              ,   /* number of ions per bunch                         */
       ptot                 ,   /* total pressure in the vacuum chamber [nTorr]     */
       bsp                  ,   /* [ns] bunch spacing                               */
       sz0                  ,   /* [ns] rms bunch length                            */
       emxN0                ,   /* [nm] rms hor. normalized emittance               */
       emyN0                ,   /* [nm] rms vert. normalized emittance              */
       gammain              ,   /* relativistic gamma at beginning of line          */
       gammafin             ,   /* relativistic gamma at end of line                */
       gammanext            ,   /* stores a local value of gamma along the line     */
       gammaprev            ,   /* stores a local value of gamma along the line     */
       fodol                ,   /* [m] length of the FODO cell                      */
       mu                   ,   /* [degrees] phase advance per cell                 */
       nxkick0              ,   /* Amplitude of the horizontal kick in t=0 in number of sigmas   */
       nykick0              ,   /* Amplitude of the vertical kick in t=0 in number of sigmas     */
       rpipey               ,    /* y-size of the chamber for resistive wall        */
       elec_thresh          ,    /* threshold value of electric field to have 
                                    field ionization                                */
       freqr                ,    /* resonant frequency of the broad band resonator
                                    to be expressed in GHz (x and y assumed equal)  */
       merit                ,    /* quality factor of the broad band impedance      */
       zetat                ,    /* shunt impedance of the braod band resonator
                                    to be expressed in MOhm/m                       */
       freqrz               ,    /* resonant frequency of the longitudinal resonator
                                    to be expressed in MHz                          */
       meritz               ,    /* quality factor of the longitudinal impedance    */
       rshz                 ,    /* shunt impedance of the braod band resonator
                                    to be expressed in MOhm                         */
       omegar               ,    /* auxiliary parameters for transverse wake        */
       omegabar             , 
       alphat                ,
       omegarz              ,    /* auxiliary parameters for longitudinal wake      */
       omegabarz            , 
       alphaz               ,
       beta_max             ,   /* [m] maximum beta function in the FODO            */
       beta_min             ,   /* [m] minimum beta function in the FODO            */
       alpha_max            ,   /* [m] maximum alfa function in the FODO            */
       alpha_min            ,   /* [m] minimum alfa function in the FODO            */
       sx_min               ,   /* [m] maximum rms hor. beam size                   */
       sx_max               ,   /* [m] minimum rms hor. beam size                  */
       sy_min               ,   /* [m] maximum rms vert. beam size                   */
       sy_max               ,   /* [m] minimum rms vert. beam size                  */
       xion_max             ,   /* [m] maximum ion x position                       */
       yion_max             ,   /* [m] minimum ion y position                       */
       factor2              ,   /* auxiliary parameter for the force calculation    */
       fscale2              ,   /* auxiliary parameter for the force calculation    */

       eel                  ,   /* auxiliary parameter for the force calculation    */
       epr                  ,   /* auxiliary parameter for the force calculation    */

       wakefac              ,   /* factor for the wake field kick                   */
       tstepp               ,   /* time step between two interactions               */
       tstep                ,   /* time step during the interaction bunch-cloud     */
       dstep                ,   /* s-distance between two kick points               */
       qe0                  ,   /* number of electrons in one macroelectron         */
       qion                 ,   /* number of ions in one macroion                   */
       cxya                 ,   /* cxya -> sxya: transport parameters for the phase */
       sxya                 ,   /* space coordinates of the bunch particles         */
       cxa                  ,   /* cxya -> sxya: transport parameters for the phase */
       sxa                  ,   /* space coordinates of the bunch particles         */
       cya                  ,   /* cxya -> sxya: transport parameters for the phase */
       sya                  ,   /* space coordinates of the bunch particles         */

       xave                 ,   /* average horizontal offset of the bunch           */
       yave                 ,   /* average vertical offset of the bunch             */
       zave                 ,   /* average longitudinal offset of the bunch         */
       xpave                ,   /* average bunch offset in dx/ds                    */
       ypave                ,   /* average bunch offset in dy/ds                    */
       sxsize               ,   /* square of horizontal rms bunch size              */
       sysize               ,   /* square of vertical rms bunch size                */
       sxpsize              ,
       sypsize              ,
       xemit                ,   /* horizontal emittance                             */
       yemit                ,   /* vertical emittance                               */
  //invarx               ,   /* Courant-Snyder horizontal invariant              */
  //invary               ,   /* Courant-Snyder vertical invariant                */

       xoff                 ,   /* local variables for bunch offsets and sizes      */
       yoff                 ,
       sx                   ,
       sy                   , 
       sxp2                 ,
       syp2                 ,
       sxppr                ,
       syppr                ,
       localemitx           ,
       localemity           ,   /* --------------------------------------------     */
       sx0                  ,   /* all variables to be used for transport when the  */
       sy0                  ,   /* beta's are loaded from an external twiss file    */
       sx1                  , 
       sy1                  , 
       betax0               , 
       betay0               , 
       betax1               , 
       betay1               ,
       mux0                 , 
       muy0                 , 
       alfax0               , 
       alfay0               , 
       alfax1               , 
       alfay1               ,
       ss0                  , 
       ss1                  ,   /* --------------------------------------------     */ 
       benergy0             ,   /* --------------------------------------------     */ 
       bsol                 ,
       bgrad                , 
       bdip                 ,
       lengthfi             ,   /* starting point for field ionization              */
       kick_xi              ,   /* velocity x-kick from each bunch on each ion      */ 
       kick_yi              ,   /* velocity y-kick from each slice on each ion      */
       kick_xe              ,   /* x-kick on each bunch electron                    */
       kick_ye              ,   /* y-kick on each bunch electron                    */
       xel0                 ,   /* x-coordinate of a bunch particle after having 
                                   been kicked, before transport to the next kick   */
       xpel0                ,   /* dx/ds of a bunch particle after having been 
                                   kicked and before transport to the next kick     */
       yel0                 ,   /* y-coordinate of a bunch particle after having 
                                   been kicked, before transport to the next kick   */
       ypel0                ,   /* dy/ds of a bunch particle after having been 
                                   kicked and before transport                      */
       coeff1               ,    /* Auxiliary parameters for transport of electrons */
       coeff2               ,    /* through a FODO (or line) */
       coeff3               ,    /* Auxiliary parameters for transport of electrons */
       coeff4               ,    /* through a FODO (or line) */
       coeffga              ,    /* Auxiliary parameters for transport of electrons */
       coeff1x              ,    /* through a FODO (or line) */
       coeff1xbis           ,    /* Auxiliary parameters for transport of electrons */
       coeff2x              ,    /* through a FODO (or line) */
       coeff3x              ,    /* Auxiliary parameters for transport of electrons */
       coeff3xbis           ,    /* Auxiliary parameters for transport of electrons */
       coeff4x              ,    /* through a FODO (or line) */
       coeff4xbis           ,    /* through a FODO (or line) */
       coeff4xter           ,    /* through a FODO (or line) */
       coeff1y              ,    /* Auxiliary parameters for transport of electrons */
       coeff1ybis           ,    /* Auxiliary parameters for transport of electrons */
       coeff2y              ,    /* through a FODO (or line) */
       coeff3y              ,    /* Auxiliary parameters for transport of electrons */
       coeff3ybis           ,    /* Auxiliary parameters for transport of electrons */
       coeff4y              ,    /* through a FODO (or line) */
       coeff4ybis           ,    /* through a FODO (or line) */
       coeff4yter           ,    /* through a FODO (or line) */
       *long_pos            ,    /* From here to the end: vectors to load energies, */
       *benergy             ,    /* positions, beta, alfa functions along a generic */
       *betarrayx           ,    /* line (twiss.dat file used)                      */
       *betarrayy           ,
       *muarrayy            ,
       *muarrayx            ,   
       *alfarrayx           ,   
       *alfarrayy           ;

           
/* For sake of clarity we specify here that the broad band model that we've assumed 
   considers an impedance (omegar = 2*PI*freqr*10^(-9)): 

      Z = omegar/omega * zetat/[1+i merit (omegar/omega - omega/omegar)].
   
   If the impedance is: 

      Z = c/omega * R_s/[1+i merit (omegar/omega - omega/omegar)], 
   
   then one has to convert R_s into our input parameter by means of: 

   zetat=c/omegar*R_s  
   
   zetat [Ohm/m] and R_s [Ohm/m^2]                                                    

   In the longitudinal plane:

      Z|| = rshz/[1+i meritz (omegarz/omega - omega/omegarz)] 
                                                                                     */
/*  Global declarations... 
    function that evaluates the electric mutual force between bunch particles
    and electrons has to be included                                                 */

extern void ffrank (double*, double*, double*, double*, double*, double*, 
                    double*, double*, long*);

/* begin   gio-May04             */

void reverse(char s[])
{
  int c, i, l;

  for (i=0, l=strlen(s)-1; i<l; i++, l--) {
    c = s[i];
    s[i] = s[l];
    s[l] = c;
  }
}


void itoa(int n, char s[])
{
  int i, sign;

  if ((sign = n) < 0)
    n = -n;
  i = 0;
  do { s[i++] =n %10 + '0';
  }  while ((n /= 10) > 0);
  if (sign <0)
    s[i++] = '-';
  s[i] = '\0';
      reverse(s);
}
void nrerror(char error_text[])

{
  fprintf(stderr,"ERROR...\n");
  fprintf(stderr,"%s\n",error_text);
  fprintf(stderr,"...now exiting...\n");
  exit(1);
}

double *dvector(long n1, long nh)
     /* allocate a double vector with subscript range v[n1..nh] */
{
  double *v;
  
  v=(double *)malloc((size_t) ((nh-n1+1+NR_END)*sizeof(double)));
  if (!v) nrerror("allocation failure in dvector()");
  return v-n1+NR_END;
}

void free_dvector(double *v, long n1, long nh)
     /* free a double vector allocated with dvector() */
{
  free((FREE_ARG) (v+n1-NR_END));
}


long *ivector(long n1, long nh)
     /* allocate an integer vector with subscript range v[n1..nh] */
{
  long *v;
  
  v=(long *)malloc((size_t) ((nh-n1+1+NR_END)*sizeof(long)));
  if (!v) nrerror("allocation failure in ivector()");
  return v-n1+NR_END;
}

void free_ivector(long *v, long n1, long nh)
     /* free a long vector allocated with dvector() */
{
  free((FREE_ARG) (v+n1-NR_END));
}

/******************************************************************************
 *** Subroutine     : read_data                                             ***
 *** Effect         : Define main parameters by reading data from specified ***
 ***                  configuration-file                                    ***
 *** Parameters     : none                                                  ***
 *** Gbl var   used :                                                       ***
 *** Gbl var effect :                                                       ***
 *** Constants used : none                                                  ***
 *** Subrout.  used : none                                                  ***
 ******************************************************************************/

read_data ()

/* gio-May04; reads different beta  */
{ char  cfg_filename[100], dummy_string[64], twiss_filename[10];
  FILE  *cfg_file_ptr, *twiss_file;
  long  dummy_int, k_cou, condition, iread;
  double dummy[8];

  printf ("\n\n    Please specify name (without extension) of desired") ;
  printf ("\n    configuration-file : ") ;
  scanf ("%s",filename) ;

  strcpy (cfg_filename,filename);
  strcat (cfg_filename,".cfg");

  cfg_file_ptr = fopen(cfg_filename,"r");
  if (cfg_file_ptr == NULL)
     { printf ("\n\n    Configuration-file does not exist.\n") ;

     //       end () ; }
       exit(2); }

  fscanf (cfg_file_ptr,"%s",dummy_string);
  fscanf (cfg_file_ptr,"%ld",&iion);
  fscanf (cfg_file_ptr,"%s",dummy_string);
  fscanf (cfg_file_ptr,"%ld",&nspe);
  printf ("number of species = %ld\n",nspe);
  pss=dvector(0,nspe-1);
  am=dvector(0,nspe-1);
  crsec=dvector(0,nspe-1);
  fscale1=dvector(0,nspe-1);
  prob=dvector(0,nspe-1);
  prob2=dvector(0,nspe-1);
  fscanf (cfg_file_ptr,"%s",dummy_string);
  for(iread=0;iread<nspe;iread++)
    fscanf (cfg_file_ptr,"%lf",&pss[iread]);
  fscanf (cfg_file_ptr,"%s",dummy_string);
  for(iread=0;iread<nspe;iread++)
    fscanf (cfg_file_ptr,"%lf",&am[iread]);
  fscanf (cfg_file_ptr,"%s",dummy_string);
  for(iread=0;iread<nspe;iread++)
    fscanf (cfg_file_ptr,"%lf",&crsec[iread]);
  fscanf (cfg_file_ptr,"%s",dummy_string);
  fscanf (cfg_file_ptr,"%lf",&nele);
  printf ("number of electrons per bunch = %lg\n",nele);
  fscanf (cfg_file_ptr,"%s",dummy_string);
  fscanf (cfg_file_ptr,"%lf",&bsp);
  fscanf (cfg_file_ptr,"%s",dummy_string);
  fscanf (cfg_file_ptr,"%lf",&emxN0);
  fscanf (cfg_file_ptr,"%s",dummy_string);
  fscanf (cfg_file_ptr,"%lf",&emyN0);
  fscanf (cfg_file_ptr,"%s",dummy_string);
  fscanf (cfg_file_ptr,"%lf",&sz0);
  fscanf (cfg_file_ptr,"%s",dummy_string);
  fscanf (cfg_file_ptr,"%lf",&gammain);
  fscanf (cfg_file_ptr,"%s",dummy_string);  
  fscanf (cfg_file_ptr,"%lf",&gammafin);
  fscanf (cfg_file_ptr,"%s",dummy_string);  
  fscanf (cfg_file_ptr,"%lf",&lengthfi);
  fscanf (cfg_file_ptr,"%s",dummy_string);
  fscanf (cfg_file_ptr,"%ld",&i_lattice);
  printf ("i_lattice=%ld\n",i_lattice);
  fscanf (cfg_file_ptr,"%s",dummy_string);
  fscanf (cfg_file_ptr,"%ld",&nfodo);
  fscanf (cfg_file_ptr,"%s",dummy_string);
  fscanf (cfg_file_ptr,"%lf",&fodol);
  fscanf (cfg_file_ptr,"%s",dummy_string);
  fscanf (cfg_file_ptr,"%lf",&mu);
  fscanf (cfg_file_ptr,"%s",dummy_string);
  fscanf (cfg_file_ptr,"%ld",&i_kick);
  printf ("i_kick=%ld\n",i_kick);
  fscanf (cfg_file_ptr,"%s",dummy_string);
  fscanf (cfg_file_ptr,"%lf",&nxkick0);
  fscanf (cfg_file_ptr,"%s",dummy_string);
  fscanf (cfg_file_ptr,"%lf",&nykick0);
  fscanf (cfg_file_ptr,"%s",dummy_string);
  fscanf (cfg_file_ptr,"%ld",&h_pertx);
  fscanf (cfg_file_ptr,"%s",dummy_string);
  fscanf (cfg_file_ptr,"%ld",&h_perty);
  fscanf (cfg_file_ptr,"%s",dummy_string);
  fscanf (cfg_file_ptr,"%lf",&rpipey);
  printf ("rpipey=%lf\n",rpipey);
  //the following two lines added today 8 Jan 2009 for field ionization
  fscanf (cfg_file_ptr,"%s",dummy_string);
  fscanf (cfg_file_ptr,"%lf",&elec_thresh);
  printf ("elec_thresh=%lf\n",elec_thresh);
  //----------------------------------------------
  fscanf (cfg_file_ptr,"%s",dummy_string);
  fscanf (cfg_file_ptr,"%ld",&iwake_flag);
  //printf ("iwake_flag=%ld\n",iwake_flag);
  fscanf (cfg_file_ptr,"%s",dummy_string);
  fscanf (cfg_file_ptr,"%ld",&i_pipe);
  fscanf (cfg_file_ptr,"%s",dummy_string);
  fscanf (cfg_file_ptr,"%lf",&freqr);
  printf ("freqr=%lf\n",freqr);
  fscanf (cfg_file_ptr,"%s",dummy_string);
  fscanf (cfg_file_ptr,"%lf",&merit);
  fscanf (cfg_file_ptr,"%s",dummy_string);
  fscanf (cfg_file_ptr,"%lf",&zetat);
  fscanf (cfg_file_ptr,"%s",dummy_string);
  fscanf (cfg_file_ptr,"%lf",&freqrz);
  fscanf (cfg_file_ptr,"%s",dummy_string);
  fscanf (cfg_file_ptr,"%lf",&meritz);
  fscanf (cfg_file_ptr,"%s",dummy_string);
  fscanf (cfg_file_ptr,"%lf",&rshz);
  fscanf (cfg_file_ptr,"%s",dummy_string);
  fscanf (cfg_file_ptr,"%ld",&n_diag);
  fscanf (cfg_file_ptr,"%s",dummy_string);
  fscanf (cfg_file_ptr,"%ld",&nturn);
  printf ("n_diag=%ld\n",n_diag);

  fclose(cfg_file_ptr);
  printf ("\n\n\t\tConfiguration-file loaded successfully.\n") ;

  ifi = 0;

  if(i_lattice==1) {

    strcpy (twiss_filename,"twiss.dat");
    //strcat (twiss_filename,".1");
    
    twiss_file = fopen(twiss_filename,"r");
    
    if (twiss_file == NULL)
      { printf ("\n\n    File with Twiss parameters does not exist.\n") ;
      exit(2); }
    
    k_cou = 0;
    
    do {condition = fscanf (twiss_file, "%ld  %lg  %lg  %lg  %lg  %lg  %lg  %lg  %lg\n",
			    &dummy_int,&dummy[0],&dummy[1],&dummy[2],&dummy[3],&dummy[4],&dummy[5],&dummy[6],
			    &dummy[7]);
    ++k_cou;	 
    }
    while (condition != EOF);

    fclose(twiss_file);
    
    nstep=k_cou-1; 
    long_pos=dvector(0,nstep-1);
    benergy=dvector(0,nstep-1);
    betarrayx=dvector(0,nstep-1);
    betarrayy=dvector(0,nstep-1);
    alfarrayx=dvector(0,nstep-1);
    alfarrayy=dvector(0,nstep-1);
    muarrayx=dvector(0,nstep-1);
    muarrayy=dvector(0,nstep-1);
      
    twiss_file = fopen(twiss_filename,"r");
    k_cou = 0;
    
    while (k_cou < nstep) {
      fscanf (twiss_file, "%ld  %lg  %lg  %lg  %lg  %lg  %lg  %lg  %lg\n",
	      &dummy_int,&long_pos[k_cou],&benergy[k_cou],&betarrayx[k_cou],&betarrayy[k_cou],
	      &alfarrayx[k_cou],&alfarrayy[k_cou],&muarrayx[k_cou],&muarrayy[k_cou]);
      printf ("\n %ld  %f %f %f %f %f %f %f %f\n", k_cou, long_pos[k_cou], benergy[k_cou], betarrayx[k_cou],
	      betarrayy[k_cou], alfarrayx[k_cou],alfarrayy[k_cou], muarrayx[k_cou],muarrayy[k_cou]);
      muarrayx[k_cou] *= 2.*PI;
      muarrayy[k_cou] *= 2.*PI;
      ++k_cou;	 
    }
      
    betax0 = betarrayx[0];
    betay0 = betarrayy[0];
    betax1 = betarrayx[1];
    betay1 = betarrayy[1];
    alfax0 = alfarrayx[0];
    alfay0 = alfarrayy[0];
    alfax1 = alfarrayx[1];
    alfay1 = alfarrayy[1];
    mux0 = muarrayx[1];
    muy0 = muarrayy[1];
    ss0 = long_pos[0];
    ss1 = long_pos[1];
    benergy0 = benergy[0];

    fclose(twiss_file);
  }

}

/******************************************************************************
 *** Subroutine     : init_values                                           ***
 *** Effect         : Initialize all parameter values, necessary for        ***
 ***                  calculation.                                          ***
 *** Parameters     : none                                                  ***
 *** Gbl var   used :                                                       ***
 *** Gbl var effect :                                                       ***
 *** Constants used :                                                       ***
 *** Subrout.  used : none                                                  ***
 ******************************************************************************/

init_values ()

{ int im;
 double accum;

 if(i_lattice==0)
   dstep = fodol/2.;
 else
   dstep = ss1-ss0;

 tstepp = dstep/C;
 tstep = 4*sz0*1.e-9;
 ptot = 0.;
 accum = 0.;

 ili = 0;
 i_cool = 0;
 
 for(im=0;im<nspe;im++) {
   factor1 = RP/am[im]*C*sqrt(2*PI)/tstep;
   fscale1[im] = factor1/sqrt(2*PI);
   ptot+=pss[im];
   accum+=crsec[im]*pss[im];
 }
 
 qe0 = nele/(double)NELB;   
 
 nionpb = 0.;
 
 for(im=0;im<nspe;im++) {
   nionpb += crsec[im]*pss[im];
   //accum += pss[im];
   //prob[im] = accum/ptot;
   prob[im] = nionpb/accum;
   fprintf(ion_inphase,"\n Probability ion type %d (mass %lf) = %lf\n ",im,am[im],prob[im]);
 }
 
 nionpb*=3.226e-9*nele*dstep;
 
 qion = nionpb/(double)NION;
 
 omegar = 2*PI*freqr*1e9;
 alphat = omegar/(2*merit);
 omegabar = sqrt(omegar*omegar - alphat*alphat);
 
 omegarz = 2*PI*freqrz*1e6;
 alphaz = omegarz/(2*meritz);
 omegabarz = sqrt(omegarz*omegarz - alphaz*alphaz);
 
 mu *= PI/180.;

 //f_max = fodol/(4*sin(mu/2.));
 //f_min = -fodol/(4*sin(mu/2.));
  
 if(i_lattice==0) {
   beta_max = fodol*(1+sin(mu/2.))/sin(mu);
   beta_min = fodol*(1-sin(mu/2.))/sin(mu);
   alpha_max = (-1-sin(mu/2.))/cos(mu/2.);
   alpha_min = (1-sin(mu/2.))/cos(mu/2.);
   sx_min = sqrt(beta_min*emxN0*1.e-9/gammain);
   sy_min = sqrt(beta_min*emyN0*1.e-9/gammain);
   sx_max = sqrt(beta_max*emxN0*1.e-9/gammain);
   sy_max = sqrt(beta_max*emyN0*1.e-9/gammain);

   sx0 = sx_min;
   sy0 = sy_max;
   betax0 = beta_min;
   betay0 = beta_max;

   coeff1 = sqrt(beta_max/beta_min);
   coeff2 = sqrt(beta_max*beta_min);
   coeff3 = 1/coeff1;
   coeff4 = 1/coeff2;
   
   nstep = nfodo*2;
   gammaprev = gammain;
 
   cxya = cos(mu/2.);
   sxya = sin(mu/2.);
 }
 else {
   coeff1x = sqrt(betax1/betax0);
   coeff1xbis = alfax0;
   coeff2x = sqrt(betax1*betax0);
   coeff3x = 1/coeff1x;
   coeff3xbis = alfax1;
   coeff4x = 1/coeff2x;
   coeff4xbis = 1. + alfax1*alfax0;
   coeff4xter = alfax0 - alfax1;
   coeff1y = sqrt(betay1/betay0);
   coeff1ybis = alfay0;
   coeff2y = sqrt(betay1*betay0);
   coeff3y = 1/coeff1y;
   coeff3ybis = alfay1;
   coeff4y = 1/coeff2y;
   coeff4ybis = 1. + alfay1*alfay0;
   coeff4yter = alfay0 - alfay1;

   gammaprev = benergy0*1.e3/0.511;   
   sx0 = sqrt(betax0*emxN0*1.e-9/gammaprev);
   sx1 = sqrt(betax1*emxN0*1.e-9/gammaprev);
   sy0 = sqrt(betay0*emyN0*1.e-9/gammaprev);
   sy1 = sqrt(betay1*emyN0*1.e-9/gammaprev);
   
   cxa = cos(mux0);
   sxa = sin(mux0);
   cya = cos(muy0);
   sya = sin(muy0);
   
 }

 
 n_diag2 = nstep/10;
 if(nstep<10) n_diag2=1;

 wakefac = - E*E/(ME*gammaprev*C*C);

 factor2 = RE*sqrt(2*PI)/gammaprev/tstep;
 fscale2 = factor2/sqrt(2*PI);
  
 i_dip = 0;
 bsol = 0.;
 bgrad = 0.;
 bdip = 0.;
 
 eel = 1.0;
 epr = 1.0;

 sx0fion = 7.e-6;
 sy0fion = 4.e-6;
 drrefill = 1.e-7;
 
}

/******************************************************************************
 *** Subroutine     : update_values                                         ***
 *** Effect         : update the parameter values that change with energy   ***
 *** Parameters     :                                                       ***
 *** Gbl var   used :                                                       ***
 *** Gbl var effect :                                                       ***
 *** Constants used :                                                       ***
 *** Subrout.  used : none                                                  ***
 ******************************************************************************/

update_values ()

{  double gamma_ave, dstep0;
 
 //using Jean-Bernard's model of field ionization: after the first 4km, full ionization of the beam volume within 10 bunches
 //this is determined by setting ifi=1

 //if(long_pos[it-1]>=11430.)
 if(long_pos[it-1]>=lengthfi)
   ifi=1;

 if(i_lattice==0){
   gammanext = gammain + (gammafin - gammain)*(double)it/(double)nstep;
   
   sx_min = sqrt(beta_min*emxN0*1.e-9/gammanext);
   sy_min = sqrt(beta_min*emyN0*1.e-9/gammanext);
   sx_max = sqrt(beta_max*emxN0*1.e-9/gammanext);
   sy_max = sqrt(beta_max*emyN0*1.e-9/gammanext);
 
   coeff1 = sqrt(gammaprev/gammanext*beta_max/beta_min);
   coeff2 = sqrt(gammaprev/gammanext*beta_max*beta_min);
   coeff3 = sqrt(gammaprev/gammanext*beta_min/beta_max);
   coeff4 = sqrt(gammaprev/gammanext/beta_max/beta_min);

   fprintf(sigma_pr,"%13.8e  %13.8e  %13.8e  %13.8e  %13.8e\n",gammanext,sx_min,sy_max,sx_max,sy_min);

 }
 else{ 
   gammanext = benergy[it]*1.e3/0.511;
   dstep0 = dstep;
   dstep = long_pos[it]-long_pos[it-1];
   tstepp = dstep/C;
   qion *=dstep/dstep0;

   if (i_cool==0){
     coeffga = sqrt(gammaprev/gammanext);
     coeff1x = sqrt(betarrayx[it]/betarrayx[it-1]);
     coeff1xbis = alfarrayx[it-1];
     coeff2x = sqrt(betarrayx[it]*betarrayx[it-1]);
     coeff3x = 1./coeff1x;
     coeff3xbis = alfarrayx[it];
     coeff4x = 1./coeff2x;
     coeff4xbis = 1. + alfarrayx[it]*alfarrayx[it-1];
     coeff4xter = alfarrayx[it-1] - alfarrayx[it];
     coeff1y = sqrt(betarrayy[it]/betarrayy[it-1]);
     coeff1ybis = alfarrayy[it-1];
     coeff2y = sqrt(betarrayy[it]*betarrayy[it-1]);
     coeff3y = 1./coeff1y;
     coeff3ybis = alfarrayy[it];
     coeff4y = 1./coeff2y;
     coeff4ybis = 1. + alfarrayy[it]*alfarrayy[it-1];
     coeff4yter = alfarrayy[it-1] - alfarrayy[it];   
   }
   if (i_cool==1){
     coeff1x = sqrt(betarrayx[it]/betarrayx[it-1]);
     coeff2x = sqrt(betarrayx[it]*betarrayx[it-1]);
     coeff3x = gammaprev/gammanext*sqrt(1./betarrayx[it]*betarrayx[it-1]);
     coeff4x = gammaprev/gammanext*sqrt(1./betarrayx[it]/betarrayx[it-1]);
     coeff1y = sqrt(betarrayy[it]/betarrayy[it-1]);
     coeff2y = sqrt(betarrayy[it]*betarrayy[it-1]);
     coeff3y = gammaprev/gammanext*sqrt(1./betarrayy[it]*betarrayy[it-1]);
     coeff4y = gammaprev/gammanext*sqrt(1./betarrayy[it]/betarrayy[it-1]);
   }

   sx0 = sqrt(betarrayx[it-1]*emxN0*1.e-9/gammaprev);
   sx1 = sqrt(betarrayx[it]*emxN0*1.e-9/gammanext);
   sy0 = sqrt(betarrayy[it-1]*emyN0*1.e-9/gammaprev);
   sy1 = sqrt(betarrayy[it]*emyN0*1.e-9/gammanext);

   cxa = cos(muarrayx[it] - muarrayx[it-1]);
   sxa = sin(muarrayx[it] - muarrayx[it-1]);
   cya = cos(muarrayy[it] - muarrayy[it-1]);
   sya = sin(muarrayy[it] - muarrayy[it-1]);

   fprintf(sigma_pr,"%13.8e  %13.8e  %13.8e  %13.8e  %13.8e  %13.8e\n",long_pos[it-1],gammanext,sx0,sy0,sx1,sy1);

 }

 //cut field ionization...

 /*
 int im;

 if(ifi==1){ 
   nionpb = 0.;
   
   for(im=0;im<nspe;im++) {
     nionpb += (4.*PI*sx0*sy0*1.e22)*pss[im];
     prob2[im] = pss[im]/ptot;
     fprintf(ion_inphase,"\n Field ionization: step %d \n ", it);
     fprintf(ion_inphase,"\n Probability ion type %d (mass %lf) = %lf\n ",im,am[im],prob2[im]);
   }
   
   nionpb*=3.226e-9*dstep;
   //ion charge is changed to that coming from field ionization
   qionfi[0] = nionpb/(double)NION;
   
 }
 */
 /* 

 more detailed module to calculate also the exact area with field ionization
 dropped for the moment because it's too complicated....

    elec_norm_p= 0.;
    for(pp=0;pp<=1000;pp++){
    xtst = (double)pp*3./1000.*sx0;
    ytst = 0.; 
    xoffc = 0.;
    yoffc = 0.;
    ffrank_(&xtst,&ytst,&xoffc,&yoffc,&sx0,&sy0,&pre,&pim,&ili);
    elec_norm = E*nb/4./PI/EPS0/sz0*sqrt(pre*pre+pim*pim)/5.e9;
    if(elec_norm>1. && elec_norm_p<1.){
    xfi1 = xtst;
    ifi = 1;
    }
    if(elec_norm<1. && elec_norm_p>1.)   
    xfi2 = xtst;
    elec_norm_p = elec_norm;
    
 }
 
 elec_norm_p= 0.;
 for(pp=0;pp<=1000;pp++){
 ytst = (double)pp*3./1000.*sy0;
 xtst = 0.; 
 xoffc = 0.;
 yoffc = 0.;
 ffrank_(&xtst,&ytst,&xoffc,&yoffc,&sx0,&sy0,&pre,&pim,&ili);
 elec_norm = E*nb/4./PI/EPS0/sz0*sqrt(pre*pre+pim*pim)/1.e11;
 if(elec_norm>1 && elec_norm>elec_norm_p)
 yfi1 = ytst;
 else
 yfi2 = ytst;
 elec_norm_p = elec_norm;
 }
 
 if(ifi==1){ 
 nionpb2 = 0.;
 
 for(im=0;im<nspe;im++) {
     nionpb2 += (PI*(yfi1-yfi2)*(xfi1-xfi2)*1.e22)*pss[im];
     prob2[im] = pss[im]/ptot;
     fprintf(ion_inphase,"\n Gas ionization \n ");
     fprintf(ion_inphase,"\n Probability ion type %d (mass %lf) = %lf\n ",im,am[im],prob2[im]);
     }
     
     nionpb2*=3.226e-9*dstep;
     //ion charge is enhanced by the contribution coming from field ionization
     qion += nionpb2/(double)NION;

     }
 */
 //end field ionization check
 
 gamma_ave = (gammanext + gammaprev)/2.;
 factor2 = RE*sqrt(2*PI)/gamma_ave/tstep;
 fscale2 = factor2/sqrt(2*PI);

 /* if 'crsec' is function of energy (gammanext) use these lines to update ... 
    
 nionpb = 0.;
 
 for(im=0;im<nspe;im++) {
 crsec[im]*=1.;
    nionpb += crsec[im]*pss[im];
    }
    
    the loop that follows is not necessary if the cross section dependence on energy is linear!
    
    for(im=0;im<nspe;im++) {
    accum += crsec[im]*pss[im];
    prob[im] = accum/nionpb;
    }
    
    nionpb*=3.226e-9*nele*dstep;
    
    qion = nionpb/(double)NION;
    
 */

 wakefac = - E*E/(ME*gamma_ave*C*C);
 
}

/******************************************************************************
 *** Subroutine     : build_distrib                                         ***
 *** Effect         : Calculate arbitrary distribution of bunch.            ***
 *** Parameters     : vector containing the quantity whose distribution     ***
 ***                  we want to plot, number of elements of this vector    ***
 *** Gbl var effect :                                                       ***
 *** Constants used :                                                       ***
 *** Subrout.  used : none                                                  ***
 ******************************************************************************/


void verteil(double data[], unsigned long nele, double distri[], unsigned long ncell,
             double index[])

{ long j, channel1, channel2;

  double  scaled, mean_data, mean_sqr_data, fraction1, fraction2, data_center,
          data_width ;

  mean_data     = 0 ;
  mean_sqr_data = 0 ;

  for (j=0 ; j<nele ; ++j)
         { 
           mean_data     += data[j] ;
           mean_sqr_data += data[j]*data[j] ; }

  mean_data     /= nele ;
  mean_sqr_data /= nele ;

  data_center = mean_data ;
  data_width  = 10.0*sqrt(mean_sqr_data-mean_data*mean_data) ;

  for (j=0 ; j<ncell ; ++j)
      distri[j] = 0.0 ;

  for (j=0 ; j<nele ; ++j)
         { scaled = (data[j]-data_center+data_width/2.0)/data_width*ncell ;
           channel1 = (long) scaled ;
           if (scaled-channel1 > 0.5)
              { fraction1 = 1.5 - scaled + channel1 ;
                channel2  = channel1 + 1 ;
                fraction2 = 1.0 - fraction1 ; }
            else
              { fraction1 = scaled - channel1 + 0.5 ;
                channel2  = channel1 - 1 ;
                fraction2 = 1.0 - fraction1 ; }
           if (channel1 >= 0 && channel1 < ncell)
              distri[channel1] += fraction1 ;
           if (channel2 >= 0 && channel2 < ncell)
              distri[channel2] += fraction2 ; }

  for (j=0 ; j<ncell ; ++j)
    {
      index[j] = data_center - data_width/2.0 + j*data_width/(ncell-1); 
      distri[j] /= data_width/ncell ;
    }

}

/******************************************************************************
 *** Subroutine     : build_distrib                                         ***
 *** Effect         : Calculate arbitrary distribution of bunch.            ***
 *** Parameters     : vector containing the quantity whose distribution     ***
 ***                  we want to plot, number of elements of this vector    ***
 *** Gbl var effect :                                                       ***
 *** Constants used :                                                       ***
 *** Subrout.  used : none                                                  ***
 ******************************************************************************/

/* creates distributions for ensembles of macroparticles with different chanrges  */

void verteilq(double data[], double data2[], double step, unsigned long nele, double distri[], unsigned long ncell,
             double index[])

{ long   j, channel1, channel2, ind_charge ;

  double  scaled, mean_data, mean_sqr_data, fraction1, fraction2, data_center, 
          data_width ;

  mean_data     = 0 ;
  mean_sqr_data = 0 ;

  for (j=0 ; j<nele ; ++j)
         { 
	   ind_charge = j/NION;
           mean_data     += data[j]*data2[ind_charge]/step ;
           mean_sqr_data += data[j]*data[j]*data2[ind_charge]*data2[ind_charge]/step/step ; }

  mean_data     /= nele ;
  mean_sqr_data /= nele ;

  data_center = mean_data ;
  data_width  = 10.*sqrt(mean_sqr_data-mean_data*mean_data);

  for (j=0 ; j<ncell ; ++j)
      distri[j] = 0.0 ;

  for (j=0 ; j<nele ; ++j)
         { ind_charge = j/NION;
	   scaled = (data[j]*data2[ind_charge]/step-data_center+data_width/2.0)/data_width*ncell ;
           channel1 = (long) scaled ;
           if (scaled-channel1 > 0.5)
              { fraction1 = 1.5 - scaled + channel1 ;
                channel2  = channel1 + 1 ;
                fraction2 = 1.0 - fraction1 ; }
            else
              { fraction1 = scaled - channel1 + 0.5 ;
                channel2  = channel1 - 1 ;
                fraction2 = 1.0 - fraction1 ; }
           if (channel1 >= 0 && channel1 < ncell)
              distri[channel1] += fraction1 ;
           if (channel2 >= 0 && channel2 < ncell)
              distri[channel2] += fraction2 ; }

  for (j=0 ; j<ncell ; ++j)
    {
      index[j] = data_center - data_width/2.0 + j*data_width/(ncell-1); 
      distri[j] /= data_width/ncell ;
    }

}

/******************************************************************************
 *** Subroutine     : build_distrib2                                        ***
 *** Effect         : Calculate arbitrary distribution.                     ***
 *** Parameters     : vector containing the quantity whose distribution     ***
 ***                  we want to plot, number of elements of this vector    ***
 *** Gbl var effect :                                                       ***
 *** Constants used :                                                       ***
 *** Subrout.  used : none                                                  ***
 ******************************************************************************/

/* specially suitable to plot rectangular distributions, having a sharp edge  */

void verteil2(double data[], unsigned long nele, double distri[], unsigned long ncell,
             double index[])

{ long   j, channel1, channel2;

  double  scaled, mean_data, mean_sqr_data, fraction1, fraction2, data_center,
          data_width;

  mean_data     = 0 ;
  mean_sqr_data = 0 ;

  for (j=0 ; j<nele ; ++j)
         { 
           mean_data     += data[j] ;
           mean_sqr_data += data[j]*data[j] ; }

  mean_data     /= nele ;
  mean_sqr_data /= nele ;

  data_center = mean_data ;
  data_width  = 5.0*sqrt(mean_sqr_data-mean_data*mean_data) ;

  for (j=0 ; j<ncell ; ++j)
      distri[j] = 0.0 ;

  for (j=0 ; j<nele ; ++j)
         { scaled = (data[j]-data_center+data_width/2.0)/data_width*ncell ;
           channel1 = (long) scaled ;
           if (scaled-channel1 > 0.5)
              { fraction1 = 1.5 - scaled + channel1 ;
                channel2  = channel1 + 1 ;
                fraction2 = 1.0 - fraction1 ; }
            else
              { fraction1 = scaled - channel1 + 0.5 ;
                channel2  = channel1 - 1 ;
                fraction2 = 1.0 - fraction1 ; }
           if (channel1 >= 0 && channel1 < ncell)
              distri[channel1] += fraction1 ;
           if (channel2 >= 0 && channel2 < ncell)
              distri[channel2] += fraction2 ; }

  for (j=0 ; j<ncell ; ++j)
    {
      index[j] = data_center - data_width/2.0 + j*data_width/(ncell-1); 
      distri[j] /= data_width/ncell ;
    }

}

/******************************************************************************
 *** Subroutine     : generate                                              ***
 *** Effect         : generate and initialize electrons turn by turn        ***
 *** Parameters     : none                                                  ***
 *** Gbl var   used : qe,xe,ye,xpe,ype,rpipex,rpipey,it,nelpc               ***
 *** Gbl var used   : press,crse,circ,npr0,pleff,plppm,elflag               ***
 *** Constants used : PI                                                    ***
 *** Subrout.  used : none                                                  ***
 ******************************************************************************/

int generate ()
     
{ long i,j,*accum,im,pp,qq;
 double cx0,cy0,sxx0,syy0,casual,vtherm,accumul,xtst,ytst,xoffc,yoffc,maxe,pre,pim,elec_norm,xionlim,yionlim,
   xionlim2,ymax,randradius,randangle,r1,r2;

 cx0=xoff;
 cy0=yoff;
 sxx0=sx;
 syy0=sy;
 accum=ivector(0,nspe-1);

 //start bunch by bunch field ionization check
 //limit for field ionization to be set in the input file: elec_thresh [GV/m]

 //ifi = 0;
 maxe = 0.;
 ymax = 0.;
 xionlim = 0.;
 yionlim = 0.;
 xionlim2 = 0.;

 //in the following we scan a region [-4sx,4sx] around the bunch axis (x=0 first, and y=0 then) 

 if(jmain==0) {
   xtst = 0.;
   for(qq=1000;qq>=0;qq--){
     ytst = (double)qq*10./1000.*sxx0; 
     xoffc = 0.;
     yoffc = 0.;
     ffrank_(&xtst,&ytst,&xoffc,&yoffc,&sxx0,&syy0,&pre,&pim,&ili);
     //elec_norm = E*nele/4./PI/EPS0/sz0/0.3*sqrt(pre*pre+pim*pim)/elec_thresh/1.e9;
     elec_norm = E*nele/4./PI/EPS0/sz0/0.3*sqrt(pre*pre)/elec_thresh/1.e9;
     if(elec_norm>1.0){
       //ifi = 1;
       //printf("qqy = %d \n", qq);
       qq = -1;
       yionlim=ytst;
     }
   }  

   xtst = 0.;
   for(qq=1000;qq>=0;qq--){
     ytst = (double)qq*10./1000.*sxx0; 
     xoffc = 0.;
     yoffc = 0.;
     ffrank_(&xtst,&ytst,&xoffc,&yoffc,&sxx0,&syy0,&pre,&pim,&ili);
     //elec_norm = E*nele/4./PI/EPS0/sz0/0.3*sqrt(pre*pre+pim*pim)/elec_thresh/1.e9;
     elec_norm = E*nele/4./PI/EPS0/sz0/0.3*sqrt(pre*pre)/elec_thresh/1.e9;
     if(elec_norm>maxe) {
       maxe = elec_norm;
       ymax = ytst;
     }     
   }  
  
   ytst = 0.;
   for(pp=1000;pp>=0;pp--){
     xtst = (double)pp*10./1000.*sxx0; 
     xoffc = 0.;
     yoffc = 0.;
     ffrank_(&xtst,&ytst,&xoffc,&yoffc,&sxx0,&syy0,&pre,&pim,&ili);
     //elec_norm = E*nele/4./PI/EPS0/sz0/0.3*sqrt(pre*pre+pim*pim)/elec_thresh/1.e9;
     elec_norm = E*nele/4./PI/EPS0/sz0/0.3*sqrt(pim*pim)/elec_thresh/1.e9;
     if(elec_norm>maxe) maxe = elec_norm;
     
     if(elec_norm>1.0){
       //ifi = 1;
       //printf("qqx = %d \n", pp);
       pp = -1;
       xionlim=xtst;
     }

   }

   ytst = ymax;
   for(pp=1000;pp>=0;pp--){
     xtst = (double)pp*10./1000.*sxx0; 
     xoffc = 0.;
     yoffc = 0.;
     ffrank_(&xtst,&ytst,&xoffc,&yoffc,&sxx0,&syy0,&pre,&pim,&ili);
     elec_norm = E*nele/4./PI/EPS0/sz0/0.3*sqrt(pre*pre+pim*pim)/elec_thresh/1.e9;
     
     if(elec_norm>1.0){
       //ifi = 1;
       //printf("qqx = %d \n", pp);
       pp = -1;
       xionlim2=xtst;
     }
   
   }  
   
   fprintf(wakefield_pr,"%13.8e  %13.8e  %13.8e  %13.8e  %13.8e  %13.8e\n",long_pos[it-1],sxx0,syy0,xionlim,yionlim,xionlim2);

 //printf("\n step %d, elec_norm_max = %lf \n",it,maxe);

 // change by ADRIAN OEFTIGER 15 Dec '10: use calculated area for field ionization
    sx0fion = xionlim2;
    sy0fion = yionlim;
    if ( xionlim2>1.e-7 ) { ifi = 1; } else { ifi = 0; }
 // end change
  
 }

 for (j=0;j<nspe;j++)
   accum[j]=0;
 
 switch(ifi) {
  
 case 0:
   
   //if(jmain==0) printf("\n Scattering ionization \n");
   if(jmain==0){ 
     fprintf(ion_inphase,"\n Scattering ionization: turn %ld step %ld bunch %ld \n ",nt,it,jmain);     
     fprintf(wakefield2_pr,"%13.8e  %13.8e\n",qion/dstep,qion/dstep);
   }

   qionfi[jmain] = qion;
   for (i=jmain*NION; i<(jmain+1)*NION; i++)
     {  
       xion[i] = cx0 + sxx0*sqrt(-2*log(rand ()/(double)RAND_MAX))*cos(2*PI*rand ()/(double)RAND_MAX);
       yion[i] = cy0 + syy0*sqrt(-2*log(rand ()/(double)RAND_MAX))*cos(2*PI*rand ()/(double)RAND_MAX);
       vprionx[i-jmain*NION] = xion[i];
       vpriony[i-jmain*NION] = yion[i];
       vprionr[i-jmain*NION] = sqrt(xion[i]*xion[i]+yion[i]*yion[i]);
       xpion[i] = 1.e-10*(rand ()/(double)RAND_MAX - 0.5);
       ypion[i] = 1.e-10*(rand ()/(double)RAND_MAX - 0.5);
       casual = rand ()/(double)RAND_MAX;
       for (j=0;j<nspe;j++)
	 if (casual < prob[j]) {
	   mion[i]=j;
	   accum[j]+=1;
	   j=nspe+1;
	 }
       
       if(fabs(xion[i])>xion_max) xion_max=fabs(xion[i]);
       if(fabs(yion[i])>yion_max) yion_max=fabs(yion[i]);
       //printf("%d  %lf  %lf  %d  %lf  %lf \n", 
       //  i, xion[i]*1.e9, yion[i]*1.e9, mion[i], xion_max*1.e9, yion_max*1.e9);
     }
   /*
     for (j=0;j<nspe;j++) {
     if(j<(nspe-1)) printf("\t species %d: %d ", j, accum[j]);
     if(j==(nspe-1)) printf("\t species %d: %d\n ", j, accum[j]);
     }
   */
   free_ivector(accum,0,nspe-1);
   break;
   
 case 1:   
   if(jmain<=10) {
     //printf("\n Field ionization \n");
     nionpb = 0.;
     accumul = 0.;

     fprintf(ion_inphase,"\n Field ionization: turn %ld step %ld bunch %ld \n ",nt,it,jmain);     
     for(im=0;im<nspe;im++) {
       nionpb += (PI*sx0fion*sy0fion*1.e22)*pss[im]/10.;
       accumul +=pss[im];
       prob2[im] = accumul/ptot;
       fprintf(ion_inphase,"\n Probability ion type %ld (mass %lf) = %lf\n ",im,am[im],prob2[im]);
     }
     
     nionpb*=3.226e-9*dstep;
     //ion charge is changed to that coming from field ionization
     qionfi[jmain] = nionpb/(double)NION;
     if(jmain==10) fprintf(wakefield2_pr,"%13.8e  ",qionfi[jmain]/dstep);
     
     for (i=jmain*NION; i<(jmain+1)*NION; i++)
       {  
	 /*
	 do { xion[i] = cx0 + sx0fion*2.*(rand ()/(double)RAND_MAX - 0.5);
	 yion[i] = cy0 + sy0fion*2.*(rand ()/(double)RAND_MAX - 0.5);
	 vprionx[i-jmain*NION] = xion[i];
	 vpriony[i-jmain*NION] = yion[i];
	 vprionr[i-jmain*NION] = sqrt(xion[i]*xion[i]+yion[i]*yion[i]);
	 }
	 while (xion[i]*xion[i]/sx0fion/sx0fion + yion[i]*yion[i]/sy0fion/sy0fion >= 1);
	 */
	 
	 randangle = 2*PI*rand ()/(double)RAND_MAX;
	 randradius = sqrt(rand ()/(double)RAND_MAX);
	 xion[i] = cx0 + sx0fion*randradius*cos(randangle);
	 yion[i] = cy0 + sy0fion*randradius*sin(randangle);
    
	 vprionx[i-jmain*NION] = xion[i];
	 vpriony[i-jmain*NION] = yion[i];
	 vprionr[i-jmain*NION] = sqrt(xion[i]*xion[i]+yion[i]*yion[i]);
	 xpion[i] = 1.e-10*(rand ()/(double)RAND_MAX - 0.5);
	 ypion[i] = 1.e-10*(rand ()/(double)RAND_MAX - 0.5);
	 casual = rand ()/(double)RAND_MAX;
	 for (j=0;j<nspe;j++)
	   if (casual < prob2[j]) {
	     mion[i]=j;
	     accum[j]+=1;
	     j=nspe+1;
	   }
	 
	 if(fabs(xion[i])>xion_max) xion_max=fabs(xion[i]);
	 if(fabs(yion[i])>yion_max) yion_max=fabs(yion[i]);
	 //printf("%d  %lf  %lf  %d  %lf  %lf \n", 
	 //  i, xion[i]*1.e9, yion[i]*1.e9, mion[i], xion_max*1.e9, yion_max*1.e9);
       }
     /*
       for (j=0;j<nspe;j++) {
       if(j<(nspe-1)) printf("\t species %d: %d ", j, accum[j]);
       if(j==(nspe-1)) printf("\t species %d: %d\n ", j, accum[j]);
       }
     */
     free_ivector(accum,0,nspe-1);
   }

   if(jmain>10) { 
     nionpb = 0.;
     accumul = 0.;
     if(jmain==11) fprintf(ion_inphase,"\n Field ionization: turn %ld step %ld bunch %ld \n ", nt, it, jmain);     

     for(im=0;im<nspe;im++) {
       vtherm = sqrt(8.*4.14e-21/PI/am[im]/MU);
       accumul += (PI/4.*(sx0fion+sy0fion)*vtherm)*pss[im];
     }

     for(im=0;im<nspe;im++) {
       vtherm = sqrt(8.*4.14e-21/PI/am[im]/MU);
       nionpb += (PI/4.*(sx0fion+sy0fion)*vtherm)*pss[im];
       prob2[im] = nionpb/accumul;
       if(jmain==11) fprintf(ion_inphase,"\n Probability ion type %ld (mass %lf) = %lf\n ",im,am[im],prob2[im]);
     }
     
     nionpb*=32125.6*dstep*bsp;
     //ion charge is changed to that coming from field ionization
     qionfi[jmain] = nionpb/(double)NION;
     if(jmain==11) fprintf(wakefield2_pr,"%13.8e\n",qionfi[jmain]/dstep);
     
     for (i=jmain*NION; i<(jmain+1)*NION; i++)
       {  
	 randangle = 2*PI*rand ()/(double)RAND_MAX;
	 randradius = rand ()/(double)RAND_MAX;
	 /*
	 xion[i] = cx0 + (sx0fion - rand ()/(double)RAND_MAX*drrefill)*cos(randangle);
	 yion[i] = cy0 + (sy0fion - rand ()/(double)RAND_MAX*drrefill)*sin(randangle);
	 */
	 r1 = sqrt((sx0fion - drrefill)*(sx0fion - drrefill) + randradius*drrefill*(2.*sx0fion - drrefill));
	 r2 = sqrt((sy0fion - drrefill)*(sy0fion - drrefill) + randradius*drrefill*(2.*sy0fion - drrefill));
	 xion[i] = cx0 + r1*cos(randangle);
	 yion[i] = cy0 + r2*sin(randangle);
	 vprionx[i-jmain*NION] = xion[i];
	 vpriony[i-jmain*NION] = yion[i];
	 vprionr[i-jmain*NION] = sqrt(xion[i]*xion[i]+yion[i]*yion[i]);
	 xpion[i] = 1.e-10*(rand ()/(double)RAND_MAX - 0.5);
	 ypion[i] = 1.e-10*(rand ()/(double)RAND_MAX - 0.5);
	 casual = rand ()/(double)RAND_MAX;
	 for (j=0;j<nspe;j++)
	   if (casual < prob2[j]) {
	     mion[i]=j;
	     accum[j]+=1;
	     j=nspe+1;
	   }
	 
	 if(fabs(xion[i])>xion_max) xion_max=fabs(xion[i]);
	 if(fabs(yion[i])>yion_max) yion_max=fabs(yion[i]);
	 //printf("%d  %lf  %lf  %d  %lf  %lf \n", 
	 //  i, xion[i]*1.e9, yion[i]*1.e9, mion[i], xion_max*1.e9, yion_max*1.e9);
       }
     /*
       for (j=0;j<nspe;j++) {
       if(j<(nspe-1)) printf("\t species %d: %d ", j, accum[j]);
       if(j==(nspe-1)) printf("\t species %d: %d\n ", j, accum[j]);
       }
     */
     free_ivector(accum,0,nspe-1);
     
   }
   break;
 }
 
 return 0;
}


/******************************************************************************
 *** Subroutine     : bunch properties                                      ***
 *** Effect         : Calculates bunch per bunch centroids and sizes        ***
 *** Parameters     : none                                                  ***
 *** Gbl var   used :                                                       ***
 *** Gbl var effect :                                                       ***
 *** Constants used :                                                       ***
 *** Subrout.  used : none                                                  ***
 ******************************************************************************/

bunch_properties ()

{ double xsum, ysum, xpsum, ypsum, x2sum, y2sum, xp2sum, yp2sum, xpprsum, ypprsum, xel_max_loc, yel_max_loc;
 long i,j;  
  
 for(i=0; i<NBUN; i++) {

   xsum = 0.0;
   ysum = 0.0; 
   xpsum = 0.0;
   ypsum = 0.0; 
   x2sum = 0.0;
   y2sum = 0.0;
   xp2sum = 0.0;
   yp2sum = 0.0;
   xpprsum = 0.0;
   ypprsum = 0.0;
   xel_max_loc = 0.0;
   yel_max_loc = 0.0;


   for (j=i*NELB; j<(i+1)*NELB; j++) {
     xsum += xel[j];
     ysum += yel[j];
     xpsum += xpel[j];
     ypsum += ypel[j];
     if(fabs(xel[j])>xel_max_loc) xel_max_loc = fabs(xel[j]);
     if(fabs(yel[j])>yel_max_loc) yel_max_loc = fabs(yel[j]);
   }
 
   xsum /= (double)NELB;
   ysum /= (double)NELB;
   xpsum /= (double)NELB;
   ypsum /= (double)NELB;

   for (j=i*NELB; j<(i+1)*NELB; j++)
     { x2sum += (xel[j]-xsum)*(xel[j]-xsum);
     y2sum += (yel[j]-ysum)*(yel[j]-ysum);
     xp2sum += (xpel[j]-xpsum)*(xpel[j]-xpsum);
     yp2sum += (ypel[j]-ypsum)*(ypel[j]-ypsum);
     xpprsum += (xel[j]-xsum)*(xpel[j]-xpsum);
     ypprsum += (yel[j]-ysum)*(ypel[j]-ypsum);
   }

   xs[i] = xsum;
   ys[i] = ysum;
   sxs[i] = sqrt(x2sum/(double)NELB);
   sys[i] = sqrt(y2sum/(double)NELB);
   sxp2s[i] = xp2sum/(double)NELB;
   syp2s[i] = yp2sum/(double)NELB;
   xppr[i] = xpprsum/(double)NELB;
   yppr[i] = ypprsum/(double)NELB;
   xel_max[i] = xel_max_loc;
   yel_max[i] = yel_max_loc;

 }

}

/******************************************************************************
 *** Subroutine     : ion properties                                        ***
 *** Effect         : Calculates ion distribution features                  ***
 *** Parameters     : none                                                  ***
 *** Gbl var   used :                                                       ***
 *** Gbl var effect :                                                       ***
 *** Constants used :                                                       ***
 *** Subrout.  used : none                                                  ***
 ******************************************************************************/

ion_properties ()

{ double xsum, ysum, x2sum, y2sum;
 long j, ioncount;  

 /* 07.06.2010 
    we decide to calculate the ion properties only over those ions 
    which are within a radius of 10sigma from the current bunch
*/
  
 xsum = 0.0;
 ysum = 0.0; 
 x2sum = 0.0;
 y2sum = 0.0;
 ioncount = 0;

 for (j=0; j<ntemp; j++) {
     xsum += xion[j];
     ysum += yion[j];
     ioncount++;
 }
 
 xsum /= (double)ioncount;
 ysum /= (double)ioncount;

 for (j=0; j<ntemp; j++){
     x2sum += (xion[j]-xsum)*(xion[j]-xsum);
     y2sum += (yion[j]-ysum)*(yion[j]-ysum);
 }

 xsion = xsum;
 ysion = ysum;
 sxsion = sqrt(x2sum/(double)ioncount);
 sysion = sqrt(y2sum/(double)ioncount);
 
}

/******************************************************************************
 *** Subroutine     : initialization                                        ***
 *** Effect         : Initialize electron position and velocity             ***
 *** Parameters     : none                                                  ***
 *** Gbl var   used :                                                       ***
 *** Gbl var effect :                                                       ***
 *** Constants used :                                                       ***
 *** Subrout.  used : none                                                  ***
 ******************************************************************************/

initialization ()

{ double xkick, ykick, xsum, ysum, xpsum, ypsum, u, v, s;
  long i,j;  

 /*
 if(i_lattice==0) { 
   sx0 = sx_min;
   sy0 = sy_max;
   betax0 = beta_min;
   betay0 = beta_max;
 }
 */

 for (j=0; j<NBUN; j++) {
   xsum = 0.0;
   ysum = 0.0;
   xpsum = 0.0;
   ypsum = 0.0;
   xkick = 0.0;
   ykick = 0.0;

   switch (i_kick){

   case 0:
     break;

   case 1:
     xkick = nxkick0*sx0;
     ykick = nykick0*sy0;
     break;

   case 2:
     xkick = nxkick0*sx0*sin(2*PI*(double)j*(double)h_pertx/(double)(NBUN-1));
     ykick = nykick0*sy0*sin(2*PI*(double)j*(double)h_perty/(double)(NBUN-1));
     break;

   case 3:
     xkick = nxkick0*sx0*(2.*rand () / (double)RAND_MAX - 1.);
     ykick = nykick0*sy0*(2.*rand () / (double)RAND_MAX- 1.);
     break;

   default:
     break;

   }

   switch(i_lattice){

   case 0:
     for (i=j*NELB; i<(j+1)*NELB; i++){ 
       do { u = 2.0 * rand () / (double)RAND_MAX - 1.0 ;
       v = 2.0 * rand () / (double)RAND_MAX - 1.0 ;
       s = u*u + v*v ; }
       while (s >= 1.0) ;
     
       xel[i] = sx0 * u * sqrt(-2.0*log(s)/s);
       xpel[i] = sx0/betax0 * v * sqrt(-2.0*log(s)/s);
     
     
       do { u = 2.0 * rand () / (double)RAND_MAX - 1.0 ;
       v = 2.0 * rand () / (double)RAND_MAX - 1.0 ;
       s = u*u + v*v ; }
       while (s >= 1.0) ;
     
       yel[i] = sy0 * u * sqrt(-2.0*log(s)/s);
       ypel[i] = sy0/betay0 * v * sqrt(-2.0*log(s)/s);
     
       xsum += xel[i];
       ysum += yel[i];
       xpsum += xpel[i];
       ypsum += ypel[i];
     }
     break;
     
   case 1:
     for (i=j*NELB; i<(j+1)*NELB; i++){
       
       u = sqrt(-2.0 * log( rand()/(double)(RAND_MAX) ));
       v = 2.0*PI * ( rand()/(double)(RAND_MAX) );

       xel[i] = sx0 * u * cos(v);
       xpel[i] = -sx0/betax0 * u * (sin(v) + alfax0 * cos(v));

       u = sqrt(-2.0 * log( rand()/(double)(RAND_MAX) ));
       v = 2.0*PI * ( rand()/(double)(RAND_MAX) );

       yel[i] = sy0 * u * cos(v);
       ypel[i] = -sy0/betay0 * u * (sin(v) + alfay0 * cos(v));

       xsum += xel[i];
       ysum += yel[i];
       xpsum += xpel[i];
       ypsum += ypel[i];
     }
     break;
   } 
   
   xave = xsum/(double)NELB;
   yave = ysum/(double)NELB;
   xpave = xpsum/(double)NELB;
   ypave = ypsum/(double)NELB;
   /*   
   for (i=j*NELB; i<(j+1)*NELB; i++)
     { xel[i] += xkick - xave;
     xpel[i] += -xpave;
     yel[i] += ykick - yave;
     ypel[i] += -ypave;
     }
   */
   for (i=j*NELB; i<(j+1)*NELB; i++)
     { xel[i] += xkick;
     yel[i] += ykick;
     }
 }

}

/******************************************************************************
 *** Subroutine     : wake_function                                         ***
 *** Effect         : Function that gives the wake function at a given      ***
 ***                  location z < 0 (z is 'pos')                           ***
 *** Parameters     : pos                                                   ***
 *** Gbl var   used :                                                       ***
 *** Gbl var effect :                                                       ***
 *** Constants used : none                                                  ***
 *** Subrout.  used : none                                                  ***
 ******************************************************************************/

double wake_func (double pos)

{ double field;

  field = (zetat*1e6*omegar*omegar/merit/omegabar)*exp(alphat*pos/C)*sin(omegabar*pos/C);

  return field;

}

/******************************************************************************
 *** Subroutine     : wake_function2                                        ***
 *** Effect         : Function that gives the wake function at a given      ***
 ***                  location z < 0 (z is 'pos') for a resistive wall      ***
 *** Parameters     : pos                                                   ***
 *** Gbl var   used :                                                       ***
 *** Gbl var effect :                                                       ***
 *** Constants used : none                                                  ***
 *** Subrout.  used : none                                                  ***
 ******************************************************************************/

double wake_reswall (double pos, double piper, double length)

{ double sigma, lambdas, mur, wfield ;

 sigma = 1.4e6; //stainless steel in SI [1/Ohm/m] 
 mur=1;
 
 lambdas = 1/(Z0*sigma);
 wfield = -C*Z0*length/PI/pow(piper,3.)*sqrt(-lambdas*mur/PI/pos);
 
 /* field = -1/PI/piper/piper/piper/2/PI/EPS0*sqrt(C/sigma)*1/sqrt(-pos)*length;
    this formula can be correclty used to get the wake field in SI but giving the
    conductivity in CGS units (1/s) */
 
 return wfield;

}
   
/******************************************************************************
 *** Subroutine     : wake_function2                                        ***
 *** Effect         : Function that gives the wake function at a given      ***
 ***                  location z < 0 (z is 'pos') for a resistive wall      ***
 ***                  with inductive by-pass (A. Koschik)                   ***
 *** Parameters     : pos                                                   ***
 *** Gbl var   used :                                                       ***
 *** Gbl var effect :                                                       ***
 *** Constants used : none                                                  ***
 *** Subrout.  used : none                                                  ***
 ******************************************************************************/

double wake_reswall_ib (double pos, double piper, double length)

{  double sigma, lambdas, mur, wfield ;

 sigma = 1.4e6; //stainless steel in SI [1/Ohm/m] 
 mur=1.;
 
 lambdas = 1/(Z0*sigma);
 wfield = -C*Z0*length/PI/pow(piper,3.)*sqrt(-lambdas*mur/PI/pos);
 // wfield += 
 //2.*C*length*mur/pow(piper,4.)/PI/sigma*exp(-4.*lambdas*mur*pos/piper/piper)*
 //(1. - gsl_sf_erf (sqrt(-4.*lambdas*mur*pos/piper/piper)));
 
 return wfield;
 
}

/******************************************************************************
 *** Subroutine     : longitudinal wake_function                            ***
 *** Effect         : Function that gives the z-wake function at a given    ***
 ***                  location z < 0 (z is 'pos') - single bunch in cavity  ***
 *** Parameters     : pos                                                   ***
 *** Gbl var   used :                                                       ***
 *** Gbl var effect :                                                       ***
 *** Constants used : none                                                  ***
 *** Subrout.  used : none                                                  ***
 ******************************************************************************/

double wake_funcz (double pos)

{ double field ;

  field = (rshz*1e6*2*alphaz)*exp(alphaz*pos/C)*(cos(omegabarz*pos/C) + alphaz/omegabarz*sin(omegabarz*pos/C));
  return field;

}

/******************************************************************************
 *** Subroutine     : open_files                                            ***
 *** Effect         : Open files for storing calculation data.              ***
 *** Parameters     : none                                                  ***
 *** Gbl var   used :                                                       ***
 *** Gbl var effect :                                                       ***
 *** Constants used : none                                                  ***
 *** Subrout.  used : none                                                  ***
 ******************************************************************************/

open_files ()

{ char   chap[3]              ,
    ions_filename[80]         ,            /* ion trajectories over whole run    */
    cntr_filename[80]         ,            /* train properties                   */ 
    trainhdtl_fname[80]       ,            /* train head-tail distribution       */ 
    elephase_fname[80]        ,            /* ....                               */
    ion_filename[80]          ,            /* some general information           */
    iini_filename[80]         ,            /* ions                               */
    eleini_filename[80]       ,
    samp_fname[80]            ,
    samp2_fname[80]           ,
    sigma_filename[80]        ,
    track_filename[80]        ,       
    grid_filename[80]         ;       
 long i                       ;
 
 strcpy (cntr_filename,filename) ;
 strcat (cntr_filename,"_prt.dat") ;
 centr_pr = fopen(cntr_filename,"w");
 
 strcpy (trainhdtl_fname,filename) ;
 strcat (trainhdtl_fname,"_hdtl.dat") ;
 trainhdtl_pr = fopen(trainhdtl_fname,"w");
 
 strcpy (elephase_fname,filename) ;
 strcat (elephase_fname,"_eleb.dat") ;
 ele_phase = fopen(elephase_fname,"w");
 
 strcpy (ion_filename,filename) ;
 strcat (ion_filename,"_ionph.dat") ;
 ion_inphase = fopen(ion_filename,"w");
 
 for(i=0;i<15;i++){
   itoa(i,chap);
   strcpy (iini_filename,filename) ;
   strcat (iini_filename,chap) ;
   strcat (iini_filename,"_iinix.dat") ;
   ion_inix[i] = fopen(iini_filename,"w");
   
   strcpy (eleini_filename,filename) ;
   strcat (eleini_filename,chap) ;
   strcat (eleini_filename,"_iiniy.dat") ;
   ion_iniy[i] = fopen(eleini_filename,"w");
 
   strcpy (samp_fname,filename) ;
   strcat (samp_fname,chap) ;
   strcat (samp_fname,"_sample.dat") ;
   samp_ion[i] = fopen(samp_fname,"w");
   
   strcpy (samp2_fname,filename) ;
   strcat (samp2_fname,chap) ;
   strcat (samp2_fname,"_sample2.dat") ;
   samp2_ion[i] = fopen(samp2_fname,"w");
 }

 strcpy (track_filename,filename) ;
 strcat (track_filename,"_trk.dat") ;
 wakefield_pr = fopen(track_filename,"w");

 strcpy (track_filename,filename) ;
 strcat (track_filename,"_trk2.dat") ;
 wakefield2_pr = fopen(track_filename,"w");

 strcpy (sigma_filename,filename) ;
 strcat (sigma_filename,"_sgm.dat") ;
 sigma_pr = fopen(sigma_filename,"w");

 strcpy(ions_filename, filename);
 strcat(ions_filename, "_ions.dat");
 ions_pr = fopen(ions_filename, "w");

 strcpy(grid_filename, filename);
 strcat(grid_filename, "_grid.dat");
 grid_pr = fopen(grid_filename, "w");
 
}

/******************************************************************************
 ****                             PIC module                               ****
 ******************************************************************************/

void*
xmalloc (size_t size) {
  void* tmp;
  tmp=malloc(size);
  if (tmp!=NULL) {
    return tmp;
  }
  fprintf(stderr,"ERROR: not enough memory\n");
  exit(1);
}

/*
  Distributes the beam particles for field calculation
*/

void particles_distribute(GRID *grid,double x[],double y[],int n, int i_which,
                          double n_macro)
{
  int i,j,i1,i2;
  double d_x,d_y,delta_x_i,delta_y_i,offset_x,offset_y;
  double cut_x,cut_y,max_x,max_y;
  double h_x,h_y;
  double *rho;
  int n_x,n_y;
  int distribute_in=0,distribute_out=0;

  n_x=grid->n_cell_x;
  n_y=grid->n_cell_y;
  cut_x=grid->cut_x;
  max_x=grid->max_x;
  cut_y=grid->cut_y;
  max_y=grid->max_y;

  if (i_which==1) {
    rho=grid->rho1;
  }
  else {
    rho=grid->rho2;
  }
  d_x=grid->delta_x;
  d_y=grid->delta_y;
  delta_x_i=1.0/grid->delta_x;
  delta_y_i=1.0/grid->delta_y;
  offset_x=grid->offset_x;
  offset_y=grid->offset_y;
  for (i=0;i<n_x*n_y;i++){
    rho[i]=0.0;
  }
  switch (grid->distribution) {
  case 1:
    for (i=0;i<n;i++){
      i1=(int)floor(x[i]*delta_x_i+offset_x);
      i2=(int)floor(y[i]*delta_y_i+offset_y);
      if ((i1>=0)&&(i1<n_x)&&(i2>=0)&&(i2<n_y)){
        rho[i1*n_y+i2] += n_macro;
        distribute_in++;
      }
      else{
        distribute_out++;
      }
    }
    break;
  case 2:
    for (i=0;i<n;i++){
      if ((x[i]>=-cut_x+0.5*d_x)&&(x[i]<max_x+0.5*d_x)&&(y[i]>=-cut_y+0.5*d_y)&&(y[i]<max_y+0.5*d_y)){
	h_x=(x[i]*delta_x_i+offset_x-0.5);
	i1=(int)h_x;
	h_x -= (float)i1;
	h_y=(y[i]*delta_y_i+offset_y-0.5);
	i2=(int)h_y;
	h_y -= (float)i2;
	  j=i1*n_y+i2;
	  rho[j] += (1.0-h_x)*(1.0-h_y)*n_macro;
	  j=(i1+1)*n_y+i2;
	  rho[j] += h_x*(1.0-h_y)*n_macro;
	  j=i1*n_y+i2+1;
	  rho[j] += (1.0-h_x)*h_y*n_macro;
	  j=(i1+1)*n_y+i2+1;
	  rho[j] += h_x*h_y*n_macro;
	  distribute_in++;
	}
	else{
	  distribute_out++;
	}
    }
    break;
  }
  /*
    double densita;
  if (it==1){
  if (i_which==2){
     for (i1=0;i1<n_x;i1++)
     {  densita=0.;
      for (i2=0;i2<n_y;i2++)
	densita+=rho[i1*n_y+i2];
  	fprintf(ele_ini,"%d %13.8e\n",i1,densita); 
	//	printf("%d %13.8e\n",i1,densita); 
       
     }
      fprintf(ele_ini,"\n");
     fprintf(ele_ini,"\n");
  }}
  */
  /*
  if (it==1){
  if (i_which==2){
    for (i2=0;i2<n_y;i2++)
     {  densita=0.;
      for (i1=0;i1<n_x;i1++)
	densita+=rho[i1*n_y+i2];
  	fprintf(prot_ini,"%d %13.8e\n",i2,densita); 
	
     }
    fprintf(prot_ini,"\n");
    fprintf(prot_ini,"\n");


  }}
  */
 
  /*
     if (i_which==2){
     for (i1=0;i1<n_x;i1++)
      for (i2=0;i2<n_y;i2++)
  	fprintf(ele_ini,"%d %d %13.8e\n",i1,i2,rho[i1*n_y+i2]); //, sqrt((i1-n_x/2)*(i1-n_x/2)+(i2-n_y/2)*(i2-n_y/2))); 
	//	printf("%d %d %13.8e\n",i1,i2,rho[i1*n_y+i2]); }
	  
     // fprintf(ele_ini,"\n");
     // fprintf(ele_ini,"\n");
     }
  */

}

void particles_distribute2(GRID *grid,double x[],double y[],int n,int i_which, 
                           double q[])
{
  int i,j,i1,i2,ind_charge;
  double d_x,d_y,delta_x_i,delta_y_i,offset_x,offset_y;
  double cut_x,cut_y,max_x,max_y;
  double h_x,h_y;
  double *rho;
  int n_x,n_y;
  int distribute_in=0,distribute_out=0;

  n_x=grid->n_cell_x;
  n_y=grid->n_cell_y;
  cut_x=grid->cut_x;
  max_x=grid->max_x;
  cut_y=grid->cut_y;
  max_y=grid->max_y;

  if (i_which==1) {
    rho=grid->rho1;
  }
  else {
    rho=grid->rho2;
  }
  d_x=grid->delta_x;
  d_y=grid->delta_y;
  delta_x_i=1.0/grid->delta_x;
  delta_y_i=1.0/grid->delta_y;
  offset_x=grid->offset_x;
  offset_y=grid->offset_y;
  for (i=0;i<n_x*n_y;i++){
    rho[i]=0.0;
  }
  switch (grid->distribution) {
  case 1:
    for (i=0;i<n;i++){
      ind_charge=i/NION;
      i1=(int)floor(x[i]*delta_x_i+offset_x);
      i2=(int)floor(y[i]*delta_y_i+offset_y);
      if ((i1>=0)&&(i1<n_x)&&(i2>=0)&&(i2<n_y)){
        rho[i1*n_y+i2] += q[ind_charge];
        distribute_in++;
      }
      else{
        distribute_out++;
      }
    }
    break;
  case 2:
    for (i=0;i<n;i++){
      if ((x[i]>=-cut_x+0.5*d_x)&&(x[i]<max_x+0.5*d_x)&&(y[i]>=-cut_y+0.5*d_y)&&(y[i]<max_y+0.5*d_y)){
	ind_charge=i/NION;
	h_x=(x[i]*delta_x_i+offset_x-0.5);
	i1=(int)h_x;
	h_x -= (float)i1;
	h_y=(y[i]*delta_y_i+offset_y-0.5);
	i2=(int)h_y;
	h_y -= (float)i2;
        j=i1*n_y+i2;
        rho[j] += (1.0-h_x)*(1.0-h_y)*q[ind_charge];
        j=(i1+1)*n_y+i2;
        rho[j] += h_x*(1.0-h_y)*q[ind_charge];
        j=i1*n_y+i2+1;
        rho[j] += (1.0-h_x)*h_y*q[ind_charge];
        j=(i1+1)*n_y+i2+1;
        rho[j] += h_x*h_y*q[ind_charge];
        distribute_in++;
      }
      else{
        distribute_out++;
      }
    }
    break;
  }
}

/**********************************************/
/* Routines for the calculation of the fields */
/**********************************************/

/* This routine calculates the potentials on the outer boundary of the grid
   by direct calculation. */

void foldfronts (double *rho,double *dist,double *phi,int n_x,int n_y)
{
  double s;
  int i1,i2,i3,i4;

  i2=0;
  for (i1=0;i1<n_x;i1++)
    {
      s=0.0;
      for (i3=0;i3<n_x;i3++)
        {
          for (i4=0;i4<n_y;i4++)
            {
              s+=rho[i3*n_y+i4]
                *dist[abs(i1-i3)*n_y+abs(i2-i4)];
            }
        }
      phi[i1*n_y+i2]=s;
    }
  i2=n_y-1;
  for (i1=0;i1<n_x;i1++)
    {
      s=0.0;
      for (i3=0;i3<n_x;i3++)
        {
          for (i4=0;i4<n_y;i4++)
            {
              s+=rho[i3*n_y+i4]
                *dist[abs(i1-i3)*n_y+abs(i2-i4)];
            }
        }
      phi[i1*n_y+i2]=s;
    }
  i1=0;
  for (i2=1;i2<n_y-1;i2++)
    {
      s=0.0;
      for (i3=0;i3<n_x;i3++)
        {
          for (i4=0;i4<n_y;i4++)
            {
              s+=rho[i3*n_y+i4]
                *dist[abs(i1-i3)*n_y+abs(i2-i4)];
            }
        }
      phi[i1*n_y+i2]=s;
    }
  i1=n_x-1;
  for (i2=1;i2<n_y-1;i2++)
    {
      s=0.0;
      for (i3=0;i3<n_x;i3++)
        {
          for (i4=0;i4<n_y;i4++)
            {
              s+=rho[i3*n_y+i4]
                *dist[abs(i1-i3)*n_y+abs(i2-i4)];
            }
        }
      phi[i1*n_y+i2]=s;
    }
}


 #include "fourtrans22.c"
// #include "../../../fourtrans2.c"


/* Fast fourier transformation routine used by foldfft */

void fourtrans (double* data,int nn[],int isign)
{
  /* System generated locals */
  int i__1, i__2, i__3, i__4, i__5, i__6;
  double d__1;
  
  /* Builtin functions */
  //  double sin();
  
  /* Local variables */
  int idim, ibit, nrem, ntot, i2rev, i3rev, n;
  double theta, tempi, tempr;
  int nprev, i2;
  double wtemp;
  int i1, i3, k1, k2;
  double wi, wr;
  int ip1, ip2, ip3;
  double wpi, wpr;
  int ifp1, ifp2;
  int ndim=2;
  /* Parameter adjustments */
  --nn;
  --data;
  
  /* Function Body */
  ntot = 1;
  i__1 = ndim;
  for (idim = 1; idim <= i__1; ++idim)
    {
      ntot *= nn[idim];
    }
  nprev = 1;
  i__1 = ndim;
  for (idim = 1; idim <= i__1; ++idim) {
    n = nn[idim];
    nrem = ntot / (n * nprev);
    ip1 = nprev << 1;
    ip2 = ip1 * n;
    ip3 = ip2 * nrem;
    i2rev = 1;
    i__2 = ip2;
    i__3 = ip1;
    for (i2 = 1; i__3 < 0 ? i2 >= i__2 : i2 <= i__2; i2 += i__3) {
      if (i2 < i2rev) {
        i__4 = i2 + ip1 - 2;
        for (i1 = i2; i1 <= i__4; i1 += 2) {
          i__5 = ip3;
          i__6 = ip2;
          for (i3 = i1; i__6 < 0 ? i3 >= i__5 : i3 <= i__5; i3 += 
               i__6) {
            i3rev = i2rev + i3 - i2;
            tempr = data[i3];
            tempi = data[i3 + 1];
                        data[i3] = data[i3rev];
            data[i3 + 1] = data[i3rev + 1];
            data[i3rev] = tempr;
            data[i3rev + 1] = tempi;
          }
        }
      }
      ibit = ip2 / 2;
    L1:
      if (ibit >= ip1 && i2rev > ibit) {
        i2rev -= ibit;
        ibit /= 2;
        goto L1;
      }
      i2rev += ibit;
    }
    ifp1 = ip1;
  L2:
    if (ifp1 < ip2) {
      ifp2 = ifp1 << 1;
      theta = isign * 6.28318530717959 / (ifp2 / ip1);
      /* Computing 2nd power */
      d__1 = sin(theta * .5);
      wpr = d__1 * d__1 * -2.;
      wpi = sin(theta);
      wr = 1.;
      wi = 0.;
      i__3 = ifp1;
      i__2 = ip1;
      for (i3 = 1; i__2 < 0 ? i3 >= i__3 : i3 <= i__3; i3 += i__2) {
        i__4 = i3 + ip1 - 2;
        for (i1 = i3; i1 <= i__4; i1 += 2) {
          i__6 = ip3;
                    i__5 = ifp2;
          for (i2 = i1; i__5 < 0 ? i2 >= i__6 : i2 <= i__6; i2 += 
               i__5) {
            k1 = i2;
            k2 = k1 + ifp1;
            tempr = wr * data[k2] - wi * data[k2 + 1];
            tempi = wr * data[k2 + 1] + wi * data[k2];
            data[k2] = data[k1] - tempr;
            data[k2 + 1] = data[k1 + 1] - tempi;
            data[k1] += tempr;
            data[k1 + 1] += tempi;
          }
        }
        wtemp = wr;
        wr = wr * wpr - wi * wpi + wr;
        wi = wi * wpr + wtemp * wpi + wi;
      }
      ifp1 = ifp2;
      goto L2;
    }
    nprev = n * nprev;
  }
  return;
}

/* This routine calculates the potentials with the help of the fast fourier
   transform */

void fold_fft(double *rho1,double *rho2,double *dist,double *phi1,
              double *phi2,double temp[],int n_x,int n_y)
{
  double help;
  int i1,i2,i,j,j0;
  int nn[2];

  /* return no field */

  /* fill array with charge */

  for (i=0;i<n_x*n_y*8;i++) {
    temp[i]=0.0;
  }
  for (i1=0;i1<n_x;i1++) {
    for (i2=0;i2<n_y;i2++) {
      j=2*(i1*2*n_y+i2);
      j0=i1*n_y+i2;
      temp[j]=rho1[j0];
      temp[j+1]=rho2[j0];
    }
  }
  
  nn[0]=2*n_y;
  nn[1]=2*n_x;
  i1=2;
  i2=1;
  //fourtrans(temp,nn,i2);
  fourtrans3(temp,nn,i2);
  for (i=0;i<n_x*n_y*4;i++) {
    j=2*i;
    help=temp[j]*dist[j]-temp[j+1]*dist[j+1];
    temp[j+1]=temp[j+1]*dist[j]+temp[j]*dist[j+1];
    temp[j]=help;
  }
  i2=-1;
  //fourtrans(temp,nn,i2);
  fourtrans3(temp,nn,i2);
  for (i1=0;i1<n_x;i1++) {
    for (i2=0;i2<n_y;i2++) {
      phi1[i1*n_y+i2]=temp[2*(i1*2*n_y+i2)]/((double)4*n_x*n_y);
      phi2[i1*n_y+i2]=temp[2*(i1*2*n_y+i2)+1]/((double)4*n_x*n_y);
    }
  }
}

void init_sor (double *a,double *b,double *c,double *d,
double *e,GRID grid)
{
  double a0,b0,c0,d0,e0,factor;
  int i;

  factor=1.0;
  a0=factor*(grid.delta_y/grid.delta_x);
  b0=a0;
  c0=factor*(grid.delta_x/grid.delta_y);
  d0=c0;
  e0=-2.0*(a0+c0);
  for (i=0;i<grid.n_cell_x*grid.n_cell_y;i++)
    {
      a[i]=a0;
      b[i]=b0;
      c[i]=c0;
      d[i]=d0;
      e[i]=e0;
    }
}

void sor (double *a,double *b,double *c,double *d,double *e,
          double *rho, double *phi,int n_x,int n_y)
{
  const int MAXIT=1000;
  const double eps=1e-5;
  double anormf=0.0,anorm;
  double omega,resid;
  int i1,i2,n,j;
  float rjac,delta_x=1000.0,delta_y=64.0;

  rjac=(delta_y*delta_y*cos(PI/(float)n_x)+delta_x*delta_x*cos(PI/(float)n_y))
    /(delta_x*delta_x+delta_y*delta_y);
  for (i1=1;i1<n_x-1;i1++)
    {
      for (i2=1;i2<n_y-1;i2++)
        {
          anormf += fabs(rho[i1*n_y+i2]);
        }
    }
  omega=1.0;
  for (n=0;n<MAXIT;n++)
    {
      anorm=0.0;
      for (i1=1;i1<n_x-1;i1++)
        {
          for (i2=1;i2<n_y-1;i2++)
            {
              if (((i1+i2) % 2)==(n % 2))
                {
                  j=i1*n_y+i2;
                  resid=a[j]*phi[j+n_y]+b[j]*phi[j-n_y]+c[j]*phi[j+1]
                    +d[j]*phi[j-1]+e[j]*phi[j]-rho[j];
                  anorm += fabs(resid);
                  phi[j] -= omega*resid/e[j];
                }
            }
        }
      if (n==0)
        {
          omega=1.0/(1.0-0.5*rjac*rjac);
        }
      else
        {
          omega=1.0/(1.0-0.25*rjac*rjac*omega);
        }
      if ((n>0)&&(anorm<eps*anormf))
        {
          return;
        }
    }
  printf ("Mist");
  scanf ("%d",&i1);
}

/* Routine to calculate the parameters necessary for the iterative method of
   the potential calculation sor2 */

void init_sor2 (GRID grid,double *parameter)
{
  double factor;
  factor=1.0;
  parameter[0]=factor*(grid.delta_y/grid.delta_x);
  parameter[1]=parameter[0];
  parameter[2]=factor*(grid.delta_x/grid.delta_y);
  parameter[3]=parameter[2];
  parameter[4]=-2.0*(parameter[0]+parameter[2]);
  parameter[5]=(grid.delta_y*grid.delta_y*cos(PI/(float)grid.n_cell_x)
                +grid.delta_x*grid.delta_x*cos(PI/(float)grid.n_cell_y))
                        /(grid.delta_x*grid.delta_x+grid.delta_y*grid.delta_y);
}

/* Routine to calculate the potentials with an iterative method */

void sor2 (double *rho,double *phi,int n_x,int n_y,double *parameter)
{
  const int MAXIT=1000;
  const double eps=1e-5;
  double anormf=0.0,anorm;
  double omega,resid,e_inv;
  int i1,i2,n,j;
  double a,b,c,d,e,rjac;

  a=parameter[0];
  b=parameter[1];
  c=parameter[2];
  d=parameter[3];
  e=parameter[4];
  rjac=parameter[5];
  for (i1=1;i1<n_x-1;i1++)
    {
      for (i2=1;i2<n_y-1;i2++)
        {
          anormf += fabs(rho[i1*n_y+i2]);
        }
    }
  omega=1.0;
  e_inv=1.0/e;
  for (n=0;n<MAXIT;n++)
    {
      anorm=0.0;
      for (i1=1;i1<n_x-1;i1++)
        {
          for (i2=1;i2<n_y-1;i2++)
            {
              if (((i1+i2) % 2)==(n % 2))
                {
                  j=i1*n_y+i2;
                  resid=a*phi[j+n_y]+b*phi[j-n_y]+c*phi[j+1]
                    +d*phi[j-1]+e*phi[j]-rho[j];
                  anorm += fabs(resid);
                  phi[j] -= omega*resid*e_inv;
                }
            }
        }
      if (n==0)
        {
          omega=1.0/(1.0-0.5*rjac*rjac);
        }
      else
        {
          omega=1.0/(1.0-0.25*rjac*rjac*omega);
        }
      if ((n>0)&&(anorm<eps*anormf))
        {
          return;
        }
    }
  printf ("Warning: too many iterations in function sor2\n");
}

double
f_potential(double x,double y)
{
  double sum;
  sum=x*y*(log(x*x+y*y)-3.0)+x*x*atan(y/x)+y*y*atan(x/y);
  return sum;
}

double
f_potential_2(double x0,double y0,double dx,double dy)
{
  double x,y,sum;
  x=x0+dx;
  y=y0+dy;
  sum=x*y*(log(x*x+y*y)-3.0)
      +x*x*atan(y/x)+y*y*atan(x/y);
  x=x0+dx;
  y=y0-dy;
  sum-=x*y*(log(x*x+y*y)-3.0)
      +x*x*atan(y/x)+y*y*atan(x/y);
  x=x0-dx;
  y=y0+dy;
  sum-=x*y*(log(x*x+y*y)-3.0)
      +x*x*atan(y/x)+y*y*atan(x/y);
  x=x0-dx;
  y=y0-dy;
  sum+=x*y*(log(x*x+y*y)-3.0)
      +x*x*atan(y/x)+y*y*atan(x/y);
  return sum;
}

/* This routine initialises the grid for the calculation of the potentials.
   what_grid decides what distribution will be used, what_grid=2 should be
   prefered. */

void
grid_init_phys (GRID *grid,float cut_x,float cut_y)
{
   int i1,i2,j,j0,j_tl,j_bl,j_tr,j_br;
   double factor,phi0;
   double x0,y0,x,y,d_x,d_y;
   int n_x,n_y;
   int nn[2];

   n_x=grid->n_cell_x;
   n_y=grid->n_cell_y;
   grid->min_x=-((float)n_x-2)/((float)n_x)*cut_x;
   grid->max_x=((float)n_x-2)/((float)n_x)*cut_x;
   grid->min_y=-((float)n_y-2)/((float)n_y)*cut_y;
   grid->max_y=((float)n_y-2)/((float)n_y)*cut_y;
   grid->cut_x = cut_x;
   grid->cut_y = cut_y;
  
   grid->offset_x=0.5*(float)n_x;
   grid->offset_y=0.5*(float)n_y;

   grid->delta_x=2.0*cut_x/((double)n_x);
   grid->delta_y=2.0*cut_y/((double)n_y);
   d_x=grid->delta_x;
   d_y=grid->delta_y;
   factor=1.0/(d_x*d_y);
   switch (grid->integration_method) {
   case 1:
   case 3:
     phi0=factor*
       f_potential_2(0.0,0.0,0.5*grid->delta_x,0.5*grid->delta_y);
     for (i1=0;i1<grid->n_cell_x;i1++) {
       for (i2=0;i2<grid->n_cell_y;i2++) {
         j=i1*grid->n_cell_y+i2;
         x0=(double)i1;
         y0=(double)i2;
         grid->dist[j]=factor*
           f_potential_2(x0*grid->delta_x,y0*grid->delta_y,
                         0.5*grid->delta_x,0.5*grid->delta_y)
           -phi0;
       }
     }
     break;
   case 2:
     for (i1=0;i1<n_x*n_y*8;i1++) {
       grid->dist[i1]=0.0;
     }
     for (i1=0;i1<=n_x;i1++) {
       x=i1*d_x-0.5*d_x;
       for (i2=0;i2<=n_y;i2++) {
	 y=(double)i2*d_y-0.5*d_y;	 
	 j=i1*(n_y+1)+i2;
	 grid->pot[j]=f_potential(x,y);
       }
     }
     phi0=factor*(grid->pot[n_y+2]-grid->pot[n_y+1]
				-grid->pot[1]+grid->pot[0]);
     for (i1=0;i1<n_x;i1++) {
       for (i2=0;i2<n_y;i2++) {
         j0=2*(i1*2*n_y+i2);
         x0=(double)i1;
         y0=(double)i2;
	 j_tr=(i1+1)*(n_y+1)+(i2+1);
	 j_br=(i1+1)*(n_y+1)+i2;
	 j_tl=i1*(n_y+1)+(i2+1);
	 j_bl=i1*(n_y+1)+i2;
         grid->dist[j0]=factor*(grid->pot[j_tr]-grid->pot[j_br]
				-grid->pot[j_tl]+grid->pot[j_bl])-phi0;
         if (i2>=1) {
           grid->dist[2*((i1+1)*2*n_y-i2)]=grid->dist[j0];
           if (i1>=1) {
             grid->dist[2*(2*n_y*(2*n_x-i1)+i2)]=grid->dist[j0];
             grid->dist[2*(2*n_y*(2*n_x-i1)+2*n_y-i2)]=grid->dist[j0];
           }
         }
         else {
           if (i1>=1) {
             grid->dist[2*(2*n_y*(2*n_x-i1)+i2)]=grid->dist[j0];
           }
         }
       }
     }
     nn[0]=2*n_y;
     nn[1]=2*n_x;
     i1=2;
     i2=1;
     //fourtrans(grid->dist,nn,i2);
     fourtrans2(grid->dist,nn,i2);
     break;
   }
}

GRID*
grid_init_comp (int n_x,int n_y,int integration_method)
{
  GRID *grid;
  
  grid=(GRID*)xmalloc(sizeof(GRID));
  grid->n_cell_x=n_x;
  grid->n_cell_y=n_y;

  grid->integration_method=integration_method;
  grid->distribution=2;
  grid->interpolation=2;

  grid->rho1=(double*)xmalloc(sizeof(double)*n_x*n_y);
  grid->rho2=(double*)xmalloc(sizeof(double)*n_x*n_y);
  grid->phi1=(double*)xmalloc(sizeof(double)*n_x*n_y);
  grid->phi2=(double*)xmalloc(sizeof(double)*n_x*n_y);
  grid->pot=(double*)xmalloc(sizeof(double)*(n_x+1)*(n_y+1));

  if (integration_method==2) {
    grid->dist=(double*)xmalloc(sizeof(double)*n_x*n_y*8);
    grid->temp=(double*)xmalloc(sizeof(double)*n_x*n_y*8);
  }
  else {
    grid->dist=(double*)xmalloc(sizeof(double)*n_x*n_y);
  }
  return grid;
}


/* This routine folds the charge distribution and the greensfunction of
a grid to get the potentials. */

void foldfields (double *rho,double *dist,double *phi,int n_x,int n_y)
{
  double s;
  int i1,i2,i3,i4;

  for (i1=0;i1<n_x;i1++) {
    for (i2=0;i2<n_y;i2++) {
      s=0.0;
      for (i3=0;i3<n_x;i3++) {
        for (i4=0;i4<n_y;i4++) {
          s+=rho[i3*n_y+i4]*dist[abs(i1-i3)*n_y+abs(i2-i4)];
        }
      }
      phi[i1*n_y+i2]=s;
    }
  }
}


void field_calculate(GRID *grid)
{
  switch (grid->integration_method) {
  case 1:
    foldfields(grid->rho1,grid->dist,grid->phi1,grid->n_cell_x,grid->n_cell_y);
    foldfields(grid->rho2,grid->dist,grid->phi2,grid->n_cell_x,grid->n_cell_y);
    break;
  case 2:
    fold_fft(grid->rho1,grid->rho2,grid->dist,grid->phi1,grid->phi2,
             grid->temp,grid->n_cell_x,grid->n_cell_y);
    break;
  }
}

void particles_move(GRID *grid,double energy,double ics0[],
                    double yps0[],double icsp[],double ypsp[],
                    int nn,int i_which,double step,double scale)
{ int i,i1,i2;
  int elconta, rediffcont;

  double ax,ay;
  double phi1_x,phi2_x,phi3_x,phi1_y,phi2_y,phi3_y,h_x,h_y;
  int n_x,n_y;
  register int j;
  register double h,h_p;
  double *phi;
  double delta_x_i,delta_y_i,offset_x,offset_y,delta_x,delta_y;
  double xx,yy,xxp,yyp;
  double min_x,max_x,min_y,max_y;
  double k11,k12,k13,k14,k21,k22,k23,k24,k31,k32,k33,k34,k41,k42,k43,k44;
  double velprx, velpry, bux, buy, bmod;

  elconta=0;
  rediffcont=0;

  min_x=grid->min_x;
  max_x=grid->max_x;
  min_y=grid->min_y;
  max_y=grid->max_y;
  delta_x=grid->delta_x;
  delta_y=grid->delta_y;
  delta_x_i=1.0/grid->delta_x;
  delta_y_i=1.0/grid->delta_y;

  offset_x=grid->offset_x;
  offset_y=grid->offset_y;
  n_x=grid->n_cell_x;
  n_y=grid->n_cell_y;

    if (i_which==1){    //to move the electrons
      phi=grid->phi2;
    }
    else{               //to move the ions
      phi=grid->phi1;
    }

  switch (grid->interpolation){
  case 1:
    for (i=0;i<nn;i++){
      i1=(int)floor(ics0[i]*delta_x_i+offset_x);
      i2=(int)floor(yps0[i]*delta_y_i+offset_y);
      if ((i1<n_x-1)&&(i1>0)&&(i2<n_y-1)&&(i2>0)){
        j=i1*n_y+i2;
        ax = (phi[j+n_y]
              -phi[j-n_y])*step
          /(2.0*delta_x*energy);
        ay = (phi[j+1]
              -phi[j-1])*step
          /(2.0*delta_y*energy);
        icsp[i] -= ax*scale;
        ypsp[i] -= ay*scale;
      }
      else {
        fprintf(stderr,"WARNING: particle is not on the grid\n");
      }
      ics0[i] += icsp[i]*step;
      yps0[i] += ypsp[i]*step;
    }
    break;

    case 2:
      for (i=0;i<nn;i++){
        if(i_which==2 && i_dip!=2) {
          ics0[i] += icsp[i]*step;
          yps0[i] += ypsp[i]*step;
	}
	
        xx=ics0[i];
        yy=yps0[i];
        xxp=icsp[i];
        yyp=ypsp[i];
        if ((xx>min_x)&&(xx<max_x)&&(yy>min_y)&&(yy<max_y)){
          h_x=xx*delta_x_i+offset_x;
          h_y=yy*delta_y_i+offset_y;

          i1=(int)h_x;
          h=h_y-0.5;
          i2=(int)h;
          h -= i2;
          j=i1*n_y+i2;
          h_p=1.0-h;
          phi1_y=h*phi[j+n_y+1]+h_p*phi[j+n_y];
          phi2_y=h*phi[j+1]+h_p*phi[j];
          phi3_y=h*phi[j-n_y+1]+h_p*phi[j-n_y];

          i2=(int)h_y;
          h=h_x-0.5;
          i1=(int)h;
          h -= i1;
          h_p=1.0-h;
          j=i1*n_y+i2;
          phi1_x=h*phi[j+n_y+1]+h_p*phi[j+1];
          phi2_x=h*phi[j+n_y]+h_p*phi[j];
          phi3_x=h*phi[j+n_y-1]+h_p*phi[j-1];

          h_x -= (int)h_x;
          h_y -= (int)h_y;
          ax = (h_x*(phi1_y-phi2_y)+(1.0-h_x)*(phi2_y-phi3_y))/(delta_x*energy);
          ay = (h_y*(phi1_x-phi2_x)+(1.0-h_y)*(phi2_x-phi3_x))/(delta_y*energy);

          switch(i_dip){

            case 0:
                 icsp[i] -= ax*step*scale;
                 ypsp[i] -= ay*step*scale;
                 break;	 

	    case 1:
	      if(i_which==2)
                 icsp[i] -= 0.0;
	      else
		icsp[i] -= ax*step*scale;

                 ypsp[i] -= ay*step*scale;
                 break;

	    case 2:
                 if(i_which==2)
                  { k11 = xxp*step;
                    k21 = (-ax*scale + E/ELMS*bsol*yyp)*step;
                    k31 = yyp*step;
                    k41 = (-ay*scale - E/ELMS*bsol*xxp)*step;

                    k12 = (xxp+k21/2)*step;
                    k22 = (-ax*scale + E/ELMS*bsol*(yyp+k41/2))*step;
                    k32 = (yyp+k41/2)*step;
                    k42 = (-ay*scale - E/ELMS*bsol*(xxp+k21/2))*step;

                    k13 = (xxp+k22/2)*step;
                    k23 = (-ax*scale + E/ELMS*bsol*(yyp+k42/2))*step;
                    k33 = (yyp+k42/2)*step;
                    k43 = (-ay*scale - E/ELMS*bsol*(xxp+k22/2))*step;

                    k14 = (xxp+k23)*step;
                    k24 = (-ax*scale + E/ELMS*bsol*(yyp+k43))*step;
                    k34 = (yyp+k43)*step;
                    k44 = (-ay*scale - E/ELMS*bsol*(xxp+k23))*step;

                    ics0[i] += k11/6 + k12/3 + k13/3 + k14/6;
                    icsp[i] += k21/6 + k22/3 + k23/3 + k24/6;
                    yps0[i] += k31/6 + k32/3 + k33/3 + k34/6;
                    ypsp[i] += k41/6 + k42/3 + k43/3 + k44/6;

		  }

                 if(i_which==1)
                  { icsp[i] -= ax*scale*step;
   	            ypsp[i] -= ay*scale*step;
		  }
                 break;

	    case 3:

	         if(i_which==2)
                  { velprx = icsp[i] - ax*step*scale;
                    velpry = ypsp[i] - ay*step*scale;
                    bux = bgrad*yps0[i];
                    buy = bdip - bgrad*ics0[i];
                    bmod = sqrt(bux*bux+buy*buy);
                    bux *= fabs(velprx)/velprx/bmod;
                    buy *= fabs(velpry)/velpry/bmod;
                    icsp[i] = (bux*velprx + buy*velpry)*bux;
                    ypsp[i] = (bux*velprx + buy*velpry)*buy;
		  }
                 if(i_which==1)
                  { icsp[i] -= ax*scale*step;
   	            ypsp[i] -= ay*scale*step;
		  }
                 break;
	  }
	}

        else{
	  // COMMENTED OUT BY ADRIAN OEFTIGER:
          //fprintf(stderr,"WARNING: particle is not on the grid\n");
	  //fprintf(stderr,"WARNING: it's an electron \n"); 
        }
      }
      break;
  }
}

void particles_move_ion(GRID *grid,double energy,double ics0[],
                    double yps0[],double icsp[],double ypsp[],long mm[],
                    int nn,int i_which,double step)
{ static float scale;
  int i,i1,i2;
  
  int elconta, rediffcont;

  double ax,ay;
  double phi1_x,phi2_x,phi3_x,phi1_y,phi2_y,phi3_y,h_x,h_y;
  int n_x,n_y;
  register int j;
  register double h,h_p;
  double *phi;
  double delta_x_i,delta_y_i,offset_x,offset_y,delta_x,delta_y;
  double xx,yy,xxp,yyp;
  double min_x,max_x,min_y,max_y;
  double k11,k12,k13,k14,k21,k22,k23,k24,k31,k32,k33,k34,k41,k42,k43,k44;
  double velprx, velpry, bux, buy, bmod;

  elconta=0;
  rediffcont=0;

  min_x=grid->min_x;
  max_x=grid->max_x;
  min_y=grid->min_y;
  max_y=grid->max_y;
  delta_x=grid->delta_x;
  delta_y=grid->delta_y;
  delta_x_i=1.0/grid->delta_x;
  delta_y_i=1.0/grid->delta_y;

  offset_x=grid->offset_x;
  offset_y=grid->offset_y;
  n_x=grid->n_cell_x;
  n_y=grid->n_cell_y;

    if (i_which==1){    //to move the electrons
      phi=grid->phi2;
    }
    else{               //to move the ions
      phi=grid->phi1;
    }

  switch (grid->interpolation){
  case 1:
    for (i=0;i<nn;i++){
      i1=(int)floor(ics0[i]*delta_x_i+offset_x);
      i2=(int)floor(yps0[i]*delta_y_i+offset_y);
      if ((i1<n_x-1)&&(i1>0)&&(i2<n_y-1)&&(i2>0)){
        j=i1*n_y+i2;
        ax = (phi[j+n_y]
              -phi[j-n_y])*step
          /(2.0*delta_x*energy);
        ay = (phi[j+1]
              -phi[j-1])*step
          /(2.0*delta_y*energy);
        icsp[i] -= ax*scale;
        ypsp[i] -= ay*scale;
      }
      else {
        fprintf(stderr,"WARNING: particle is not on the grid\n");
      }
      ics0[i] += icsp[i]*step;
      yps0[i] += ypsp[i]*step;
    }
    break;

    case 2:
      for (i=0;i<nn;i++){
	/*
        if(i_which==2 && i_dip!=2) {
          ics0[i] += icsp[i]*step;
          yps0[i] += ypsp[i]*step;
	}
	*/

	scale=fscale1[mm[i]];
        xx=ics0[i];
        yy=yps0[i];
        xxp=icsp[i];
        yyp=ypsp[i];
        if ((xx>min_x)&&(xx<max_x)&&(yy>min_y)&&(yy<max_y)){
          h_x=xx*delta_x_i+offset_x;
          h_y=yy*delta_y_i+offset_y;
          i1=(int)h_x;
          h=h_y-0.5;
          i2=(int)h;
          h -= i2;
          j=i1*n_y+i2;
          h_p=1.0-h;
          phi1_y=h*phi[j+n_y+1]+h_p*phi[j+n_y];
          phi2_y=h*phi[j+1]+h_p*phi[j];
          phi3_y=h*phi[j-n_y+1]+h_p*phi[j-n_y];
          i2=(int)h_y;
          h=h_x-0.5;
          i1=(int)h;
          h -= i1;
          h_p=1.0-h;
          j=i1*n_y+i2;
          phi1_x=h*phi[j+n_y+1]+h_p*phi[j+1];
          phi2_x=h*phi[j+n_y]+h_p*phi[j];
          phi3_x=h*phi[j+n_y-1]+h_p*phi[j-1];
          h_x -= (int)h_x;
          h_y -= (int)h_y;
          ax = (h_x*(phi1_y-phi2_y)+(1.0-h_x)*(phi2_y-phi3_y))/(delta_x*energy);
          ay = (h_y*(phi1_x-phi2_x)+(1.0-h_y)*(phi2_x-phi3_x))/(delta_y*energy);

          switch(i_dip){

            case 0:
                 icsp[i] -= ax*step*scale;
                 ypsp[i] -= ay*step*scale;
                 break;	 

	    case 1:
	      if(i_which==2)
                 icsp[i] -= 0.0;
	      else
		icsp[i] -= ax*step*scale;

                 ypsp[i] -= ay*step*scale;
                 break;

	    case 2:
                 if(i_which==2)
                  { k11 = xxp*step;
                    k21 = (-ax*scale + E/ELMS*bsol*yyp)*step;
                    k31 = yyp*step;
                    k41 = (-ay*scale - E/ELMS*bsol*xxp)*step;

                    k12 = (xxp+k21/2)*step;
                    k22 = (-ax*scale + E/ELMS*bsol*(yyp+k41/2))*step;
                    k32 = (yyp+k41/2)*step;
                    k42 = (-ay*scale - E/ELMS*bsol*(xxp+k21/2))*step;

                    k13 = (xxp+k22/2)*step;
                    k23 = (-ax*scale + E/ELMS*bsol*(yyp+k42/2))*step;
                    k33 = (yyp+k42/2)*step;
                    k43 = (-ay*scale - E/ELMS*bsol*(xxp+k22/2))*step;

                    k14 = (xxp+k23)*step;
                    k24 = (-ax*scale + E/ELMS*bsol*(yyp+k43))*step;
                    k34 = (yyp+k43)*step;
                    k44 = (-ay*scale - E/ELMS*bsol*(xxp+k23))*step;

                    ics0[i] += k11/6 + k12/3 + k13/3 + k14/6;
                    icsp[i] += k21/6 + k22/3 + k23/3 + k24/6;
                    yps0[i] += k31/6 + k32/3 + k33/3 + k34/6;
                    ypsp[i] += k41/6 + k42/3 + k43/3 + k44/6;

		  }

                 if(i_which==1)
                  { icsp[i] -= ax*scale*step;
   	            ypsp[i] -= ay*scale*step;
		  }
                 break;

	    case 3:

	         if(i_which==2)
                  { velprx = icsp[i] - ax*step*scale;
                    velpry = ypsp[i] - ay*step*scale;
                    bux = bgrad*yps0[i];
                    buy = bdip - bgrad*ics0[i];
                    bmod = sqrt(bux*bux+buy*buy);
                    bux *= fabs(velprx)/velprx/bmod;
                    buy *= fabs(velpry)/velpry/bmod;
                    icsp[i] = (bux*velprx + buy*velpry)*bux;
                    ypsp[i] = (bux*velprx + buy*velpry)*buy;
		  }
                 if(i_which==1)
                  { icsp[i] -= ax*scale*step;
   	            ypsp[i] -= ay*scale*step;
		  }
                 break;
	  }
	}

        else{
          //fprintf(stderr,"WARNING: particle is not on the grid\n");
	  //fprintf(stderr,"WARNING: it's an ion \n"); 
	  ionoutgrid += 1;
        }
      }
      break;
  }
}

void field_poissonbc(GRID *grid)

{
  double    *rho11,
            *rho22,
            *work, 
            *rcos, 
            ddx, 
            ddy;

  int       nn, mm, kloc, nloc;
  
  ddx     = 1.0;
  ddy     = grid->delta_y/grid->delta_x;
  mm     = grid->n_cell_x  ; 
  nn     = grid->n_cell_y ; 

  rho11=dvector(0,mm*nn-1);
  rho22=dvector(0,mm*nn-1);
  work=dvector(0,mm+nn+2);
  rcos=dvector(0,mm*nn-1);

  for (kloc=0; kloc<mm; kloc++)
     for (nloc=0; nloc<nn; nloc++)
       {
	 /* BUG, fixed Gio'05
	 if (kloc!=0 && nloc!=0) 
	   rcos[kloc*mm + nloc]=1.0/(2.0*((cos((double)(kloc)*PI/mm)-1.0)/(ddx*ddx)+(cos((double)(nloc)*PI/nn)-1.0)/(ddy*ddy)))*4.0/(mm*nn);
	 else
	   rcos[0]=0.;

        rho11[kloc*mm + nloc] = (4*PI)*(*(grid->rho1 + nloc*nn + kloc));
        rho22[kloc*mm + nloc] = (4*PI)*(*(grid->rho2 + nloc*nn + kloc));
	*/
	
        rcos[kloc*nn + nloc]=1.0/(2.0*((cos((double)(kloc)*PI/mm)-1.0)/(ddx*ddx)+(cos((double)(nloc)*PI/nn)-1.0)/(ddy*ddy)))*4.0/(mm*nn);

        rho11[kloc*nn + nloc] = (4*PI)*(*(grid->rho1 + kloc*nn + nloc));
        rho22[kloc*nn + nloc] = (4*PI)*(*(grid->rho2 + kloc*nn + nloc));
	
      }

  potfft_(rho22,rcos,work,&mm,&nn); 
  potfft_(rho11,rcos,work,&mm,&nn);



  for (kloc=0; kloc<mm; kloc++){
    for (nloc=0; nloc<nn; nloc++)    
      { 
	/* BUG, fixed Gio'05 
	*(grid->phi1 + nloc*nn + kloc) = 1/ddy*rho11[kloc*mm + nloc];
        *(grid->phi2 + nloc*nn + kloc) = 1/ddy*rho22[kloc*mm + nloc];
		*/
	
        *(grid->phi1 + kloc*nn + nloc) = 1/ddy*rho11[kloc*nn + nloc];
        *(grid->phi2 + kloc*nn + nloc) = 1/ddy*rho22[kloc*nn + nloc];
      }}

  free_dvector(rho11,0,mm*nn-1);
  free_dvector(rho22,0,mm*nn-1);
  //  free_dvector(work,0,mm*nn+2);
  free_dvector(work,0,mm+nn+2);
  free_dvector(rcos,0,mm*nn-1);

}

/******************************************************************************
 ****                            Main Program                              ****
 ******************************************************************************/

int main (int argc, char *argv[])

{ 
  // change by ADRIAN OEFTIGER: if arguments are given, seed random generator with first argument
  if ( argc > 1 ) 
    srand(atoi(argv[1]));

  GRID *grid;
  long dg_index,abs_index;
  double timep, wake_store[NBUN], gridextx, gridexty, forkickx, forkicky;

  printf ("             Program for studying the fast beam-ion instability          \n");
  printf ("                                                                       \n\n");

  printf ("--> Reading data from configuration file.\n") ;
  read_data () ;

  printf ("--> Open files.\n") ;
  open_files ();
  printf ("\t\t Complete\n") ;

  printf ("--> Init values.\n") ;
  init_values () ;
  printf ("\t\t Complete\n") ;

  printf ("--> Generation of bunch distributions.\n");
  initialization ();
  printf ("\t\t Complete\n") ;
  
  printf ("--> Generation of the grid for the particle-in-cell.\n");
  printf ("    calculation of the fields \n");
    
  if(iwake_flag==1){
    switch(i_pipe) {
    case 1:
    case 11:
      for(imain=1;imain<NBUN;imain++){
	wake_store[imain]=wake_func(-(double)imain*bsp);
	//wakez_store[imain]=wake_funcz(-(double)imain*bsp);
	fprintf(wakefield_pr,"%13.8e  %13.8e\n",-(double)imain*bsp*0.3,wake_store[imain]);
      }
      printf("\n\n The wake field is broad-band impedance \n\n");
      break;
      
    case 2:
    case 22:
      for(imain=1;imain<NBUN;imain++){
	wake_store[imain]=wake_reswall(-(double)imain*bsp,rpipey,dstep);
	//wakez_store[imain]=wake_funcz(-(double)imain*bsp);
	fprintf(wakefield_pr,"%13.8e  %13.8e\n",-(double)imain*bsp*0.3,wake_store[imain]);
     }
      //change the longitudinal field!!!!
     printf("\n\n The wake field is resistive wall with classical theory \n\n");
     break;
     
   case 3:
   case 33:
     for(imain=1;imain<NBUN;imain++){
       wake_store[imain]=wake_reswall_ib(-(double)imain*bsp,rpipey,dstep);
       //wakez_store[imain]=wake_funcz(-(double)imain*bsp);
       fprintf(wakefield_pr,"%13.8e  %13.8e\n",-(double)imain*bsp*0.3,wake_store[imain]);
     }
     printf("\n\n The wake field is resistive wall with inductive by-pass \n\n");
     break;

   default:
     printf("\n\n The wake field is not defined, either you switch iwake_flag to 0 or you choose another i_pipe \n\n");
     exit(2);
     break;
   }
 }

  printf ("\t\t Complete\n") ;

  if(iion==1)  
    grid=grid_init_comp(128,128,2);
  
  fprintf (ion_inphase, "electrons per bunch = %d\t  ions per bunch = %d\t number of bunches = %d\n",NELB,NION,NBUN);
  fprintf (ion_inphase, "beta function values    -> beta_max = %g  beta_min = %g\n",beta_max,beta_min);  //ele-Mar05  
  //fprintf (ion_inphase, "initial grid sizes -> gridextx = %g  gridexty = %g\n",gridextx,gridexty);


  /* Loop over turns starts */

  for (nt=1; nt<=nturn; nt++) {

    /* here the loop over the total number of steps (interaction points) 
       starts...                                                             */
  
    for(it=1; it<nstep; it++) { 
      xion_max = 0.0;
      yion_max = 0.0;      
	
      timep = (double)(it - 1)*tstepp;
      printf("%ld %ld %f \n",nt,it,timep);
      
      for(lmain=0;lmain<3;lmain++) {
	xave = 0.0;
	yave = 0.0;
	zave = 0.0;
	xpave = 0.0;
	ypave = 0.0;
	xemit = 0.0;
	yemit = 0.0;
	//invary = 0.0;
	//invarx = 0.0;
	
	for(imain=lmain*NELB*NBUN/3; imain<(lmain+1)*NELB*NBUN/3; imain++) { 
	  xave += xel[imain];
	  yave += yel[imain];
	  xpave += xpel[imain];
	  ypave += ypel[imain];
	}
	  
	xave /= (double)(NELB*NBUN/3);
	yave /= (double)(NELB*NBUN/3);
	xpave /= (double)(NELB*NBUN/3);
	ypave /= (double)(NELB*NBUN/3);
	
	for(imain=lmain*NELB*NBUN/3; imain<(lmain+1)*NELB*NBUN/3; imain++) {
	  sxsize += (xel[imain] - xave)*(xel[imain] - xave);
	  sysize += (yel[imain] - yave)*(yel[imain] - yave);
	  sxpsize += (xpel[imain] - xpave)*(xpel[imain] - xpave);
	  sypsize += (ypel[imain] - ypave)*(ypel[imain] - ypave); 
	  xemit += (xel[imain] - xave)*(xpel[imain] - xpave);
	  yemit += (yel[imain] - yave)*(ypel[imain] - ypave);
	  //invary += (ypr[imain] - yave)*(ypr[imain] - yave) + betay0*betay0*(ypp[imain] - ypave)*(ypp[imain]-ypave);
	  //invarx += (xpr[imain] - xave)*(xpr[imain] - xave) + betax0*betax0*(xpp[imain] - xpave)*(xpp[imain]-xpave);
	}
	  
	sxsize = sqrt(sxsize/(double)(NELB*NBUN/3));
	sysize = sqrt(sysize/(double)(NELB*NBUN/3));
	sxpsize /= (double)(NELB*NBUN/3);
	sypsize /= (double)(NELB*NBUN/3);
	xemit /= (double)(NELB*NBUN/3);
	yemit /= (double)(NELB*NBUN/3);
	//invarx /= (double)(NELB*NBUNCH/3);
	//invary /= (double)(NELB*NBUNCH/3);
	  
	if(lmain==0)
	  fprintf(centr_pr,"%13.8e  %13.8e  %13.8e  %13.8e  %13.8e  %13.8e  %13.8e  %13.8e  ",
		  timep,long_pos[it-1],xave,yave,sxsize,sysize,sqrt(sxsize*sxsize*sxpsize-xemit*xemit),
		  sqrt(sysize*sysize*sypsize-yemit*yemit));
	if(lmain==1)
	  fprintf(centr_pr,"%13.8e  %13.8e  %13.8e  %13.8e  %13.8e  %13.8e  ",
		  xave,yave,sxsize,sysize,sqrt(sxsize*sxsize*sxpsize-xemit*xemit),
		  sqrt(sysize*sysize*sypsize-yemit*yemit));
	if(lmain==2)
	  fprintf(centr_pr,"%13.8e  %13.8e  %13.8e  %13.8e  %13.8e  %13.8e\n",
		  xave,yave,sxsize,sysize,sqrt(sxsize*sxsize*sxpsize-xemit*xemit),
		  sqrt(sysize*sysize*sypsize-yemit*yemit));
      
      }
      bunch_properties();	//defines vectors sxs,sys,xs,ys, xel_max, yel_max
    
      update_values ();
	
      for(jmain=0; jmain<NBUN; jmain++) { 
	sx = sxs[jmain];
	sy = sys[jmain];
	xoff = xs[jmain];
	yoff = ys[jmain];
	sxp2 = sxp2s[jmain];
	syp2 = syp2s[jmain];
	sxppr = xppr[jmain];
	syppr = yppr[jmain];
	localemitx = sqrt(sx*sx*sxp2 - sxppr*sxppr);
	localemity = sqrt(sy*sy*syp2 - syppr*syppr);
	    
	/* offsets and rms sizes of the bunches are printed to file          */
	if(i_lattice==0) {
	  if(it % n_diag == 0 || it == 1) { 
	    if(it%2==1)
	      fprintf(trainhdtl_pr, "%13.8e  %13.8e  %13.8e  %13.8e  %13.8e  %13.8e  %13.8e  %13.8e  %13.8e\n",
		      (double)jmain*bsp, xoff, yoff, sx, sy, 
		      sx*sx/beta_min*gammaprev*1.e9, sy*sy/beta_max*gammaprev*1.e9, 
		      localemitx*gammaprev*1.e9, localemity*gammaprev*1.e9);
	    else
	      fprintf(trainhdtl_pr, "%13.8e  %13.8e  %13.8e  %13.8e  %13.8e  %13.8e  %13.8e  %13.8e  %13.8e\n",
		      (double)jmain*bsp, xoff, yoff, sx, sy, 
		      sx*sx/beta_max*gammaprev*1.e9, sy*sy/beta_min*gammaprev*1.e9, 
		      localemitx*gammaprev*1.e9, localemity*gammaprev*1.e9);
	  }
	}
	else {
	  if(it % n_diag == 0 || it == 1) 
	    fprintf(trainhdtl_pr, "%13.8e  %13.8e  %13.8e  %13.8e  %13.8e  %13.8e  %13.8e  %13.8e  %13.8e\n",
		    (double)jmain*bsp, xoff, yoff, sx, sy, 
		    sx*sx/betarrayx[it-1]*gammaprev*1.e9, sy*sy/betarrayy[it-1]*gammaprev*1.e9, 
		    localemitx*gammaprev*1.e9, localemity*gammaprev*1.e9);
	}
	    
	//printf ("\n bunch number = %d, sx = %lf, sy = %lf, xoff = %lf, yoff = %lf\n",jmain,sx,sy,xoff,yoff);
	    
	if (iion==1) {
	      
	  ionoutgrid = 0;
	  generate();
	      
	  ntemp = (jmain+1)*NION;
	  
	  ion_properties();
	  //printf ("\n it = %d, n_diag2 = %d, ntemp = %d \n",it,n_diag2,ntemp);	
	  
	  abs_index = (it-1)/n_diag2;
	  dg_index = (it-1)%n_diag2;
	  //printf ("\n abs_index = %d, dg_index = %d, ntemp = %d \n",abs_index,dg_index,ntemp);	
	  if (dg_index==0) {
	    verteilq (xion,qionfi,dstep,ntemp,hdp,NDC,indxp);
	    verteil (vprionx,NION,hdp2,NDC,indxp2);	
	    for (lmain=0; lmain<NDC; lmain++)
	      fprintf(ion_inix[abs_index],"%13.8e  %13.8e  %13.8e  %13.8e\n",
		      indxp[lmain],hdp[lmain],indxp2[lmain],hdp2[lmain]);
	    fprintf(ion_inix[abs_index],"\n\n");
		
	    verteilq (yion,qionfi,dstep,ntemp,hdp,NDC,indxp);	
	    verteil (vpriony,NION,hdp2,NDC,indxp2);	
	    for (lmain=0; lmain<NDC; lmain++)
	      fprintf(ion_iniy[abs_index],"%13.8e  %13.8e  %13.8e  %13.8e\n",
		      indxp[lmain],hdp[lmain],indxp2[lmain],hdp2[lmain]);
	    fprintf(ion_iniy[abs_index],"\n\n");
	  }
	      
	  //gridextx = 1.1*fmax(xion_max,6*sxs[jmain]);
	  //gridexty = 1.1*fmax(yion_max,6*sys[jmain]);
	  //printf ("\n bunch number = %d, gridextx = %lf, gridexty = %lf \n",jmain,gridextx*1.e9,gridexty*1.e9);
	  
	  if(xion_max < fabs(xsion)+4*sxsion && yion_max < fabs(ysion)+4*sysion) {    
	    gridextx = 1.1*fmax(xion_max,10.*xel_max[jmain]);
	    gridexty = 1.1*fmax(yion_max,10.*yel_max[jmain]); 
	  }
	  else {
	    gridextx = 1.1*fmax(fabs(xsion)+5.*sxsion,10.*xel_max[jmain]);
	    gridexty = 1.1*fmax(fabs(ysion)+5.*sysion,10.*yel_max[jmain]);
	  }
	  
	  fprintf (grid_pr, "%ld  %d  %13.5e  %13.5e  %13.5e  %13.5e  %13.5e  %13.5e  %13.5e  %13.5e\n",
		   it,jmain,gridextx*1.e9,gridexty*1.e9,1.1e9*xion_max,1.1e9*yion_max,1.1e9*(fabs(xsion)+5.*sxsion),
		   1.1e9*(fabs(ysion)+5.*sysion),1.1e9*10.*xel_max[jmain],1.1e9*10.*yel_max[jmain]);
	      
	  grid_init_phys(grid,gridextx,gridexty);
	  
	  for(kmain=0; kmain<NELB; kmain++) { 
	    vprelx[kmain] = xel[jmain*NELB+kmain];
	    vprely[kmain] = yel[jmain*NELB+kmain];
	    vprelxp[kmain] = xpel[jmain*NELB+kmain];
	    vprelyp[kmain] = ypel[jmain*NELB+kmain];
	  }
	      
	  particles_distribute(grid,vprelx,vprely,NELB,1,qe0);
	  
	  switch(ifi) {
	  case 0:
	    particles_distribute(grid,xion,yion,ntemp,2,qion);
	    break;
	  case 1:
	    particles_distribute2(grid,xion,yion,ntemp,2,qionfi);
	    break;
	  }
	  //printf ("\n bunch number = %d, qe0 = %lf, qion = %lf \n",jmain,qe0,qion);
	  
	  /* calculate the forces */
	  
	  field_calculate(grid);
	  
	  /* move the ions - kick the beam electrons */
	  
	  particles_move(grid,epr,vprelx,vprely,vprelxp,vprelyp,NELB,1,tstep,fscale2);
	  particles_move_ion(grid,eel,xion,yion,xpion,ypion,mion,ntemp,2,tstep);
	  fprintf(ele_phase, "%ld  %ld  %ld  %lg\n ", it, jmain, ionoutgrid, (double)ionoutgrid/(double)ntemp*100.);
	  
	  //have to initialize mion[], a vector containing a flag for the ion type -i.e. mion[k] is the integer
	  //corresponding to which ion type the k-th ion is 
	  
	  for(kmain=0; kmain<NELB; kmain++) { 
	    xel[jmain*NELB+kmain] = vprelx[kmain] ;
	    yel[jmain*NELB+kmain] = vprely[kmain] ;
	    xpel[jmain*NELB+kmain] = vprelxp[kmain] ;
	    ypel[jmain*NELB+kmain] = vprelyp[kmain] ;
	  }
	  
	  
	  for (kmain=0; kmain<ntemp; kmain++) {
	    xion[kmain]+=xpion[kmain]*bsp*1e-9;
	    yion[kmain]+=ypion[kmain]*bsp*1e-9;
	    if(fabs(xion[kmain])>xion_max) xion_max=fabs(xion[kmain]);
	    if(fabs(yion[kmain])>yion_max) yion_max=fabs(yion[kmain]);
	  }
	      
	  if (dg_index==0)
	    fprintf(samp_ion[abs_index], "%ld  %ld  %lg  %lg  %lg  %lg  %ld  %lg  %lg  %lg  %lg\n", 
		    jmain, mion[10], xion[10], yion[10], 
		    xpion[10], ypion[10], mion[100], xion[100], yion[100], xpion[100], ypion[100]);
	  if (dg_index==1)
	    fprintf(samp2_ion[abs_index], "%ld  %ld  %lg  %lg  %lg  %lg  %ld  %lg  %lg  %lg  %lg\n", 
		    jmain, mion[10], xion[10], yion[10], 
		    xpion[10], ypion[10], mion[101], xion[101], yion[101], xpion[101], ypion[101]);

	  //        fprintf(ions_pr, "%ld  %ld  %lg  %lg  %lg  %lg\n", jmain, mion[42], xion[42], yion[42], xpion[42], ypion[42]);
	      
	}
	    
	/*  and  now the kick due to the wake-field is added
	    to the x' and y' of each proton in order to study the combined
	    effect ion + wake-field                                         */ 
	
	if (iwake_flag == 1) { 
	  forkickx = 0.;
	  forkicky = 0.;
	  //forkickz = 0.;
	  
	  switch(i_pipe) {
	    
	  case 1:
	  case 2:
	  case 3:
	    for (mmain=0; mmain<jmain; mmain++) { 
	      forkickx += wakefac*nele*xs[mmain]*wake_store[jmain-mmain];
	      forkicky += wakefac*nele*ys[mmain]*wake_store[jmain-mmain];
	      //forkickz += wakefac*nele*wakez_store[jmain-mmain];
	    }
	    break;
	    
	  default:
	    break;
	  }
	  
	  for (kmain=0; kmain<NELB; kmain++) { 
	    kick_xe = 0;
	    kick_ye = 0;
	    //kick_z = wakefac*alphaz*rshz*1.e6*nele;
	    switch(i_pipe) {
	      
	    case 1:
	    case 2:
	    case 3:
	      kick_xe += forkickx;
	      kick_ye += forkicky;
	      //kick_z += forkickz;
	      break;
	      
	    case 11:
	    case 22:
	    case 33:
	      for (mmain=0; mmain<jmain; mmain++) { 
		kick_xe += wakefac*nele*PI*PI/24*(xs[mmain]-xel[(jmain-1)*NELB+kmain])*wake_store[jmain-mmain];
		kick_ye += wakefac*nele*PI*PI/12*(ys[mmain]+0.5*yel[(jmain-1)*NELB+kmain])*wake_store[jmain-mmain];
		//kick_z += wakefac*nele*wakez_store[jmain-mmain];
	      }
	      break;
	      
	    default:
	      break;
	    }
	    
	    xpel[jmain*NELB+kmain] += kick_xe;
	    ypel[jmain*NELB+kmain] += kick_ye;
	    //dp[(jmain-1)*NELB+kmain] += kick_z;
	  }
	  
	}
	
      }     //closes the loop on bunches ->jmain
	
      /* loop over all the bunch particles starts, for the transformation of
	   their phase space coordinates all through the piece of ring
	   between two subsequent kicks...                                      */
	
      for(kmain=0; kmain<NELB*NBUN; kmain++) {
	xel0 = xel[kmain];
	xpel0 = xpel[kmain];          
	yel0 = yel[kmain];
	ypel0 = ypel[kmain];
	  
	/*
	  xx1 = x0;
	  xxp1 = xp0;
	  yy1 = ips0;
	  yyp1 = yp0;
	  
	  if(i_oct==1)
	  { epsx = (x0*x0 + xp0*xp0*betax0*betax0)/betax0;
	  epsy = (ips0*ips0 + yp0*yp0*betay0*betay0)/betay0;
	  
	  // Remember that epsx = 2*actx;  epsy = 2*acty   
	    
	  cxa = cos(2*PI*(htune + qpxx*epsx + qpxy*epsy)/(double)nkick);
	  sxa = sin(2*PI*(htune + qpxx*epsx + qpxy*epsy)/(double)nkick);
	  cya = cos(2*PI*(vtune + qpxy*epsx + qpyy*epsy)/(double)nkick);
	  sya = sin(2*PI*(vtune + qpxy*epsx + qpyy*epsy)/(double)nkick);
	  
	  }
	    
	  if(i_coupl==1)
	  { xxp1 -= key_xy*yy1;
	  yyp1 -= key_xy*xx1;
	  }
	*/
	if(i_lattice==0) {      
	  if(it%2==1) {		
	    /* horizontal coordinates transformation - linear rotation          */
	    xel[kmain] = coeff1 * cxya * xel0 + coeff2 * sxya * xpel0;
	    xpel[kmain] = -coeff4 * sxya * xel0 + coeff3 * cxya * xpel0;
	    
	    /* vertical coordinates transformation - linear rotation            */
	    
	    yel[kmain] = coeff3 * cxya * yel0 + coeff2 * sxya * ypel0;
	    ypel[kmain] = -coeff4 * sxya * yel0 + coeff1 * cxya * ypel0;
	  }
	  else {
	    xel[kmain] = coeff3 * cxya * xel0 + coeff2 * sxya * xpel0;
	    xpel[kmain] = -coeff4 * sxya * xel0 + coeff1 * cxya * xpel0;
	    yel[kmain] = coeff1 * cxya * yel0 + coeff2 * sxya * ypel0;
	    ypel[kmain] = -coeff4 * sxya * yel0 + coeff3 * cxya * ypel0;
	  }
	}
	else {
	  xel[kmain] = coeffga * (coeff1x * (cxa + coeff1xbis * sxa) * xel0 + coeff2x * sxa * xpel0);
	  xpel[kmain] = coeffga * (coeff4x * (-coeff4xbis * sxa + coeff4xter * cxa) * xel0 + coeff3x * (cxa - coeff3xbis * sxa) * xpel0);
	  yel[kmain] = coeffga * (coeff1y * (cya + coeff1ybis * sya) * yel0 + coeff2y * sya * ypel0);
	  ypel[kmain] = coeffga * (coeff4y * (-coeff4ybis * sya + coeff4yter * cya) * yel0 + coeff3y * (cya - coeff3ybis * sya) * ypel0);
	}
	  
	/*
	  if (kmain % 500==0 && (it-1)%n_diag==0)
	  fprintf(ion_phase,"%d \t%13.8e\t%13.8e\t%13.8e\t%13.8e\t%13.8e\t%13.8e\n",kmain,xel[kmain],xpel[kmain],yel[kmain],ypel[kmain]); 
	  }
	  
	  if (it %nkick==1) //(it % (n_diag/2)==0)
	  {
	  fprintf(track_pr,"\n");
	  fprintf(track_pr,"\n");
	  }
	  
	  
	*/
      }
      
      gammaprev = gammanext;
      
      fprintf(trainhdtl_pr,"\n");
      fprintf(trainhdtl_pr,"\n");
      fprintf(ions_pr, "\n\n");
      //fprintf(prot_phase,"\n");
      //fprintf(prot_phase,"\n");
      //fprintf(prot_ini,"\n");
      //fprintf(prot_ini,"\n");
      
      /*
	for (lmain=0;lmain<499;lmain++)
	{ enedis[lmain] /= qt_imp;
	fprintf (elwall_ene, "%13.8e  %13.8e\n",(double)lmain*enestep+enestep/2,enedis[lmain]);
	}
	fprintf (elwall_ene,"\n");
	fprintf (elwall_ene,"\n");
      */
      
      
      
    }
    
  }
  
  free_dvector(pss,0,nspe-1);
  free_dvector(am,0,nspe-1);
  free_dvector(crsec,0,nspe-1);
  free_dvector(fscale1,0,nspe-1);
  free_dvector(betarrayx,0,nstep);
  free_dvector(betarrayy,0,nstep);
  free_dvector(benergy,0,nstep);
  free_dvector(alfarrayx,0,nstep);
  free_dvector(alfarrayy,0,nstep);
  free_dvector(muarrayx,0,nstep);
  free_dvector(muarrayy,0,nstep);
  free_dvector(prob,0,nspe-1);
  
  fclose(centr_pr);
  fclose(trainhdtl_pr);
  fclose(ele_phase);
  for(lmain=0;lmain<10;lmain++) {  
    fclose(ion_inix[lmain]);
    fclose(ion_iniy[lmain]);
    fclose(samp_ion[lmain]);
    fclose(samp2_ion[lmain]);
  }
  fclose(ion_inphase);
  fclose(wakefield_pr);
  fclose(wakefield2_pr);
  fclose(sigma_pr);
  fclose(ions_pr);
  fclose(grid_pr);

  return EXIT_SUCCESS;

}

