/*
Spherically symmetric solver for the Einstein equations
coupled to hydrodynamics.
Copyright William East, 2009
*/

#define NUM_HYDRO_VAR 3
#define CONSTANT_PI 3.1415926535897932384626433
//=============================================================================
// application interface functions
//=============================================================================
#include <stdlib.h>
#include <stdio.h>
#include <math.h>
#include <mpi.h>
#include "pamr.h"
#include "amrd.h"
#include "fluidEvolve.h"
#include "initialData.h"
#include "gravityElliptic.h"
#include "fluxCorrection.h"
#include "ss_gr_hydro.h"
#include "blackHoleFinder.h"
#include "flux.h"

//=============================================================================
// id parameters 
//=============================================================================

real Pressure_central; //Central pressure of TOV initial star
real U_amp;  //Coordinate velocity profile amplitude

//Not a parameter, but derived from id parameters
real initialRadius; //radius of initial TOV star

//=============================================================================
// some convenient, "local" global variables
//=============================================================================

real aMax=1.0;
real TMax = 0.0;
real TGlobalMax = 0.0;
real rMax;
real *consVar_n, *consVar_np1, *consVar, *primVar_n, *primVar_np1, *primVar, *a_n, *a_np1, *a, *alpha_n, *alpha_np1, *alpha; // in AMRH n/np1
real *consVar_v, *primVar_v;
real *fluxCorr;
real *q1_n, *q1_np1, *q1, *q2_n, *q2_np1, *q2, *q3_n, *q3_np1, *q3, *p1_n, *p1_np1, *p1, *p2_n, *p2_np1, *p2, *p3_n, *p3_np1, *p3;
real *q1_v, *q2_v, *q3_v;
real *p1_v, *p2_v, *p3_v;
real *fcs_1, *fcs_2, *fcs_3;
real *phi_n, *phi_np1, *phi;
real *phi_res, *phi_lop, *phi_rhs;
real *a_res;
real *T_trace;
real *mask_c, *mask_v, *mask_mg;
real *wavg, *wavg_mg;

real *rVertex, *rCell;

int shape[3],shape_c[3],ghost_width[6],Nr,phys_bdy[6],numCells,g_rank,dim;
real base_bbox[6],bbox[6],dr,dt;
real t;
int g_L;

int q1_n_gfn, q2_n_gfn, q3_n_gfn, q1_np1_gfn, q2_np1_gfn, q3_np1_gfn; 
int q1_gfn, q2_gfn, q3_gfn;
int q1_v_gfn, q2_v_gfn, q3_v_gfn;
int fcs_1_gfn, fcs_2_gfn, fcs_3_gfn;
int p1_n_gfn, p2_n_gfn, p3_n_gfn, p1_np1_gfn, p2_np1_gfn, p3_np1_gfn; 
int p1_gfn, p2_gfn, p3_gfn;
int p1_v_gfn, p2_v_gfn, p3_v_gfn;
int a_n_gfn, a_np1_gfn;
int a_gfn;
int a_res_gfn;
int T_trace_gfn;
int alpha_n_gfn, alpha_np1_gfn;
int alpha_gfn;
int phi_n_gfn, phi_np1_gfn, phi_gfn;
int phi_res_gfn, phi_lop_gfn, phi_rhs_gfn;
int mask_mg_gfn, mask_v_gfn, mask_c_gfn;
int wavg_gfn, wavg_mg_gfn;

#if 0
int Krr, Kthth;
#endif

//=============================================================================
// call after variables have been defined
//=============================================================================
void set_gfns(void)
{
    if ((q1_n_gfn   = PAMR_get_gfn("q1",PAMR_AMRH,2))<0) AMRD_stop("set_gnfs error",0);
    if ((q1_np1_gfn = PAMR_get_gfn("q1",PAMR_AMRH,1))<0) AMRD_stop("set_gnfs error",0);
    if ((q1_gfn=PAMR_get_gfn("q1",PAMR_MGH,0))<0) AMRD_stop("set_gnfs error",0);
    if ((q1_v_gfn=PAMR_get_gfn("q1_v",PAMR_MGH,0))<0) AMRD_stop("set_gnfs error",0);

    if ((q2_n_gfn   = PAMR_get_gfn("q2",PAMR_AMRH,2))<0) AMRD_stop("set_gnfs error",0);
    if ((q2_np1_gfn = PAMR_get_gfn("q2",PAMR_AMRH,1))<0) AMRD_stop("set_gnfs error",0);
    if ((q2_gfn=PAMR_get_gfn("q2",PAMR_MGH,0))<0) AMRD_stop("set_gnfs error",0);
    if ((q2_v_gfn=PAMR_get_gfn("q2_v",PAMR_MGH,0))<0) AMRD_stop("set_gnfs error",0);
        
    if ((q3_n_gfn   = PAMR_get_gfn("q3",PAMR_AMRH,2))<0) AMRD_stop("set_gnfs error",0);
    if ((q3_np1_gfn = PAMR_get_gfn("q3",PAMR_AMRH,1))<0) AMRD_stop("set_gnfs error",0);
    if ((q3_gfn=PAMR_get_gfn("q3",PAMR_MGH,0))<0) AMRD_stop("set_gnfs error",0);
    if ((q3_v_gfn=PAMR_get_gfn("q3_v",PAMR_MGH,0))<0) AMRD_stop("set_gnfs error",0);
    
 
    if ((fcs_1_gfn= PAMR_get_gfn("q1_fcs",PAMR_AMRH,1))<0) AMRD_stop("set_gnfs error",0);
    if ((fcs_2_gfn= PAMR_get_gfn("q2_fcs",PAMR_AMRH,1))<0) AMRD_stop("set_gnfs error",0);
    if ((fcs_3_gfn= PAMR_get_gfn("q3_fcs",PAMR_AMRH,1))<0) AMRD_stop("set_gnfs error",0);

    if ((p1_n_gfn   = PAMR_get_gfn("p1",PAMR_AMRH,2))<0) AMRD_stop("set_gnfs error",0);
    if ((p1_np1_gfn = PAMR_get_gfn("p1",PAMR_AMRH,1))<0) AMRD_stop("set_gnfs error",0);
    if ((p1_gfn=PAMR_get_gfn("p1",PAMR_MGH,0))<0) AMRD_stop("set_gnfs error",0);
    if ((p1_v_gfn=PAMR_get_gfn("p1_v",PAMR_MGH,0))<0) AMRD_stop("set_gnfs error",0);
    
    if ((p2_n_gfn   = PAMR_get_gfn("p2",PAMR_AMRH,2))<0) AMRD_stop("set_gnfs error",0);
    if ((p2_np1_gfn = PAMR_get_gfn("p2",PAMR_AMRH,1))<0) AMRD_stop("set_gnfs error",0);
    if ((p2_gfn=PAMR_get_gfn("p2",PAMR_MGH,0))<0) AMRD_stop("set_gnfs error",0);
    if ((p2_v_gfn=PAMR_get_gfn("p2_v",PAMR_MGH,0))<0) AMRD_stop("set_gnfs error",0);
    
    if ((p3_n_gfn   = PAMR_get_gfn("p3",PAMR_AMRH,2))<0) AMRD_stop("set_gnfs error",0);
    if ((p3_np1_gfn = PAMR_get_gfn("p3",PAMR_AMRH,1))<0) AMRD_stop("set_gnfs error",0);
    if ((p3_gfn=PAMR_get_gfn("p3",PAMR_MGH,0))<0) AMRD_stop("set_gnfs error",0);
    if ((p3_v_gfn=PAMR_get_gfn("p3_v",PAMR_MGH,0))<0) AMRD_stop("set_gnfs error",0);
    
    if ((a_n_gfn   = PAMR_get_gfn("a_v",PAMR_AMRH,2))<0) AMRD_stop("set_gnfs error",0);
    if ((a_np1_gfn = PAMR_get_gfn("a_v",PAMR_AMRH,1))<0) AMRD_stop("set_gnfs error",0);
    if ((a_gfn=PAMR_get_gfn("a_v",PAMR_MGH,0))<0) AMRD_stop("set_gnfs error",0);
   
    if ((a_res_gfn=PAMR_get_gfn("a_res",PAMR_AMRH,1))<0) AMRD_stop("set_gnfs error",0);
    
    if ((T_trace_gfn=PAMR_get_gfn("T",PAMR_AMRH,1))<0) AMRD_stop("set_gnfs error",0);

    if ((alpha_n_gfn   = PAMR_get_gfn("alpha_v",PAMR_AMRH,2))<0) AMRD_stop("set_gnfs error",0);
    if ((alpha_np1_gfn = PAMR_get_gfn("alpha_v",PAMR_AMRH,1))<0) AMRD_stop("set_gnfs error",0);
    if ((alpha_gfn=PAMR_get_gfn("alpha_v",PAMR_MGH,0))<0) AMRD_stop("set_gnfs error",0);

    if ((phi_n_gfn   = PAMR_get_gfn("phi_v",PAMR_AMRH,2))<0) AMRD_stop("set_gnfs error",0);
    if ((phi_np1_gfn = PAMR_get_gfn("phi_v",PAMR_AMRH,1))<0) AMRD_stop("set_gnfs error",0);
    if ((phi_gfn=PAMR_get_gfn("phi_v",PAMR_MGH,0))<0) AMRD_stop("set_gnfs error",0);
    

    if ((phi_res_gfn=PAMR_get_gfn("phi_res_v",PAMR_MGH,0))<0) AMRD_stop("set_gnfs error",0);
    if ((phi_lop_gfn=PAMR_get_gfn("phi_lop_v",PAMR_MGH,0))<0) AMRD_stop("set_gnfs error",0);
    if ((phi_rhs_gfn=PAMR_get_gfn("phi_rhs_v",PAMR_MGH,0))<0) AMRD_stop("set_gnfs error",0);
    if ((mask_mg_gfn=PAMR_get_gfn("cmask_v",PAMR_MGH,0))<0) AMRD_stop("set_gnfs error",0);
    if ((mask_v_gfn=PAMR_get_gfn("cmask_v",PAMR_AMRH,1))<0) AMRD_stop("set_gnfs error",0);
    if ((mask_c_gfn=PAMR_get_gfn("cmask",PAMR_AMRH,1))<0) AMRD_stop("set_gnfs error",0);
   
    #if 0 
    if ((Krr_gfn=PAMR_get_gfn("Krr_v",PAMR_MGH,0))<0) AMRD_stop("set_gnfs error",0);
    if ((Kthth_gfn=PAMR_get_gfn("Kthth_v",PAMR_MGH,0))<0) AMRD_stop("set_gnfs error",0);
    #endif 1
    if (AMRD_num_inject_wavg_vars) if ((wavg_mg_gfn = PAMR_get_gfn("wavg",PAMR_MGH,0))<0) AMRD_stop("set_gnfs error",0);
    if (AMRD_num_inject_wavg_vars) if ((wavg_gfn    = PAMR_get_gfn("wavg",PAMR_AMRH,1))<0) AMRD_stop("set_gnfs error",0);
}

//=============================================================================
// call with valid iter to set up globals:
//=============================================================================
void ldptr(void)
{
   int ngfs;
   real dx0[3];
   real *x0[3], *x0_c[3];
   real *gfs[PAMR_MAX_GFNS];
   static int first=1;

   if (first) 
   {
      first=0; 
      set_gfns();
      PAMR_get_global_bbox(base_bbox);
   }

    if (!(PAMR_get_g_attribs(&g_rank,&dim,shape,shape_c,bbox,ghost_width,&t,&ngfs,x0,x0_c,gfs))) 
      AMRD_stop("ldptr: PAMR_get_g_attribs failed\n","");  

   PAMR_get_g_level(&g_L);
   PAMR_get_dxdt(g_L,dx0,&dt);
   dr=dx0[0];

   if ((bbox[0]-base_bbox[0])<dr/2) phys_bdy[0]=1; else phys_bdy[0]=0;
   if ((base_bbox[1]-bbox[1])<dr/2) phys_bdy[1]=1; else phys_bdy[1]=0;
   Nr=shape[0];
   numCells=shape_c[0];

   PAMR_get_g_x(x0);
   PAMR_get_g_x_c(x0_c);

   rVertex=x0[0];
   rCell=x0_c[0];

   q1_n   = gfs[q1_n_gfn-1];
   q1_np1 = gfs[q1_np1_gfn-1];
   q1 = gfs[q1_gfn-1];
   q1_v = gfs[q1_v_gfn-1];
   q2_n   = gfs[q2_n_gfn-1];
   q2_np1 = gfs[q2_np1_gfn-1];
   q2 = gfs[q2_gfn-1];
   q2_v = gfs[q2_v_gfn-1];
   q3_n   = gfs[q3_n_gfn-1];
   q3_np1 = gfs[q3_np1_gfn-1];
   q3 = gfs[q3_gfn-1];
   q3_v = gfs[q3_v_gfn-1];

   fcs_1 = gfs[fcs_1_gfn-1];
   fcs_2 = gfs[fcs_2_gfn-1];
   fcs_3 = gfs[fcs_3_gfn-1];

   p1_n   = gfs[p1_n_gfn-1];
   p1_np1 = gfs[p1_np1_gfn-1];
   p1 = gfs[p1_gfn-1];
   p1_v = gfs[p1_v_gfn-1];
   p2_n   = gfs[p2_n_gfn-1];
   p2_np1 = gfs[p2_np1_gfn-1];
   p2 = gfs[p2_gfn-1];
   p2_v = gfs[p2_v_gfn-1];
   p3_n   = gfs[p3_n_gfn-1];
   p3_np1 = gfs[p3_np1_gfn-1];
   p3 = gfs[p3_gfn-1];
   p3_v = gfs[p3_v_gfn-1];
   a_n   = gfs[a_n_gfn-1];
   a_np1 = gfs[a_np1_gfn-1];
   a   = gfs[a_gfn-1];
   a_res = gfs[a_res_gfn-1];
   T_trace = gfs[T_trace_gfn-1];
   alpha_n   = gfs[alpha_n_gfn-1];
   alpha_np1 = gfs[alpha_np1_gfn-1];
   alpha = gfs[alpha_gfn-1];
   phi_n = gfs[phi_n_gfn-1];
   phi_np1 = gfs[phi_np1_gfn-1];
   phi = gfs[phi_gfn-1];

   phi_res = gfs[phi_res_gfn-1];
   phi_lop = gfs[phi_lop_gfn-1];
   phi_rhs = gfs[phi_rhs_gfn-1];

   mask_mg = gfs[mask_mg_gfn-1];
   mask_v = gfs[mask_v_gfn-1];
   mask_c = gfs[mask_c_gfn-1];
  
   #if 0
   Krr = gfs[Krr_gtn-1];
   Kthth = gfs[Kthth_gtn-1];
   #endif	
 
   if (AMRD_num_inject_wavg_vars) wavg = gfs[wavg_gfn-1];
   if (AMRD_num_inject_wavg_vars) wavg_mg = gfs[wavg_mg_gfn-1];

   consVar_n = malloc(numCells*NUM_HYDRO_VAR*sizeof(real));
   consVar_np1 = malloc(numCells*NUM_HYDRO_VAR*sizeof(real));
   consVar = malloc(numCells*NUM_HYDRO_VAR*sizeof(real));
   consVar_v = malloc(Nr*NUM_HYDRO_VAR*sizeof(real));
   
   fluxCorr = malloc(numCells*NUM_HYDRO_VAR*sizeof(real));

   primVar_n = malloc(numCells*NUM_HYDRO_VAR*sizeof(real));
   primVar_np1 = malloc(numCells*NUM_HYDRO_VAR*sizeof(real));
   primVar = malloc(numCells*NUM_HYDRO_VAR*sizeof(real));
   primVar_v = malloc(Nr*NUM_HYDRO_VAR*sizeof(real));

    /*Copy hydro variables into 'vector' format.
     Note that this memory must be deallocated by call to deallocVec 
   */
   int i;
   for(i=0; i<numCells; i++){	

	if(q1_n!=NULL && q2_n!=NULL && q3_n!=NULL){
		consVar_n[i*NUM_HYDRO_VAR+0] = q1_n[i];
		consVar_n[i*NUM_HYDRO_VAR+1] = q2_n[i];
		consVar_n[i*NUM_HYDRO_VAR+2] = q3_n[i];
	}

	if(q1_np1!=NULL && q2_np1!=NULL && q3_np1!=NULL){
		consVar_np1[i*NUM_HYDRO_VAR+0] = q1_np1[i];
		consVar_np1[i*NUM_HYDRO_VAR+1] = q2_np1[i];
		consVar_np1[i*NUM_HYDRO_VAR+2] = q3_np1[i];
	}

	if(q1!=NULL && q2!=NULL && q3!=NULL){
		consVar[i*NUM_HYDRO_VAR+0] = q1[i];
        	consVar[i*NUM_HYDRO_VAR+1] = q2[i];
		consVar[i*NUM_HYDRO_VAR+2] = q3[i];
	}

	if(fcs_1!=NULL && fcs_2!=NULL && fcs_3!=NULL){
		fluxCorr[i*NUM_HYDRO_VAR+0] = fcs_1[i];
        	fluxCorr[i*NUM_HYDRO_VAR+1] = fcs_2[i];
		fluxCorr[i*NUM_HYDRO_VAR+2] = fcs_3[i];
	}
	
	if(p1_n!=NULL && p2_n!=NULL && p3_n!=NULL){
		primVar_n[i*NUM_HYDRO_VAR+0] = p1_n[i];
		primVar_n[i*NUM_HYDRO_VAR+1] = p2_n[i];
		primVar_n[i*NUM_HYDRO_VAR+2] = p3_n[i];
	}

	if(p1_np1!=NULL && p2_np1!=NULL && p3_np1!=NULL){
		primVar_np1[i*NUM_HYDRO_VAR+0] = p1_np1[i];
		primVar_np1[i*NUM_HYDRO_VAR+1] = p2_np1[i];
		primVar_np1[i*NUM_HYDRO_VAR+2] = p3_np1[i];
	}	

	if(p1!=NULL && p2!=NULL && p3!=NULL){
		primVar[i*NUM_HYDRO_VAR+0] = p1[i];
		primVar[i*NUM_HYDRO_VAR+1] = p2[i];
		primVar[i*NUM_HYDRO_VAR+2] = p3[i];
	}
   }

   for(i=0; i<Nr; i++){
	
	if(q1_v!=NULL && q2_v!=NULL && q3_v!=NULL){
		consVar_v[i*NUM_HYDRO_VAR+0] = q1_v[i];
		consVar_v[i*NUM_HYDRO_VAR+1] = q2_v[i];
		consVar_v[i*NUM_HYDRO_VAR+2] = q3_v[i];
	}
	if(p1_v!=NULL && p2_v!=NULL && p3_v!=NULL){
		primVar_v[i*NUM_HYDRO_VAR+0] = p1_v[i];
		primVar_v[i*NUM_HYDRO_VAR+1] = p2_v[i];
		primVar_v[i*NUM_HYDRO_VAR+2] = p3_v[i];
	}
   }

}

void deallocVec(void){
   
   int i;
   for(i=0; i<numCells; i++){

	if(q1_n!=NULL && q2_n!=NULL && q3_n!=NULL){
		q1_n[i] = consVar_n[i*NUM_HYDRO_VAR+0];
		q2_n[i] = consVar_n[i*NUM_HYDRO_VAR+1];
		q3_n[i] = consVar_n[i*NUM_HYDRO_VAR+2];
	}

	if(q1_np1!=NULL && q2_np1!=NULL && q3_np1!=NULL){
		q1_np1[i] = consVar_np1[i*NUM_HYDRO_VAR+0];
		q2_np1[i] = consVar_np1[i*NUM_HYDRO_VAR+1];
		q3_np1[i] = consVar_np1[i*NUM_HYDRO_VAR+2];
	}
	
	if(fcs_1!=NULL && fcs_2!=NULL && fcs_3!=NULL){
		fcs_1[i] = fluxCorr[i*NUM_HYDRO_VAR+0];
		fcs_2[i] = fluxCorr[i*NUM_HYDRO_VAR+1];
		fcs_3[i] = fluxCorr[i*NUM_HYDRO_VAR+2];
	}
/*
	if(q1!=NULL && q2!=NULL && q3!=NULL && q1!=q1_n && q1!=q1_np1){
		q1[i] = consVar[i*NUM_HYDRO_VAR+0];
		q2[i] = consVar[i*NUM_HYDRO_VAR+1];
		q3[i] = consVar[i*NUM_HYDRO_VAR+2];
	}
*/
	if(p1_n!=NULL && p2_n!=NULL && p3_n!=NULL){
		p1_n[i] = primVar_n[i*NUM_HYDRO_VAR+0];
		p2_n[i] = primVar_n[i*NUM_HYDRO_VAR+1];
		p3_n[i] = primVar_n[i*NUM_HYDRO_VAR+2];
	}

	if(p1_np1!=NULL && p2_np1!=NULL && p3_np1!=NULL){
		p1_np1[i] = primVar_np1[i*NUM_HYDRO_VAR+0];
		p2_np1[i] = primVar_np1[i*NUM_HYDRO_VAR+1];
		p3_np1[i] = primVar_np1[i*NUM_HYDRO_VAR+2];
   	}
/*
	if(p1!=NULL && p2!=NULL && p3!=NULL){
		p1[i] = primVar[i*NUM_HYDRO_VAR+0];
		p2[i] = primVar[i*NUM_HYDRO_VAR+1];
		p3[i] = primVar[i*NUM_HYDRO_VAR+2];
	}
*/  
 }
	
	free(consVar_v);
	consVar_v = NULL;
	free(primVar_v);
	primVar_v = NULL;	
	free(consVar_n);
	consVar_n = NULL;
	free(consVar_np1);
	consVar_np1 = NULL;
	free(consVar);
	consVar = NULL;
	free(primVar_n);
	primVar_n = NULL;
	free(primVar_np1);
	primVar_np1 = NULL;
	free(primVar);
	primVar = NULL;
	free(fluxCorr);
	fluxCorr = NULL;
}

//=============================================================================
// utility routines
//=============================================================================
void const_f(real *f, real c, int length)
{
   int i;

   for (i=0; i<length; i++) f[i]=c;
}

void zero(real *f, int length)
{
   const_f(f,0, length);
}

/*
Calculate spherical volume element.  For use with wavg transfer
ldptr() should have already been called.
*/
void calculateVolumeElement(real *volume){
	int i;
	for(i=0; i<numCells; i++){
		volume[i] = 4.0/3.0*CONSTANT_PI*(rVertex[i+1]*rVertex[i+1]*rVertex[i+1]-rVertex[i]*rVertex[i]*rVertex[i]);
	}
}

/*
Set the regrid_transfer parameter of a variable
*/
int setVarRegridTransfer(char *varName, int regridTransfer){
   
   int in_amrh;
   int in_mgh;
   int num_tl;
   int amr_inject;
   int amr_interp;
   int amr_bdy_interp;
   int amr_sync;
   int mg_inject;
   int mg_interp;
   int mg_sync;
   int mg_noinj_to_amr;
   int regrid_transfer;
   int c_to_v;
   int v_to_c;
   int phys_bdy_type[6];

   if( !PAMR_get_var_attribs(varName, &in_amrh, &in_mgh, &num_tl, &amr_inject,
                   &amr_interp, &amr_bdy_interp, &amr_sync, &mg_inject, 
                   &mg_interp, &mg_sync, &mg_noinj_to_amr, &regrid_transfer, 
                   &c_to_v, &v_to_c, phys_bdy_type)) return 0;

   regrid_transfer = regridTransfer;

   if( !PAMR_set_var_attribs(varName, in_amrh, in_mgh, num_tl, amr_inject,
                   amr_interp, amr_bdy_interp, amr_sync, mg_inject, 
                   mg_interp, mg_sync, mg_noinj_to_amr, regrid_transfer, 
		   c_to_v, v_to_c, phys_bdy_type)) return 0;

   return 1;
}

/*   
Set the regrid_transfer parameter of a variable
*/
int setMGInject(char *varName, int mgInject){
   
   int in_amrh;
   int in_mgh;
   int num_tl;
   int amr_inject;
   int amr_interp;
   int amr_bdy_interp;
   int amr_sync;
   int mg_inject;
   int mg_interp;
   int mg_sync;
   int mg_noinj_to_amr;
   int regrid_transfer;
   int c_to_v;
   int v_to_c;
   int phys_bdy_type[6];

   if( !PAMR_get_var_attribs(varName, &in_amrh, &in_mgh, &num_tl, &amr_inject,
                   &amr_interp, &amr_bdy_interp, &amr_sync, &mg_inject, 
                   &mg_interp, &mg_sync, &mg_noinj_to_amr, &regrid_transfer, 
                   &c_to_v, &v_to_c, phys_bdy_type)) return 0;

   mg_inject = mgInject;

   if( !PAMR_set_var_attribs(varName, in_amrh, in_mgh, num_tl, amr_inject,
                   amr_interp, amr_bdy_interp, amr_sync, mg_inject, 
                   mg_interp, mg_sync, mg_noinj_to_amr, regrid_transfer, 
		   c_to_v, v_to_c, phys_bdy_type)) return 0;

   return 1;
   
}

//=============================================================================
// Routines required by amrd:
//=============================================================================

//=============================================================================
// Returns 0 to use default mechanism, or is expected to calculate
// the correct initial hierarchy and return 1:
//=============================================================================
int ssgrhydro_id(void)
{
   return 0;
}

//=============================================================================
// Sets custom parameters, variables, etc. Split up into two segments,
// one called before the pamr context is initialized and standard
// parameters are read, and the other afterwards
//=============================================================================
void ssgrhydro_var_pre_init(char *pfile)
{
   return;
}

void ssgrhydro_var_post_init(char *pfile)
{
   setVarRegridTransfer("q1",PAMR_FIRST_ORDER_CONS);
   setVarRegridTransfer("q2",PAMR_FIRST_ORDER_CONS);
   setVarRegridTransfer("q3",PAMR_FIRST_ORDER_CONS);

   setVarRegridTransfer("p1",PAMR_FIRST_ORDER_CONS);
   setVarRegridTransfer("p2",PAMR_FIRST_ORDER_CONS);
   setVarRegridTransfer("p3",PAMR_FIRST_ORDER_CONS);
   
   int i,j;
   char buf[64];

   if (my_rank==0)
   {
      printf("===================================================================\n");
      printf("Reading SS hydro-GR parameters:\n\n");
   }

   Pressure_central = 22.0e-4;
   U_amp = 0.4; 

   AMRD_real_param(pfile,"Pressure_c",&Pressure_central,1);
   AMRD_real_param(pfile,"U_amp",&U_amp,1);
   
   int eos_flag = 0;
   real gamma_poly = 2.0;
   real Ye = 0.2;
   AMRD_int_param(pfile, "eos_flag", &eos_flag,1);
   AMRD_real_param(pfile, "gamma_poly", &gamma_poly,1);
   AMRD_real_param(pfile, "Ye", &Ye,1);
  
   # include "equationOfState.h" 
   if(eos_flag==GAMMA_LAW){
	setEquationOfState(eos_flag, gamma_poly);
   } else if(eos_flag==ADIABATIC){
	setEquationOfState(eos_flag, gamma_poly);
   } else if(eos_flag==SHEN_EOS){
	setEquationOfState(eos_flag, Ye);
   }

   if (my_rank==0) printf("===================================================================\n");
   return;
}

//=============================================================================
// Sets all variables to their 'zero' values:
//=============================================================================
void ssgrhydro_AMRH_var_clear(void)
{ 
   ldptr();
   zero(q1_n, numCells); 
   zero(q2_n, numCells); 
   zero(q3_n, numCells); 
   zero(p1_n, numCells); 
   zero(p2_n, numCells); 
   zero(p3_n, numCells); 
   zero(a_n, Nr);
   zero(alpha_n, Nr);
   zero(phi_n, Nr);
   zero(q1_np1, numCells); 
   zero(q2_np1, numCells); 
   zero(q3_np1, numCells); 
   zero(p1_np1, numCells); 
   zero(p2_np1, numCells); 
   zero(p3_np1, numCells);
   zero(a_np1, Nr);
   zero(alpha_np1, Nr);
   zero(phi_np1, Nr);


   deallocVec();
   return;
}

//=============================================================================
// Initial data for free fields: (at tn=2)
//
// currently, we only allow for time-symmetric initial data
//=============================================================================
void ssgrhydro_free_data(void)
{
   ldptr();
   calculateVolumeElement(wavg);

   if(phys_bdy[0] && phys_bdy[1]){
	initialRadius = getInitialData(Pressure_central, U_amp, consVar_n, primVar_n, a_n, alpha_n, rVertex, rCell, Nr);
   } else{
	printf("Partial domain initial data.\n");
	initialRadius = getInitialDataPartialDomain(Pressure_central, U_amp, base_bbox[1], consVar_n, primVar_n, a_n, alpha_n, rVertex, rCell, Nr);
   }
   rMax = initialRadius;

   deallocVec();

   return;
}  

//=============================================================================
// Initial constraint data --- called after each MG iteration.
//=============================================================================
void ssgrhydro_t0_cnst_data(void)
{
   return;
}

//=============================================================================
// Calculations prior to saving info to disk.
//
// NOTE: at this point, the time sequence is: n,nm1,np1
//=============================================================================
void ssgrhydro_pre_io_calc(void)
{
   return;
}

//=============================================================================
// Returns some norm of the residual for the evolution variables,
// called after an evolution iteration.
// We're using an explicit scheme to solve for q, hence return 0
//=============================================================================
real ssgrhydro_evo_residual(void)
{
   return 0;
}

//=============================================================================
// Returns some norm of the residual for the MG variables, *AND* 
// stores the point-wise residual in "f_res" for each MG variable "f" (for
// computing new RHS's)
//=============================================================================
real ssgrhydro_MG_residual(void)
{
   real norm;   

   ldptr();

   residualPhi(phi, phi_rhs, phys_bdy, consVar_v, primVar_v, a, rVertex, mask_mg, Nr, &norm, phi_res);
   
   deallocVec();
   
   return norm;
}

//=============================================================================
// Performs 1 iteration of the evolution equations 
//=============================================================================
void ssgrhydro_evolve(int iter, int *ifc_mask)
{
   ldptr();
   
   timeStep(iter, phys_bdy[0],phys_bdy[1], dt, Nr, rVertex, consVar_n, consVar_np1, primVar_n, primVar_np1, a_n, a_np1, phi_n, phi_np1, fluxCorr, ifc_mask);

	
   //Update T_trace
   if(iter==3) {
	getStressEnergyTraceArray(T_trace,primVar_np1,numCells);
	int i;
	for(i=0; i<numCells; i++) {
		if(T_trace[i]>TMax) TMax=T_trace[i];
	}
   	if(TMax>TGlobalMax) TGlobalMax=TMax;
   }	

  //Look for black hole formation 
  if(iter==3){	
	real rLocalMax;
	real aLocalMax;
	if(findBlackHole(a_np1, rVertex, Nr, &rLocalMax, &aLocalMax)) AMRD_stop("Black hole detected",0);
	if(aLocalMax>aMax) {
		aMax = aLocalMax;
		rMax = rLocalMax;
	}
	if(phys_bdy[0]==1 && phys_bdy[1] ==1 && g_L>1){

		if(rMax>(2.0*initialRadius+3.0*dr) && aMax<1.4 && TMax<0.9*TGlobalMax) {
			char output[80];
			snprintf(output, 79, "No BH found Tmax=%E\n",TGlobalMax);
			outputString(output);
			AMRD_stop("No Black hole formation",0);
		} else if(rMax<initialRadius){
			initialRadius = rMax;
		}
		aMax=aLocalMax;
		rMax=rLocalMax;
		TMax=0.0;
   	}
	//if(aLocalMax>2.0) printf("BH at r=%E a=%E\n",aLocalMax,rLocalMax);
   }

   deallocVec();

}

//=============================================================================
// sets excision mask (NO ITERATOR, SO DON'T LOAD POINTERS!!!)
//=============================================================================
void ssgrhydro_fill_ex_mask(real *mask, real *mask_c, int dim, int *shape, int *shape_c, real *bbox, real excised)
{
}

//=============================================================================
void ssgrhydro_fill_bh_bboxes(real *bbox, int *num, int max_num)
{
}

//=============================================================================
void ssgrhydro_post_tstep(int L)
{
}

//=============================================================================
// Performs 1 relaxation sweep of the MG equations, and returns an estimate
// of the norm of the residual.
//=============================================================================
real ssgrhydro_MG_relax(void)
{
   real norm;
   ldptr();
   
   relaxPhi(phi, phi_rhs, phys_bdy, consVar_v, primVar_v, a, rVertex, mask_mg, Nr, &norm);

   deallocVec();
   
   return norm;
}

//=============================================================================
// Computes the differential operator L acting on the current grid,
// in a region specified by the grid function "mask". Stores the result
// in "f_lop" for each MG variable "f"
//=============================================================================
void ssgrhydro_L_op(void)
{
   ldptr();
   LPhi(phi, phys_bdy, consVar_v, primVar_v, a, rVertex, mask_mg, Nr, phi_lop);
   deallocVec();
}

//=============================================================================
// Called after calculating the TRE for all variables
//=============================================================================
void ssgrhydro_scale_tre(void)
{

   ldptr();
   int i,j;
   
   //Hack to prevent refinement at initial time
   if(t<2*dt) {
   	for(i=0; i<Nr; i++) (AMRD_f_tre[0])[i] = 0.0;
   }
/*   
   //Make tre monotonic and positive
   AMRD_f_tre[0][0] = fabs(AMRD_f_tre[0][0]);
   for(i=1; i<Nr; i++) {
   	AMRD_f_tre[0][i] = fabs(AMRD_f_tre[0][i]);
	j = i-1;
	while( (AMRD_f_tre[0])[i] >(AMRD_f_tre[0])[j] && j>=0){
		(AMRD_f_tre[0])[j]=(AMRD_f_tre[0])[i];
		j--;
	}
   }
*/
   deallocVec();

}
//=============================================================================
// post-regrid initialization of constant functions
//=============================================================================
void ssgrhydro_post_regrid(void)
{
 
   ldptr();
   calculateVolumeElement(wavg);
   deallocVec();
   int i;
   for(i=0; i<Nr; i++) {
	if(a_n[i]==0 || a_np1[i]==0) printf("zero a %E %E after regrid at i=%d r=%E\n\n\n\n\n",a_n[i],a_np1[i],i,rVertex[i]);
	}
}

//===============================================================================
// Sets flux correction variables to their 'zero' values.
// typeflag = 0:  corresponds to negative sign flag in mask. zero cells of type A
// typeflag = 1:  corresponds to positive sign flag in mask. zero cells of type B
//===============================================================================
void ssgrhyrdo_fcs_var_clear(int type_flag, int *ifc_mask){

	ldptr();
	clearFluxCorrection(type_flag, fluxCorr, ifc_mask, numCells);
	deallocVec();
}

//==============================================================================
// A hook function to apply the flux correction.
//==============================================================================
void ssgrhyrdo_flux_correct(void)
{
  ldptr();
  applyFluxCorrection(consVar_n, fluxCorr, rVertex, numCells);
  
  real a_cell[numCells];
  int i;
  for(i=0; i<numCells; i++) a_cell[i] = 0.5*(a_n[i]+a_n[i+1]);
  getPrimitiveArray(consVar_n, primVar_n, a_cell, numCells);
  deallocVec();
}

//=============================================================================
// post flux correction hook function
//=============================================================================
void ssgrhyrdo_post_flux_correct(void)
{
  return;
}


//=============================================================================
int main(int argc, char **argv)
{
   //Set flux correction hook functions
   amrd_set_app_fcs_var_clear_hook(ssgrhyrdo_fcs_var_clear);
   amrd_set_app_flux_correct_hook(ssgrhyrdo_flux_correct);
   amrd_set_app_post_flux_correct_hook(ssgrhyrdo_post_flux_correct);
   
   amrd(argc,argv,&ssgrhydro_id,&ssgrhydro_var_pre_init,
        &ssgrhydro_var_post_init, &ssgrhydro_AMRH_var_clear,
        &ssgrhydro_free_data, &ssgrhydro_t0_cnst_data,
        &ssgrhydro_evo_residual, &ssgrhydro_MG_residual,
        &ssgrhydro_evolve, &ssgrhydro_MG_relax, &ssgrhydro_L_op, 
        &ssgrhydro_pre_io_calc, &ssgrhydro_scale_tre, 
        &ssgrhydro_post_regrid, &ssgrhydro_post_tstep,
        &ssgrhydro_fill_ex_mask, &ssgrhydro_fill_bh_bboxes);
}

