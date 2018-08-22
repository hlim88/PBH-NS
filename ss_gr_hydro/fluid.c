/* This code was generated by rnpl, a numerical programming language */
/* copyright (c) 1994-1998 by Robert L. Marsa and Matthew W. Choptuik */

#include <stdlib.h>
#include <stdio.h>
#include <math.h>
#include <string.h>
#include <time.h>
#include <signal.h>
#include <sdf.h>
#include <librnpl.h>

#ifndef TRUE
#define TRUE 1
#endif
#ifndef FALSE
#define FALSE 0
#endif

void swap_levels(double **q1_np1,double **q1_n,double **q2_np1,double **q2_n,double **q3_np1,double **q3_n,double **p1_np1,double **p1_n,double **p2_np1,double **p2_n,double **p3_np1,double **p3_n)
{
  int *tmpi;
  double *tmpr;

  tmpr=*q1_n;
  *q1_n=*q1_np1;
  *q1_np1=tmpr;
  tmpr=*q2_n;
  *q2_n=*q2_np1;
  *q2_np1=tmpr;
  tmpr=*q3_n;
  *q3_n=*q3_np1;
  *q3_np1=tmpr;
  tmpr=*p1_n;
  *p1_n=*p1_np1;
  *p1_np1=tmpr;
  tmpr=*p2_n;
  *p2_n=*p2_np1;
  *p2_np1=tmpr;
  tmpr=*p3_n;
  *p3_n=*p3_np1;
  *p3_np1=tmpr;
}

double calc_resid(int g1_Nr)
{
  double rnplu;
  return 0.0;
}

void take_step(int *rnpldone,double *q1_np1,double *q1_n,double *q2_np1,double *q2_n,double *q3_np1,double *q3_n,double *p1_np1,double *p1_n,double *p2_np1,double *p2_n,double *p3_np1,double *p3_n,double *source1_n,double *source2_n,double *source3_n,double *flux1_n,double *flux2_n,double *flux3_n,double *fluxAlt1_n,double *fluxAlt2_n,double *fluxAlt3_n,double *alpha_n,double *a_n,double *test_n,int g1_Nr,double *r,double dt,double t)
{

  rungeKuttaUpdater(rnpldone,q1_np1,q1_n,q2_np1,q2_n,q3_np1,q3_n,p1_np1,p1_n,p2_np1,p2_n,p3_np1,p3_n,source1_n,source2_n,source3_n,flux1_n,flux2_n,flux3_n,fluxAlt1_n,fluxAlt2_n,fluxAlt3_n,alpha_n,a_n,test_n,g1_Nr,r,dt);
}

int update_gfuncs(int *rnpldone, int maxstep, double epsiter,double *q1_np1,double *q1_n,double *q2_np1,double *q2_n,double *q3_np1,double *q3_n,double *p1_np1,double *p1_n,double *p2_np1,double *p2_n,double *p3_np1,double *p3_n,double *source1_n,double *source2_n,double *source3_n,double *flux1_n,double *flux2_n,double *flux3_n,double *fluxAlt1_n,double *fluxAlt2_n,double *fluxAlt3_n,double *alpha_n,double *a_n,double *test_n,int g1_Nr,double *r,double dt,double t)
{
  take_step(rnpldone,q1_np1,q1_n,q2_np1,q2_n,q3_np1,q3_n,p1_np1,p1_n,p2_np1,p2_n,p3_np1,p3_n,source1_n,source2_n,source3_n,flux1_n,flux2_n,flux3_n,fluxAlt1_n,fluxAlt2_n,fluxAlt3_n,alpha_n,a_n,test_n,g1_Nr,r,dt,t);
  return 1;
}

void read_parameters(char *p_file,int *Nr0,int *set_Nr0,double *epsdis,int *set_epsdis,double *epsiter,int *set_epsiter,double *epsiterid,int *set_epsiterid,int *fout,int *set_fout,char **in_file,int *set_in_file,int *iter,int *set_iter,double *lambda,int *set_lambda,int *level,int *set_level,int *maxstep,int *set_maxstep,int *maxstepid,int *set_maxstepid,char **out_file,int *set_out_file,int **output,int *set_output,double *rmax,int *set_rmax,double *rmin,int *set_rmin,int *s_step,int *set_s_step,int *ser,int *set_ser,double *start_t,int *set_start_t,char **tag,int *set_tag,int **trace,int *set_trace)
{
  if(!*set_Nr0)
    if(get_param(p_file,"Nr0","long",1,Nr0)==1)
      *set_Nr0=1;
  if(!*set_epsdis)
    if(get_param(p_file,"epsdis","double",1,epsdis)==1)
      *set_epsdis=1;
  if(!*set_epsiter)
    if(get_param(p_file,"epsiter","double",1,epsiter)==1)
      *set_epsiter=1;
  if(!*set_epsiterid)
    if(get_param(p_file,"epsiterid","double",1,epsiterid)==1)
      *set_epsiterid=1;
  if(!*set_fout)
    if(get_param(p_file,"fout","long",1,fout)==1)
      *set_fout=1;
  if(!*set_in_file)
    if(get_param(p_file,"in_file","string",1,in_file)==1)
      *set_in_file=1;
  if(!*set_iter)
    if(get_param(p_file,"iter","long",1,iter)==1)
      *set_iter=1;
  if(!*set_lambda)
    if(get_param(p_file,"lambda","double",1,lambda)==1)
      *set_lambda=1;
  if(!*set_level)
    if(get_param(p_file,"level","long",1,level)==1)
      *set_level=1;
  if(!*set_maxstep)
    if(get_param(p_file,"maxstep","long",1,maxstep)==1)
      *set_maxstep=1;
  if(!*set_maxstepid)
    if(get_param(p_file,"maxstepid","long",1,maxstepid)==1)
      *set_maxstepid=1;
  if(!*set_out_file)
    if(get_param(p_file,"out_file","string",1,out_file)==1)
      *set_out_file=1;
  if(!*set_output)
    if(get_param(p_file,"output","ivec",5,output)==1)
      *set_output=1;
  if(!*set_rmax)
    if(get_param(p_file,"rmax","double",1,rmax)==1)
      *set_rmax=1;
  if(!*set_rmin)
    if(get_param(p_file,"rmin","double",1,rmin)==1)
      *set_rmin=1;
  if(!*set_s_step)
    if(get_param(p_file,"s_step","long",1,s_step)==1)
      *set_s_step=1;
  if(!*set_ser)
    if(get_param(p_file,"ser","long",1,ser)==1)
      *set_ser=1;
  if(!*set_start_t)
    if(get_param(p_file,"start_t","double",1,start_t)==1)
      *set_start_t=1;
  if(!*set_tag)
    if(get_param(p_file,"tag","string",1,tag)==1)
      *set_tag=1;
  if(!*set_trace)
    if(get_param(p_file,"trace","ivec",5,trace)==1)
      *set_trace=1;
}

void sread_parameters(char *p_str,int *Nr0,int *set_Nr0,double *epsdis,int *set_epsdis,double *epsiter,int *set_epsiter,double *epsiterid,int *set_epsiterid,int *fout,int *set_fout,char **in_file,int *set_in_file,int *iter,int *set_iter,double *lambda,int *set_lambda,int *level,int *set_level,int *maxstep,int *set_maxstep,int *maxstepid,int *set_maxstepid,char **out_file,int *set_out_file,int **output,int *set_output,double *rmax,int *set_rmax,double *rmin,int *set_rmin,int *s_step,int *set_s_step,int *ser,int *set_ser,double *start_t,int *set_start_t,char **tag,int *set_tag,int **trace,int *set_trace)
{
  if(!*set_Nr0)
    if(sget_param(p_str,"Nr0","long",1,Nr0,1)==1)
      *set_Nr0=1;
  if(!*set_epsdis)
    if(sget_param(p_str,"epsdis","double",1,epsdis,1)==1)
      *set_epsdis=1;
  if(!*set_epsiter)
    if(sget_param(p_str,"epsiter","double",1,epsiter,1)==1)
      *set_epsiter=1;
  if(!*set_epsiterid)
    if(sget_param(p_str,"epsiterid","double",1,epsiterid,1)==1)
      *set_epsiterid=1;
  if(!*set_fout)
    if(sget_param(p_str,"fout","long",1,fout,1)==1)
      *set_fout=1;
  if(!*set_in_file)
    if(sget_param(p_str,"in_file","string",1,in_file,1)==1)
      *set_in_file=1;
  if(!*set_iter)
    if(sget_param(p_str,"iter","long",1,iter,1)==1)
      *set_iter=1;
  if(!*set_lambda)
    if(sget_param(p_str,"lambda","double",1,lambda,1)==1)
      *set_lambda=1;
  if(!*set_level)
    if(sget_param(p_str,"level","long",1,level,1)==1)
      *set_level=1;
  if(!*set_maxstep)
    if(sget_param(p_str,"maxstep","long",1,maxstep,1)==1)
      *set_maxstep=1;
  if(!*set_maxstepid)
    if(sget_param(p_str,"maxstepid","long",1,maxstepid,1)==1)
      *set_maxstepid=1;
  if(!*set_out_file)
    if(sget_param(p_str,"out_file","string",1,out_file,1)==1)
      *set_out_file=1;
  if(!*set_output)
    if(sget_param(p_str,"output","ivec",5,output,1)==1)
      *set_output=1;
  if(!*set_rmax)
    if(sget_param(p_str,"rmax","double",1,rmax,1)==1)
      *set_rmax=1;
  if(!*set_rmin)
    if(sget_param(p_str,"rmin","double",1,rmin,1)==1)
      *set_rmin=1;
  if(!*set_s_step)
    if(sget_param(p_str,"s_step","long",1,s_step,1)==1)
      *set_s_step=1;
  if(!*set_ser)
    if(sget_param(p_str,"ser","long",1,ser,1)==1)
      *set_ser=1;
  if(!*set_start_t)
    if(sget_param(p_str,"start_t","double",1,start_t,1)==1)
      *set_start_t=1;
  if(!*set_tag)
    if(sget_param(p_str,"tag","string",1,tag,1)==1)
      *set_tag=1;
  if(!*set_trace)
    if(sget_param(p_str,"trace","ivec",5,trace,1)==1)
      *set_trace=1;
}

void read_attributes(char *p_file,int *out_gf,int *set_out_gf)
{
  if(!*set_out_gf)
    if(get_param(p_file,"out_gf","long",18,out_gf)==1)
      *set_out_gf=1;
}

int check_params_attribs(int Nr0,int set_Nr0,double epsdis,int set_epsdis,double epsiter,int set_epsiter,double epsiterid,int set_epsiterid,int fout,int set_fout,char *in_file,int set_in_file,int iter,int set_iter,double lambda,int set_lambda,int level,int set_level,int maxstep,int set_maxstep,int maxstepid,int set_maxstepid,char *out_file,int set_out_file,int *output,int set_output,double rmax,int set_rmax,double rmin,int set_rmin,int s_step,int set_s_step,int ser,int set_ser,double start_t,int set_start_t,char *tag,int set_tag,int *trace,int set_trace,int *out_gf,int set_out_gf)
{
  int all_ok=1;
  if(!set_Nr0 && !(1)){
    fprintf(stderr,"ERROR: parameter Nr0 has not been set.\n");
    all_ok=0;
  }else if(!set_Nr0 && 1)
    fprintf(stderr,"WARNING: using default for parameter Nr0.\n");
  if(!set_epsdis && !(1)){
    fprintf(stderr,"ERROR: parameter epsdis has not been set.\n");
    all_ok=0;
  }else if(!set_epsdis && 1)
    fprintf(stderr,"WARNING: using default for parameter epsdis.\n");
  if(!set_epsiter && !(1)){
    fprintf(stderr,"ERROR: parameter epsiter has not been set.\n");
    all_ok=0;
  }else if(!set_epsiter && 1)
    fprintf(stderr,"WARNING: using default for parameter epsiter.\n");
  if(!set_epsiterid && !(1)){
    fprintf(stderr,"ERROR: parameter epsiterid has not been set.\n");
    all_ok=0;
  }else if(!set_epsiterid && 1)
    fprintf(stderr,"WARNING: using default for parameter epsiterid.\n");
  if(!set_fout && !(1)){
    fprintf(stderr,"ERROR: parameter fout has not been set.\n");
    all_ok=0;
  }else if(!set_fout && 1)
    fprintf(stderr,"WARNING: using default for parameter fout.\n");
  if(!set_in_file && !(0)){
    fprintf(stderr,"ERROR: parameter in_file has not been set.\n");
    all_ok=0;
  }else if(!set_in_file && 0)
    fprintf(stderr,"WARNING: using default for parameter in_file.\n");
  if(!set_iter && !(1)){
    fprintf(stderr,"ERROR: parameter iter has not been set.\n");
    all_ok=0;
  }else if(!set_iter && 1)
    fprintf(stderr,"WARNING: using default for parameter iter.\n");
  if(!set_lambda && !(1)){
    fprintf(stderr,"ERROR: parameter lambda has not been set.\n");
    all_ok=0;
  }else if(!set_lambda && 1)
    fprintf(stderr,"WARNING: using default for parameter lambda.\n");
  if(!set_level && !(1)){
    fprintf(stderr,"ERROR: parameter level has not been set.\n");
    all_ok=0;
  }else if(!set_level && 1)
    fprintf(stderr,"WARNING: using default for parameter level.\n");
  if(!set_maxstep && !(1)){
    fprintf(stderr,"ERROR: parameter maxstep has not been set.\n");
    all_ok=0;
  }else if(!set_maxstep && 1)
    fprintf(stderr,"WARNING: using default for parameter maxstep.\n");
  if(!set_maxstepid && !(1)){
    fprintf(stderr,"ERROR: parameter maxstepid has not been set.\n");
    all_ok=0;
  }else if(!set_maxstepid && 1)
    fprintf(stderr,"WARNING: using default for parameter maxstepid.\n");
  if(!set_out_file && !(0)){
    fprintf(stderr,"ERROR: parameter out_file has not been set.\n");
    all_ok=0;
  }else if(!set_out_file && 0)
    fprintf(stderr,"WARNING: using default for parameter out_file.\n");
  if(!set_output && !(1)){
    fprintf(stderr,"ERROR: parameter output has not been set.\n");
    all_ok=0;
  }else if(!set_output && 1)
    fprintf(stderr,"WARNING: using default for parameter output.\n");
  if(!set_rmax && !(1)){
    fprintf(stderr,"ERROR: parameter rmax has not been set.\n");
    all_ok=0;
  }else if(!set_rmax && 1)
    fprintf(stderr,"WARNING: using default for parameter rmax.\n");
  if(!set_rmin && !(1)){
    fprintf(stderr,"ERROR: parameter rmin has not been set.\n");
    all_ok=0;
  }else if(!set_rmin && 1)
    fprintf(stderr,"WARNING: using default for parameter rmin.\n");
  if(!set_s_step && !(1)){
    fprintf(stderr,"ERROR: parameter s_step has not been set.\n");
    all_ok=0;
  }else if(!set_s_step && 1)
    fprintf(stderr,"WARNING: using default for parameter s_step.\n");
  if(!set_ser && !(1)){
    fprintf(stderr,"ERROR: parameter ser has not been set.\n");
    all_ok=0;
  }else if(!set_ser && 1)
    fprintf(stderr,"WARNING: using default for parameter ser.\n");
  if(!set_start_t && !(1)){
    fprintf(stderr,"ERROR: parameter start_t has not been set.\n");
    all_ok=0;
  }else if(!set_start_t && 1)
    fprintf(stderr,"WARNING: using default for parameter start_t.\n");
  if(!set_tag && !(1)){
    fprintf(stderr,"ERROR: parameter tag has not been set.\n");
    all_ok=0;
  }else if(!set_tag && 1)
    fprintf(stderr,"WARNING: using default for parameter tag.\n");
  if(!set_trace && !(1)){
    fprintf(stderr,"ERROR: parameter trace has not been set.\n");
    all_ok=0;
  }else if(!set_trace && 1)
    fprintf(stderr,"WARNING: using default for parameter trace.\n");
  if(!set_out_gf && !(1)){
    fprintf(stderr,"ERROR: attribute out_gf has not been set.\n");
    all_ok=0;
  }else if(!set_out_gf && 1)
    fprintf(stderr,"WARNING: using default for attribute out_gf.\n");
  return(all_ok);
}

void read_state(int argc, char **argv, lattice_type *lats,double *q1_n,double *q2_n,double *q3_n,double *p1_n,double *p2_n,double *p3_n,double *source1_n,double *source2_n,double *source3_n,double *flux1_n,double *flux2_n,double *flux3_n,double *fluxAlt1_n,double *fluxAlt2_n,double *fluxAlt3_n,double *alpha_n,double *a_n,double *test_n,int *Nr0,int set_Nr0,double *epsdis,int set_epsdis,double *epsiter,int set_epsiter,double *epsiterid,int set_epsiterid,int *fout,int set_fout,char **in_file,int set_in_file,int *iter,int set_iter,double *lambda,int set_lambda,int *level,int set_level,int *maxstep,int set_maxstep,int *maxstepid,int set_maxstepid,char **out_file,int set_out_file,int *output,int set_output,double *rmax,int set_rmax,double *rmin,int set_rmin,int *s_step,int set_s_step,int *ser,int set_ser,double *start_t,int set_start_t,char **tag,int set_tag,int *trace,int set_trace)
{
  char command[256];
  FILE *fp;
  int res,i,j;

  fp=fopen(*in_file,"r");
  if(fp==NULL){
    fprintf(stderr,"Can't open %s.\n",*in_file);
    fprintf(stderr,"Calling initial data generator.\n");
    j=sprintf(command,"fluid_init");
    for(i=1;i<argc;i++)
      j+=sprintf(command+j," %s",argv[i]);
    system(command);
    fp=fopen(*in_file,"r");
    if(fp==NULL){
      fprintf(stderr,"Can't open %s.\n",*in_file);
      fprintf(stderr,"Assuming updates will take care of initialization.\n");
      return;
    }else fclose(fp);
  }else fclose(fp);
  gft_set_single(*in_file);
  gft_read_id_gf("q1",lats[0].shape,&lats[0].rank,q1_n);
  gft_read_id_gf("q2",lats[0].shape,&lats[0].rank,q2_n);
  gft_read_id_gf("q3",lats[0].shape,&lats[0].rank,q3_n);
  gft_read_id_gf("p1",lats[0].shape,&lats[0].rank,p1_n);
  gft_read_id_gf("p2",lats[0].shape,&lats[0].rank,p2_n);
  gft_read_id_gf("p3",lats[0].shape,&lats[0].rank,p3_n);
  gft_read_id_gf("source1",lats[0].shape,&lats[0].rank,source1_n);
  gft_read_id_gf("source2",lats[0].shape,&lats[0].rank,source2_n);
  gft_read_id_gf("source3",lats[0].shape,&lats[0].rank,source3_n);
  gft_read_id_gf("flux1",lats[0].shape,&lats[0].rank,flux1_n);
  gft_read_id_gf("flux2",lats[0].shape,&lats[0].rank,flux2_n);
  gft_read_id_gf("flux3",lats[0].shape,&lats[0].rank,flux3_n);
  gft_read_id_gf("fluxAlt1",lats[0].shape,&lats[0].rank,fluxAlt1_n);
  gft_read_id_gf("fluxAlt2",lats[0].shape,&lats[0].rank,fluxAlt2_n);
  gft_read_id_gf("fluxAlt3",lats[0].shape,&lats[0].rank,fluxAlt3_n);
  gft_read_id_gf("alpha",lats[0].shape,&lats[0].rank,alpha_n);
  gft_read_id_gf("a",lats[0].shape,&lats[0].rank,a_n);
  gft_read_id_gf("test",lats[0].shape,&lats[0].rank,test_n);
  if(!set_Nr0)
    gft_read_id_int("Nr0",Nr0,1);
  if(!set_epsdis)
    gft_read_id_float("epsdis",epsdis,1);
  if(!set_epsiter)
    gft_read_id_float("epsiter",epsiter,1);
  if(!set_epsiterid)
    gft_read_id_float("epsiterid",epsiterid,1);
  if(!set_fout)
    gft_read_id_int("fout",fout,1);
  if(!set_in_file)
    gft_read_id_str("in_file",in_file,1);
  if(!set_iter)
    gft_read_id_int("iter",iter,1);
  if(!set_lambda)
    gft_read_id_float("lambda",lambda,1);
  if(!set_level)
    gft_read_id_int("level",level,1);
  if(!set_maxstep)
    gft_read_id_int("maxstep",maxstep,1);
  if(!set_maxstepid)
    gft_read_id_int("maxstepid",maxstepid,1);
  if(!set_out_file)
    gft_read_id_str("out_file",out_file,1);
  if(!set_output)
    gft_read_id_int("output",output,5);
  if(!set_rmax)
    gft_read_id_float("rmax",rmax,1);
  if(!set_rmin)
    gft_read_id_float("rmin",rmin,1);
  if(!set_s_step)
    gft_read_id_int("s_step",s_step,1);
  if(!set_ser)
    gft_read_id_int("ser",ser,1);
  if(!set_start_t)
    gft_read_id_float("start_t",start_t,1);
  if(!set_tag)
    gft_read_id_str("tag",tag,1);
  if(!set_trace)
    gft_read_id_int("trace",trace,5);
  gft_set_multi();
}

void dump_state(const double t, lattice_type *lats,double *q1_n,double *q2_n,double *q3_n,double *p1_n,double *p2_n,double *p3_n,double *source1_n,double *source2_n,double *source3_n,double *flux1_n,double *flux2_n,double *flux3_n,double *fluxAlt1_n,double *fluxAlt2_n,double *fluxAlt3_n,double *alpha_n,double *a_n,double *test_n,int Nr0,double epsdis,double epsiter,double epsiterid,int fout,char *in_file,int iter,double lambda,int level,int maxstep,int maxstepid,char *out_file,int *output,double rmax,double rmin,int s_step,int ser,double start_t,char *tag,int *trace)
{
  gft_set_single(out_file);
  gft_write_id_gf("q1",lats[0].shape,lats[0].rank,q1_n);
  gft_write_id_gf("q2",lats[0].shape,lats[0].rank,q2_n);
  gft_write_id_gf("q3",lats[0].shape,lats[0].rank,q3_n);
  gft_write_id_gf("p1",lats[0].shape,lats[0].rank,p1_n);
  gft_write_id_gf("p2",lats[0].shape,lats[0].rank,p2_n);
  gft_write_id_gf("p3",lats[0].shape,lats[0].rank,p3_n);
  gft_write_id_gf("source1",lats[0].shape,lats[0].rank,source1_n);
  gft_write_id_gf("source2",lats[0].shape,lats[0].rank,source2_n);
  gft_write_id_gf("source3",lats[0].shape,lats[0].rank,source3_n);
  gft_write_id_gf("flux1",lats[0].shape,lats[0].rank,flux1_n);
  gft_write_id_gf("flux2",lats[0].shape,lats[0].rank,flux2_n);
  gft_write_id_gf("flux3",lats[0].shape,lats[0].rank,flux3_n);
  gft_write_id_gf("fluxAlt1",lats[0].shape,lats[0].rank,fluxAlt1_n);
  gft_write_id_gf("fluxAlt2",lats[0].shape,lats[0].rank,fluxAlt2_n);
  gft_write_id_gf("fluxAlt3",lats[0].shape,lats[0].rank,fluxAlt3_n);
  gft_write_id_gf("alpha",lats[0].shape,lats[0].rank,alpha_n);
  gft_write_id_gf("a",lats[0].shape,lats[0].rank,a_n);
  gft_write_id_gf("test",lats[0].shape,lats[0].rank,test_n);
  start_t=t;
  gft_write_id_int("Nr0",&Nr0,1);
  gft_write_id_float("epsdis",&epsdis,1);
  gft_write_id_float("epsiter",&epsiter,1);
  gft_write_id_float("epsiterid",&epsiterid,1);
  gft_write_id_int("fout",&fout,1);
  gft_write_id_str("in_file",&in_file,1);
  gft_write_id_int("iter",&iter,1);
  gft_write_id_float("lambda",&lambda,1);
  gft_write_id_int("level",&level,1);
  gft_write_id_int("maxstep",&maxstep,1);
  gft_write_id_int("maxstepid",&maxstepid,1);
  gft_write_id_str("out_file",&out_file,1);
  gft_write_id_int("output",output,5);
  gft_write_id_float("rmax",&rmax,1);
  gft_write_id_float("rmin",&rmin,1);
  gft_write_id_int("s_step",&s_step,1);
  gft_write_id_int("ser",&ser,1);
  gft_write_id_float("start_t",&start_t,1);
  gft_write_id_str("tag",&tag,1);
  gft_write_id_int("trace",trace,5);
  gft_set_multi();
}

void handler(int sig, int *rnpldone, int *rmod, char **fname, int *out_gf)
{
  char temp[5];
  int quit,ch,i;

  if(sig==SIGQUIT){
    quit=1;
    *rnpldone=1;
  }else quit=0;
  while(!quit){
    fprintf(stdout,"1.  Change output frequency\n");
    fprintf(stdout,"2.  View out_gf\n");
    fprintf(stdout,"3.  Change out_gf\n");
    fprintf(stdout,"4.  Resume\n");
    fprintf(stdout,"5.  Exit\n");
    fprintf(stdout,"Enter choice:");
    scanf("%d",&ch);
    gets(temp);
    switch(ch){
      case 1 : fprintf(stdout,"Enter new frequency: 1/");
               scanf("%d",rmod);
               gets(temp);
               break;
      case 2 : for(i=0;i<18;i++){
                 fprintf(stdout,"%d. %s",i,fname[i]);
                 if(out_gf[i])
                   fprintf(stdout," ON\n");
                 else fprintf(stdout," OFF\n");
               }
               break;
      case 3 : fprintf(stdout,"Enter grid function number: ");
               scanf("%d",&i);
               gets(temp);
               if(i>=0 && i<18){
                 out_gf[i]=!out_gf[i];
                 if(!out_gf[i]) gft_close(fname[i]);
               }
               break;
      case 4 : quit=1;
               break;
      case 5 : quit=1; *rnpldone=1;
               break;
    }
  }
}

void output_gfuncs(int fout, int ser, double t, int *out_gf, char **fname,double *q1_n,double *q2_n,double *q3_n,double *p1_n,double *p2_n,double *p3_n,double *source1_n,double *source2_n,double *source3_n,double *flux1_n,double *flux2_n,double *flux3_n,double *fluxAlt1_n,double *fluxAlt2_n,double *fluxAlt3_n,double *alpha_n,double *a_n,double *test_n,int g1_rank, int *g1_shape, char *g1_cnames, double *g1_crds)
{
  if(out_gf[0]){
    if(fout)
      gft_out_full(fname[0],t,g1_shape,g1_cnames,g1_rank,g1_crds,q1_n);
    if(ser)
      rvsxynt(fname[0],t,g1_crds,q1_n,g1_shape[0]);
  }
  if(out_gf[1]){
    if(fout)
      gft_out_full(fname[1],t,g1_shape,g1_cnames,g1_rank,g1_crds,q2_n);
    if(ser)
      rvsxynt(fname[1],t,g1_crds,q2_n,g1_shape[0]);
  }
  if(out_gf[2]){
    if(fout)
      gft_out_full(fname[2],t,g1_shape,g1_cnames,g1_rank,g1_crds,q3_n);
    if(ser)
      rvsxynt(fname[2],t,g1_crds,q3_n,g1_shape[0]);
  }
  if(out_gf[3]){
    if(fout)
      gft_out_full(fname[3],t,g1_shape,g1_cnames,g1_rank,g1_crds,p1_n);
    if(ser)
      rvsxynt(fname[3],t,g1_crds,p1_n,g1_shape[0]);
  }
  if(out_gf[4]){
    if(fout)
      gft_out_full(fname[4],t,g1_shape,g1_cnames,g1_rank,g1_crds,p2_n);
    if(ser)
      rvsxynt(fname[4],t,g1_crds,p2_n,g1_shape[0]);
  }
  if(out_gf[5]){
    if(fout)
      gft_out_full(fname[5],t,g1_shape,g1_cnames,g1_rank,g1_crds,p3_n);
    if(ser)
      rvsxynt(fname[5],t,g1_crds,p3_n,g1_shape[0]);
  }
  if(out_gf[6]){
    if(fout)
      gft_out_full(fname[6],t,g1_shape,g1_cnames,g1_rank,g1_crds,source1_n);
    if(ser)
      rvsxynt(fname[6],t,g1_crds,source1_n,g1_shape[0]);
  }
  if(out_gf[7]){
    if(fout)
      gft_out_full(fname[7],t,g1_shape,g1_cnames,g1_rank,g1_crds,source2_n);
    if(ser)
      rvsxynt(fname[7],t,g1_crds,source2_n,g1_shape[0]);
  }
  if(out_gf[8]){
    if(fout)
      gft_out_full(fname[8],t,g1_shape,g1_cnames,g1_rank,g1_crds,source3_n);
    if(ser)
      rvsxynt(fname[8],t,g1_crds,source3_n,g1_shape[0]);
  }
  if(out_gf[9]){
    if(fout)
      gft_out_full(fname[9],t,g1_shape,g1_cnames,g1_rank,g1_crds,flux1_n);
    if(ser)
      rvsxynt(fname[9],t,g1_crds,flux1_n,g1_shape[0]);
  }
  if(out_gf[10]){
    if(fout)
      gft_out_full(fname[10],t,g1_shape,g1_cnames,g1_rank,g1_crds,flux2_n);
    if(ser)
      rvsxynt(fname[10],t,g1_crds,flux2_n,g1_shape[0]);
  }
  if(out_gf[11]){
    if(fout)
      gft_out_full(fname[11],t,g1_shape,g1_cnames,g1_rank,g1_crds,flux3_n);
    if(ser)
      rvsxynt(fname[11],t,g1_crds,flux3_n,g1_shape[0]);
  }
  if(out_gf[12]){
    if(fout)
      gft_out_full(fname[12],t,g1_shape,g1_cnames,g1_rank,g1_crds,fluxAlt1_n);
    if(ser)
      rvsxynt(fname[12],t,g1_crds,fluxAlt1_n,g1_shape[0]);
  }
  if(out_gf[13]){
    if(fout)
      gft_out_full(fname[13],t,g1_shape,g1_cnames,g1_rank,g1_crds,fluxAlt2_n);
    if(ser)
      rvsxynt(fname[13],t,g1_crds,fluxAlt2_n,g1_shape[0]);
  }
  if(out_gf[14]){
    if(fout)
      gft_out_full(fname[14],t,g1_shape,g1_cnames,g1_rank,g1_crds,fluxAlt3_n);
    if(ser)
      rvsxynt(fname[14],t,g1_crds,fluxAlt3_n,g1_shape[0]);
  }
  if(out_gf[15]){
    if(fout)
      gft_out_full(fname[15],t,g1_shape,g1_cnames,g1_rank,g1_crds,alpha_n);
    if(ser)
      rvsxynt(fname[15],t,g1_crds,alpha_n,g1_shape[0]);
  }
  if(out_gf[16]){
    if(fout)
      gft_out_full(fname[16],t,g1_shape,g1_cnames,g1_rank,g1_crds,a_n);
    if(ser)
      rvsxynt(fname[16],t,g1_crds,a_n,g1_shape[0]);
  }
  if(out_gf[17]){
    if(fout)
      gft_out_full(fname[17],t,g1_shape,g1_cnames,g1_rank,g1_crds,test_n);
    if(ser)
      rvsxynt(fname[17],t,g1_crds,test_n,g1_shape[0]);
  }
}

int got_signal;
void sig_handler(int sig)
{
  got_signal=sig;
}

#define IVEL 4
#define FVEL 4

int main(int argc, char **argv)
{
  FILE *fp;
  char param_file[60];
  int i,steps,ot;
  double t;
  double *lc;
  time_t st,ed;
  extern char *optarg;
  extern int optind, opterr,optopt;
  int opt,argerr=0;
  int rnpldone=0;
  int rmod;
  double dt;
  double dr;
  double *r;
  int Nr0;
  int set_Nr0;
  double epsdis;
  int set_epsdis;
  double epsiter;
  int set_epsiter;
  double epsiterid;
  int set_epsiterid;
  int fout;
  int set_fout;
  char *in_file;
  int set_in_file;
  int iter;
  int set_iter;
  double lambda;
  int set_lambda;
  int level;
  int set_level;
  int maxstep;
  int set_maxstep;
  int maxstepid;
  int set_maxstepid;
  int set_memsiz;
  char *out_file;
  int set_out_file;
  int *output;
  int set_output;
  double rmax;
  int set_rmax;
  double rmin;
  int set_rmin;
  int s_step;
  int set_s_step;
  int ser;
  int set_ser;
  double start_t;
  int set_start_t;
  char *tag;
  int set_tag;
  int *trace;
  int set_trace;
  int Nr;
  coords g1_0;
  lattice_type *lats;
  char **cname;

  double *q1_np1;
  double *q1_n;
  double *q2_np1;
  double *q2_n;
  double *q3_np1;
  double *q3_n;
  double *p1_np1;
  double *p1_n;
  double *p2_np1;
  double *p2_n;
  double *p3_np1;
  double *p3_n;
  double *source1_n;
  double *source2_n;
  double *source3_n;
  double *flux1_n;
  double *flux2_n;
  double *flux3_n;
  double *fluxAlt1_n;
  double *fluxAlt2_n;
  double *fluxAlt3_n;
  double *alpha_n;
  double *a_n;
  double *test_n;
  char **fname;
  int *out_gf;
  int set_out_gf;


  optarg=NULL; optind=1; opterr=1; optopt=0;
  /* initialize parameters */
  set_Nr0=0;
  Nr0=3;
  set_epsdis=0;
  epsdis=0.5;
  set_epsiter=0;
  epsiter=1e-05;
  set_epsiterid=0;
  epsiterid=1e-05;
  set_fout=0;
  fout=0;
  set_in_file=0;
  in_file=NULL;
  set_iter=0;
  iter=100;
  set_lambda=0;
  lambda=0.5;
  set_level=0;
  level=0;
  set_maxstep=0;
  maxstep=50;
  set_maxstepid=0;
  maxstepid=50;
  set_out_file=0;
  out_file=NULL;
  set_output=0;
  output=ivec_alloc_n(5,"output");
  output[0]=1;
  output[1]=-1;
  output[2]=-1;
  output[3]=1;
  output[4]=0;
  set_rmax=0;
  rmax=10;
  set_rmin=0;
  rmin=0;
  set_s_step=0;
  s_step=0;
  set_ser=0;
  ser=0;
  set_start_t=0;
  start_t=0;
  set_tag=0;
  tag=cvec_alloc(1);
  strcpy(tag,"");
  set_trace=0;
  trace=ivec_alloc_n(5,"trace");
  trace[0]=1;
  trace[1]=-1;
  trace[2]=-1;
  trace[3]=1;
  trace[4]=0;
  set_out_gf=0;
  out_gf=ivec_alloc_n(18,"out_gf");
  out_gf[0]=1;
  out_gf[1]=1;
  out_gf[2]=1;
  out_gf[3]=1;
  out_gf[4]=1;
  out_gf[5]=1;
  out_gf[6]=1;
  out_gf[7]=1;
  out_gf[8]=1;
  out_gf[9]=1;
  out_gf[10]=1;
  out_gf[11]=1;
  out_gf[12]=1;
  out_gf[13]=1;
  out_gf[14]=1;
  out_gf[15]=1;
  out_gf[16]=1;
  out_gf[17]=1;


  /* check command line parameters */
  while((opt=getopt(argc,argv,"p:"))!=EOF){
    switch(opt){
      case 'p' :
        sread_parameters(optarg,&Nr0,&set_Nr0,&epsdis,&set_epsdis,&epsiter,&set_epsiter,&epsiterid,&set_epsiterid,&fout,&set_fout,&in_file,&set_in_file,&iter,&set_iter,&lambda,&set_lambda,&level,&set_level,&maxstep,&set_maxstep,&maxstepid,&set_maxstepid,&out_file,&set_out_file,&output,&set_output,&rmax,&set_rmax,&rmin,&set_rmin,&s_step,&set_s_step,&ser,&set_ser,&start_t,&set_start_t,&tag,&set_tag,&trace,&set_trace);
        break;
      case '?' :
        argerr=1;
        break;
    }
  }
  if(argerr){
    fprintf(stderr," fluid\n");
    fprintf(stderr,"Usage: fluid\n");
    fprintf(stderr,"     [ -p \"parameter:=value\" ]\n");
    fprintf(stderr,"     [ parameter_file ]\n");
    exit(0);
  }
  if(optind<argc){
     strcpy(param_file,argv[optind]);
  }else{
     printf("fluid: Enter name of parameter file: ");
     fflush(stdout);
     scanf("%s",param_file);
  }
  if((fp=fopen(param_file,"r"))==NULL){
    fprintf(stderr,"fluid: Unable to open parameter file %s\n",param_file);
    exit(0);
  }
  fclose(fp);

  read_parameters(param_file,&Nr0,&set_Nr0,&epsdis,&set_epsdis,&epsiter,&set_epsiter,&epsiterid,&set_epsiterid,&fout,&set_fout,&in_file,&set_in_file,&iter,&set_iter,&lambda,&set_lambda,&level,&set_level,&maxstep,&set_maxstep,&maxstepid,&set_maxstepid,&out_file,&set_out_file,&output,&set_output,&rmax,&set_rmax,&rmin,&set_rmin,&s_step,&set_s_step,&ser,&set_ser,&start_t,&set_start_t,&tag,&set_tag,&trace,&set_trace);
  read_attributes(param_file,out_gf,&set_out_gf);
  read_attributes(".rnpl.attributes",out_gf,&set_out_gf);
  /* initialize coordinate differentials */
  if(level > 0 && Nr0%2==1){
    printf("WARNING: Nr0=%d is odd and level != 0\n",Nr0);
    printf("WARNING: Adjusting Nr0 to %d.  If this is part of a convergence\n",Nr0-1);
    printf("WARNING: test, rerun level 0 calculation with Nr0 = %d.\n",Nr0-1);
    Nr0-=1;
  }
  Nr=Nr0*(int)pow(2.0,(double)level)+1;
  dr=(rmax - rmin)/(Nr - 1);
  dt=lambda*sqrt((dr*dr)/1);

  /* initialize grids */
  lats=(lattice_type *)malloc(1*sizeof(lattice_type));
  if(lats==NULL){
    printf("Unable to malloc lats.\n");
    exit(1);
  }
  cname=(char **)malloc(1*sizeof(char *));
  if(cname==NULL){
    printf("Unable to malloc cname.\n");
    exit(1);
  }
  cname[0]=(char *)malloc(2*sizeof(char));
  if(cname[0]==NULL){
    printf("Unable to malloc cname[0].\n");
    exit(1);
  }
  strcpy(cname[0],"r");
  g1_0.n_coord=(Nr-1+1);
  r=vec_alloc_n(g1_0.n_coord,"r");
  rdvramp(r,g1_0.n_coord,rmin,dr);
  g1_0.coord=r;
  lats[0].shape[0]=g1_0.n_coord;
  lats[0].bounds[0]=1;
  lats[0].bounds[1]=Nr;
  lats[0].coords=vec_alloc_n(lats[0].shape[0],"lats[0].coords");
  lc=lats[0].coords;
  rdvcpy(lc,g1_0.coord,lats[0].shape[0]);
  lats[0].rank=1;
  lats[0].cs=0;

  /* initialize grid functions */
  fname=(char **)malloc(18*sizeof(char *));
  q1_np1=vec_alloc_n(Nr,"q1_np1");
  q1_n=vec_alloc_n(Nr,"q1_n");
  fname[0]=(char *)malloc(15*sizeof(char));
  sprintf(fname[0],"%sq1_%d",tag,level);
  q2_np1=vec_alloc_n(Nr,"q2_np1");
  q2_n=vec_alloc_n(Nr,"q2_n");
  fname[1]=(char *)malloc(15*sizeof(char));
  sprintf(fname[1],"%sq2_%d",tag,level);
  q3_np1=vec_alloc_n(Nr,"q3_np1");
  q3_n=vec_alloc_n(Nr,"q3_n");
  fname[2]=(char *)malloc(15*sizeof(char));
  sprintf(fname[2],"%sq3_%d",tag,level);
  p1_np1=vec_alloc_n(Nr,"p1_np1");
  p1_n=vec_alloc_n(Nr,"p1_n");
  fname[3]=(char *)malloc(15*sizeof(char));
  sprintf(fname[3],"%sp1_%d",tag,level);
  p2_np1=vec_alloc_n(Nr,"p2_np1");
  p2_n=vec_alloc_n(Nr,"p2_n");
  fname[4]=(char *)malloc(15*sizeof(char));
  sprintf(fname[4],"%sp2_%d",tag,level);
  p3_np1=vec_alloc_n(Nr,"p3_np1");
  p3_n=vec_alloc_n(Nr,"p3_n");
  fname[5]=(char *)malloc(15*sizeof(char));
  sprintf(fname[5],"%sp3_%d",tag,level);
  source1_n=vec_alloc_n(Nr,"source1_n");
  fname[6]=(char *)malloc(20*sizeof(char));
  sprintf(fname[6],"%ssource1_%d",tag,level);
  source2_n=vec_alloc_n(Nr,"source2_n");
  fname[7]=(char *)malloc(20*sizeof(char));
  sprintf(fname[7],"%ssource2_%d",tag,level);
  source3_n=vec_alloc_n(Nr,"source3_n");
  fname[8]=(char *)malloc(20*sizeof(char));
  sprintf(fname[8],"%ssource3_%d",tag,level);
  flux1_n=vec_alloc_n(Nr,"flux1_n");
  fname[9]=(char *)malloc(18*sizeof(char));
  sprintf(fname[9],"%sflux1_%d",tag,level);
  flux2_n=vec_alloc_n(Nr,"flux2_n");
  fname[10]=(char *)malloc(18*sizeof(char));
  sprintf(fname[10],"%sflux2_%d",tag,level);
  flux3_n=vec_alloc_n(Nr,"flux3_n");
  fname[11]=(char *)malloc(18*sizeof(char));
  sprintf(fname[11],"%sflux3_%d",tag,level);
  fluxAlt1_n=vec_alloc_n(Nr,"fluxAlt1_n");
  fname[12]=(char *)malloc(21*sizeof(char));
  sprintf(fname[12],"%sfluxAlt1_%d",tag,level);
  fluxAlt2_n=vec_alloc_n(Nr,"fluxAlt2_n");
  fname[13]=(char *)malloc(21*sizeof(char));
  sprintf(fname[13],"%sfluxAlt2_%d",tag,level);
  fluxAlt3_n=vec_alloc_n(Nr,"fluxAlt3_n");
  fname[14]=(char *)malloc(21*sizeof(char));
  sprintf(fname[14],"%sfluxAlt3_%d",tag,level);
  alpha_n=vec_alloc_n(Nr,"alpha_n");
  fname[15]=(char *)malloc(18*sizeof(char));
  sprintf(fname[15],"%salpha_%d",tag,level);
  a_n=vec_alloc_n(Nr,"a_n");
  fname[16]=(char *)malloc(14*sizeof(char));
  sprintf(fname[16],"%sa_%d",tag,level);
  test_n=vec_alloc_n(Nr,"test_n");
  fname[17]=(char *)malloc(17*sizeof(char));
  sprintf(fname[17],"%stest_%d",tag,level);

  read_state(argc, argv, lats,q1_n,q2_n,q3_n,p1_n,p2_n,p3_n,source1_n,source2_n,source3_n,flux1_n,flux2_n,flux3_n,fluxAlt1_n,fluxAlt2_n,fluxAlt3_n,alpha_n,a_n,test_n,&Nr0,set_Nr0,&epsdis,set_epsdis,&epsiter,set_epsiter,&epsiterid,set_epsiterid,&fout,set_fout,&in_file,set_in_file,&iter,set_iter,&lambda,set_lambda,&level,set_level,&maxstep,set_maxstep,&maxstepid,set_maxstepid,&out_file,set_out_file,output,set_output,&rmax,set_rmax,&rmin,set_rmin,&s_step,set_s_step,&ser,set_ser,&start_t,set_start_t,&tag,set_tag,trace,set_trace);
  if(!check_params_attribs(Nr0,set_Nr0,epsdis,set_epsdis,epsiter,set_epsiter,epsiterid,set_epsiterid,fout,set_fout,in_file,set_in_file,iter,set_iter,lambda,set_lambda,level,set_level,maxstep,set_maxstep,maxstepid,set_maxstepid,out_file,set_out_file,output,set_output,rmax,set_rmax,rmin,set_rmin,s_step,set_s_step,ser,set_ser,start_t,set_start_t,tag,set_tag,trace,set_trace,out_gf,set_out_gf)){
    fprintf(stderr,"fluid: Unable to continue.\n");
    exit(0);
  }

  iter*=(int)pow(2.0,(double)level);

  /* set up signal handler */
  got_signal=FALSE;
  signal(SIGQUIT,sig_handler);
  signal(SIGINT,sig_handler);
  st=time(&st);
  rmod=-1;
  t=start_t;
  i=s_step;
  rnpldone=0;
  /* fix a couple of defaults */
  if(!set_epsiterid)
    epsiterid=epsiter;
  if(!set_trace){
    if(trace[0]<output[0]){
      free(trace);
      trace=ivec_alloc(output[0]*IVEL+1);
    }
    rivcpy(trace,output,output[0]*IVEL+1);
  }
  /* fix ivec indicies */
  fixup_ivec(0,iter,level,output);
  fixup_ivec(0,iter,level,trace);

  if(i==0){
    if(do_ivec(i,iter,trace))
      fprintf(stdout,"Starting evolution. step: %d at t=%g\n",i,t);
    if(((rmod!=-1) && (i%rmod==0)) || do_ivec(i,iter,output)){
      output_gfuncs(fout,ser,t,out_gf,fname,q1_n,q2_n,q3_n,p1_n,p2_n,p3_n,source1_n,source2_n,source3_n,flux1_n,flux2_n,flux3_n,fluxAlt1_n,fluxAlt2_n,fluxAlt3_n,alpha_n,a_n,test_n,lats[0].rank,lats[0].shape,cname[lats[0].cs],lats[0].coords);
    }
  }
  for(i++;i<=iter && !rnpldone;i++){
    steps=update_gfuncs(&rnpldone,maxstep,epsiter,q1_np1,q1_n,q2_np1,q2_n,q3_np1,q3_n,p1_np1,p1_n,p2_np1,p2_n,p3_np1,p3_n,source1_n,source2_n,source3_n,flux1_n,flux2_n,flux3_n,fluxAlt1_n,fluxAlt2_n,fluxAlt3_n,alpha_n,a_n,test_n,lats[0].shape[0],r,dt,t);
    swap_levels(&q1_np1,&q1_n,&q2_np1,&q2_n,&q3_np1,&q3_n,&p1_np1,&p1_n,&p2_np1,&p2_n,&p3_np1,&p3_n);
    t+=dt;
    if(do_ivec(i,iter,trace))
      printf("step: %d  t=%g  steps=%d\n",i,t,steps);
    if(((rmod!=-1) && (i%rmod==0)) || do_ivec(i,iter,output)){
      output_gfuncs(fout,ser,t,out_gf,fname,q1_n,q2_n,q3_n,p1_n,p2_n,p3_n,source1_n,source2_n,source3_n,flux1_n,flux2_n,flux3_n,fluxAlt1_n,fluxAlt2_n,fluxAlt3_n,alpha_n,a_n,test_n,lats[0].rank,lats[0].shape,cname[lats[0].cs],lats[0].coords);
    }
    if(got_signal){
      handler(got_signal,&rnpldone,&rmod,fname,out_gf);
      got_signal=FALSE;
      signal(SIGQUIT,sig_handler);
      signal(SIGINT,sig_handler);
    }
  }
  ed=time(&ed);
  fprintf(stdout,"Elapsed: %lu sec\n",ed-st);
  s_step=i-1;
  dump_state(t, lats,q1_n,q2_n,q3_n,p1_n,p2_n,p3_n,source1_n,source2_n,source3_n,flux1_n,flux2_n,flux3_n,fluxAlt1_n,fluxAlt2_n,fluxAlt3_n,alpha_n,a_n,test_n,Nr0,epsdis,epsiter,epsiterid,fout,in_file,iter,lambda,level,maxstep,maxstepid,out_file,output,rmax,rmin,s_step,ser,start_t,tag,trace);
  free(q1_np1);
  free(q1_n);
  free(q2_np1);
  free(q2_n);
  free(q3_np1);
  free(q3_n);
  free(p1_np1);
  free(p1_n);
  free(p2_np1);
  free(p2_n);
  free(p3_np1);
  free(p3_n);
  free(source1_n);
  free(source2_n);
  free(source3_n);
  free(flux1_n);
  free(flux2_n);
  free(flux3_n);
  free(fluxAlt1_n);
  free(fluxAlt2_n);
  free(fluxAlt3_n);
  free(alpha_n);
  free(a_n);
  free(test_n);
  gft_close_all();
}

