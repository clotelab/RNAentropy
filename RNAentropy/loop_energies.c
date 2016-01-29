#include "loop_energies.h"
/*
#################################
# BEGIN OF FUNCTION DEFINITIONS #
#################################
*/
 double E_Hairpin(int size, int type, int si1, int sj1, const char *string, paramT *P){
  double energy;

  energy = (size <= 30) ? P->hairpin[size] : P->hairpin[30]+(P->lxc*log((size)/30.));
  if (P->model_details.special_hp){
    if (size == 4) { /* check for tetraloop bonus */
      char tl[7]={0}, *ts;
      strncpy(tl, string, 6);
      if ((ts=strstr(P->Tetraloops, tl)))
        return (P->Tetraloop_E[(ts - P->Tetraloops)/7]);
    }
    else if (size == 6) {
      char tl[9]={0}, *ts;
      strncpy(tl, string, 8);
      if ((ts=strstr(P->Hexaloops, tl)))
        return (energy = P->Hexaloop_E[(ts - P->Hexaloops)/9]);
    }
    else if (size == 3) {
      char tl[6]={0,0,0,0,0,0}, *ts;
      strncpy(tl, string, 5);
      if ((ts=strstr(P->Triloops, tl))) {
        return (P->Triloop_E[(ts - P->Triloops)/6]);
      }
      energy = (energy + (type>2 ? P->TerminalAU : 0));
      return energy;
    }
  }
  energy += P->mismatchH[type][si1][sj1];
  return energy;
}

 double E_IntLoop(int n1, int n2, int type, int type_2, int si1, int sj1, int sp1, int sq1, paramT *P){
  /* compute energy of degree 2 loop (stack bulge or interior) */
  int nl, ns;
  double energy;
  energy = INF;

  if (n1>n2) { nl=n1; ns=n2;}
  else {nl=n2; ns=n1;}

  if (nl == 0)
    return P->stack[type][type_2];  /* stack */

  if (ns==0) {                      /* bulge */
    energy = (nl<=MAXLOOP)?P->bulge[nl]:
      (P->bulge[30]+(int)(P->lxc*log(nl/30.)));
    if (nl==1) energy += P->stack[type][type_2];
    else {
      if (type>2) energy += P->TerminalAU;
      if (type_2>2) energy += P->TerminalAU;
    }
    return energy;
  }
  else {                            /* interior loop */
    if (ns==1) {
      if (nl==1)                    /* 1x1 loop */
        return P->int11[type][type_2][si1][sj1];
      if (nl==2) {                  /* 2x1 loop */
        if (n1==1)
          energy = P->int21[type][type_2][si1][sq1][sj1];
        else
          energy = P->int21[type_2][type][sq1][si1][sp1];
        return energy;
      }
      else {  /* 1xn loop */
        energy = (nl+1<=MAXLOOP)?(P->internal_loop[nl+1]) : (P->internal_loop[30]+(int)(P->lxc*log((nl+1)/30.)));
        energy += MIN2(MAX_NINIO, (nl-ns)*P->ninio[2]);
        energy += P->mismatch1nI[type][si1][sj1] + P->mismatch1nI[type_2][sq1][sp1];
        return energy;
      }
    }
    else if (ns==2) {
      if(nl==2)      {              /* 2x2 loop */
        return P->int22[type][type_2][si1][sp1][sq1][sj1];}
      else if (nl==3){              /* 2x3 loop */
        energy = P->internal_loop[5]+P->ninio[2];
        energy += P->mismatch23I[type][si1][sj1] + P->mismatch23I[type_2][sq1][sp1];
        return energy;
      }

    }
    { /* generic interior loop (no else here!)*/
      energy = (n1+n2<=MAXLOOP)?(P->internal_loop[n1+n2]) : (P->internal_loop[30]+(int)(P->lxc*log((n1+n2)/30.)));

      energy += MIN2(MAX_NINIO, (nl-ns)*P->ninio[2]);

      energy += P->mismatchI[type][si1][sj1] + P->mismatchI[type_2][sq1][sp1];
    }
  }
  return energy;
}

 double E_Stem(int type, int si1, int sj1, int extLoop, paramT *P){
  double energy = 0;
  double d5 = (si1 >= 0) ? P->dangle5[type][si1] : 0;
  double d3 = (sj1 >= 0) ? P->dangle3[type][sj1] : 0;

  if(type > 2)
    energy += P->TerminalAU;

  if(si1 >= 0 && sj1 >= 0)
    energy += (extLoop) ? P->mismatchExt[type][si1][sj1] : P->mismatchM[type][si1][sj1];
  else
    energy += d5 + d3;

  if(!extLoop) energy += P->MLintern[type];
  return energy;
}

 double E_ExtLoop(int type, int si1, int sj1, paramT *P){
  double energy = 0;
  if(si1 >= 0 && sj1 >= 0){
    energy += P->mismatchExt[type][si1][sj1];
  }
  else if (si1 >= 0){
    energy += P->dangle5[type][si1];
  }
  else if (sj1 >= 0){
    energy += P->dangle3[type][sj1];
  }

  if(type > 2)
    energy += P->TerminalAU;

  return energy;
}

 double E_MLstem(int type, int si1, int sj1, paramT *P){
  double energy = 0;
  if(si1 >= 0 && sj1 >= 0){
    energy += P->mismatchM[type][si1][sj1];
  }
  else if (si1 >= 0){
    energy += P->dangle5[type][si1];
  }
  else if (sj1 >= 0){
    energy += P->dangle3[type][sj1];
  }

  if(type > 2)
    energy += P->TerminalAU;

  energy += P->MLintern[type];

  return energy;
}

 double exp_E_Hairpin(int u, int type, short si1, short sj1, const char *string, pf_paramT *P){
  double q, kT;
  kT = P->kT;   /* kT in cal/mol  */

  if(u <= 30)
    q = P->exphairpin[u];
  else
    q = P->exphairpin[30] * exp( -(P->lxc*log( u/30.))*10./kT);

  if(u < 3) return q; /* should only be the case when folding alignments */

  if(P->model_details.special_hp){
    if(u==4) {
      char tl[7]={0,0,0,0,0,0,0}, *ts;
      strncpy(tl, string, 6);
      if ((ts=strstr(P->Tetraloops, tl))){
        if(type != 7)
          return (P->exptetra[(ts-P->Tetraloops)/7]);
        else
          q *= P->exptetra[(ts-P->Tetraloops)/7];
      }
    }
    if (u==6) {
      char tl[9]={0,0,0,0,0,0,0,0,0}, *ts;
      strncpy(tl, string, 8);
      if ((ts=strstr(P->Hexaloops, tl)))
        return  (P->exphex[(ts-P->Hexaloops)/9]);
    }
    if (u==3) {
      char tl[6]={0,0,0,0,0,0}, *ts;
      strncpy(tl, string, 5);
      if ((ts=strstr(P->Triloops, tl)))
        return (P->exptri[(ts-P->Triloops)/6]);
      if (type>2)
        return q * P->expTermAU;
      return q;
    }
  }
  /* no mismatches for tri-loops */
  q *= P->expmismatchH[type][si1][sj1];

  return q;
}

 double exp_E_IntLoop(int u1, int u2, int type, int type2, short si1, short sj1, short sp1, short sq1, pf_paramT *P){
  int ul, us, no_close = 0;
  double z = 0.;

  if ((no_closingGU) && ((type2==3)||(type2==4)||(type==3)||(type==4)))
    no_close = 1;

  if (u1>u2) { ul=u1; us=u2;}
  else {ul=u2; us=u1;}

  if (ul==0) /* stack */
    z = P->expstack[type][type2];
  else if(!no_close){
    if (us==0) {                      /* bulge */
      z = P->expbulge[ul];
      if (ul==1) z *= P->expstack[type][type2];
      else {
        if (type>2) z *= P->expTermAU;
        if (type2>2) z *= P->expTermAU;
      }
      return z;
    }
    else if (us==1) {
      if (ul==1){                    /* 1x1 loop */
        return P->expint11[type][type2][si1][sj1];
      }
      if (ul==2) {                  /* 2x1 loop */
        if (u1==1)
          return P->expint21[type][type2][si1][sq1][sj1];
        else
          return P->expint21[type2][type][sq1][si1][sp1];
      }
      else {  /* 1xn loop */
        z = P->expinternal[ul+us] * P->expmismatch1nI[type][si1][sj1] * P->expmismatch1nI[type2][sq1][sp1];
        return z * P->expninio[2][ul-us];
      }
    }
    else if (us==2) {
      if(ul==2) /* 2x2 loop */
        return P->expint22[type][type2][si1][sp1][sq1][sj1];
      else if(ul==3){              /* 2x3 loop */
        z = P->expinternal[5]*P->expmismatch23I[type][si1][sj1]*P->expmismatch23I[type2][sq1][sp1];
        return z * P->expninio[2][1];
      }
    }
    /* generic interior loop (no else here!)*/
    z = P->expinternal[ul+us] * P->expmismatchI[type][si1][sj1] * P->expmismatchI[type2][sq1][sp1];
    return z * P->expninio[2][ul-us];

  }
  return z;
}

 double exp_E_Stem(int type, int si1, int sj1, int extLoop, pf_paramT *P){
  double energy = 1.0;
  double d5 = (si1 >= 0) ? P->expdangle5[type][si1] : 1.;
  double d3 = (sj1 >= 0) ? P->expdangle3[type][sj1] : 1.;

  if(type > 2)
    energy *= P->expTermAU;

  if(si1 >= 0 && sj1 >= 0)
    energy *= (extLoop) ? P->expmismatchExt[type][si1][sj1] : P->expmismatchM[type][si1][sj1];
  else
    energy *= d5 * d3;

  if(!extLoop) energy *= P->expMLintern[type];
  return energy;
}

 double exp_E_MLstem(int type, int si1, int sj1, pf_paramT *P){
  double energy = 1.0;
  if(si1 >= 0 && sj1 >= 0){
    energy *= P->expmismatchM[type][si1][sj1];
  }
  else if(si1 >= 0){
    energy *= P->expdangle5[type][si1];
  }
  else if(sj1 >= 0){
    energy *= P->expdangle3[type][sj1];
  }

  if(type > 2)
    energy *= P->expTermAU;

  energy *= P->expMLintern[type];
  return energy;
}

 double exp_E_ExtLoop(int type, int si1, int sj1, pf_paramT *P){
  double energy = 1.0;
  if(si1 >= 0 && sj1 >= 0){
    energy *= P->expmismatchExt[type][si1][sj1];
  }
  else if(si1 >= 0){
    energy *= P->expdangle5[type][si1];
  }
  else if(sj1 >= 0){
    energy *= P->expdangle3[type][sj1];
  }

  if(type > 2)
    energy *= P->expTermAU;

  return energy;
}

 double     E_IntLoop_Co(int type, int type_2, int i, int j, int p, int q, int cutpoint, short si1, short sj1, short sp1, short sq1, int dangles, paramT *P){
  double energy = 0;
  if(type > 2)   energy += P->TerminalAU;
  if(type_2 > 2) energy += P->TerminalAU;

  if(!dangles) return energy;

  int ci = (i>=cutpoint)||((i+1)<cutpoint);
  int cj = ((j-1)>=cutpoint)||(j<cutpoint);
  int cp = ((p-1)>=cutpoint)||(p<cutpoint);
  int cq = (q>=cutpoint)||((q+1)<cutpoint);

  double d3    = ci  ? P->dangle3[type][si1]   : 0;
  double d5    = cj  ? P->dangle5[type][sj1]   : 0;
  double d5_2  = cp  ? P->dangle5[type_2][sp1] : 0;
  double d3_2  = cq  ? P->dangle3[type_2][sq1] : 0;

  double tmm   = (cj && ci) ? P->mismatchExt[type][sj1][si1]   : d5 + d3;
  double tmm_2 = (cp && cq) ? P->mismatchExt[type_2][sp1][sq1] : d5_2 + d3_2;

  if(dangles == 2) return energy + tmm + tmm_2;

  /* now we may have non-double dangles only */
  if(i+2 < p){
    if(q+2 < j){ energy += tmm + tmm_2;}
    else if(q+2 == j){ energy += (cj && cq) ? MIN2(tmm + d5_2, tmm_2 + d3) : tmm + tmm_2;}
    else energy += d3 + d5_2;
  }
  else if(i+2 == p){
    if(q+2 < j){ energy += (ci && cp) ? MIN2(tmm + d3_2, tmm_2 + d5) : tmm + tmm_2;}
    else if(q+2 == j){
      energy += MIN2(tmm, MIN2(tmm_2, MIN2(d5 + d5_2, d3 + d3_2)));
    }
    else energy += MIN2(d3, d5_2);
  }
  else{
    if(q+2 < j){ energy += d5 + d3_2;}
    else if(q+2 == j){ energy += MIN2(d5, d3_2);}
  }
  return energy;
}
