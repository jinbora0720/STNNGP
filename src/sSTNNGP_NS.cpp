#define USE_FC_LEN_T
#include <string>
#include <iostream>
#include "util.h"

#ifdef _OPENMP
#include <omp.h>
#endif 

#include <R.h>
#include <Rmath.h>
#include <Rinternals.h>
#include <R_ext/Linpack.h>
#include <R_ext/Lapack.h>
#include <R_ext/BLAS.h>
#ifndef FCONE
# define FCONE
#endif

void updateBF_parent_NS(double *B, double *F, double *c, double *C, double *coords,
                        int *nnIndx, int *nnIndxLU,
                        int *nnIndxParent, int *nnIndxLUParent, int *nnIndxLUAll,
                        int n, int m, int nIndx, int twomp1, int twomp1sq, int q,
                        double *sigmaSq, double *phi, double *nu,
                        double *rho, double *crossphi, int *adjvec, 
                        int covModel, double nThreads, double *nuUnifb) {
  
  int i, k, k2, l, l2, MeIndx, ParentIndx; // change: MeIndx added
  int info = 0;
  int inc = 1;
  double one = 1.0;
  double zero = 0.0;
  char lower = 'L';
  
  double nuMax = 0;
  int nb = 0;
  int threadID = 0;
  double e, crosssigma, crossnu;
  
  for (MeIndx = 0; MeIndx < q; MeIndx++) {
    if(nuUnifb[MeIndx] > nuMax){
      nuMax = nuUnifb[MeIndx];
    }
  }
  nb = 1+static_cast<int>(floor(nuMax));
  //bk must be 1+(int)floor(alpha) * nthread
  double *bk = (double *) R_alloc(nThreads*nb, sizeof(double));
  
  for (MeIndx = 0; MeIndx < q; MeIndx++) { // change: for loop for MeIndx
    if (MeIndx == 0) { // change: root
#ifdef _OPENMP
#pragma omp parallel for private(k, l, info, threadID, e)
#endif
      for(i = 0; i < n; i++){
#ifdef _OPENMP
        threadID = omp_get_thread_num();
#endif
        if (i > 0) {
          for(k = 0; k < nnIndxLU[n+i]; k++){
            e = dist2(coords[i], coords[n+i], coords[nnIndx[nnIndxLU[i]+k]], coords[n+nnIndx[nnIndxLU[i]+k]]);
            c[twomp1*threadID+k] = sigmaSq[MeIndx]*spCor(e, phi[MeIndx], nu[MeIndx], covModel, &bk[threadID*nb]);
            for(l = 0; l <= k; l++){
              e = dist2(coords[nnIndx[nnIndxLU[i]+k]], coords[n+nnIndx[nnIndxLU[i]+k]], coords[nnIndx[nnIndxLU[i]+l]], coords[n+nnIndx[nnIndxLU[i]+l]]);  
              C[twomp1sq*threadID+l*nnIndxLU[n+i]+k] = sigmaSq[MeIndx]*spCor(e, phi[MeIndx], nu[MeIndx], covModel, &bk[threadID*nb]); 
            }
          }
          F77_NAME(dpotrf)(&lower, &nnIndxLU[n+i], &C[twomp1sq*threadID], &nnIndxLU[n+i], &info FCONE);
          if(info != 0){error("c++ error: dpotrf failed\n");}
          F77_NAME(dpotri)(&lower, &nnIndxLU[n+i], &C[twomp1sq*threadID], &nnIndxLU[n+i], &info FCONE);
          F77_NAME(dsymv)(&lower, &nnIndxLU[n+i], &one, &C[twomp1sq*threadID], &nnIndxLU[n+i], &c[twomp1*threadID], &inc, &zero, &B[nIndx*MeIndx+nnIndxLU[i]], &inc FCONE); // change: B index
          F[n*MeIndx+i] = sigmaSq[MeIndx] - F77_NAME(ddot)(&nnIndxLU[n+i], &B[nIndx*MeIndx+nnIndxLU[i]], &inc, &c[twomp1*threadID], &inc); // change: B and F index
        } else {
          B[nIndx*MeIndx+i] = 0; // change: B index
          F[n*MeIndx+i] = sigmaSq[MeIndx]; // change: F index
        }
      }
    } else {
      ParentIndx = which(1, &adjvec[q*MeIndx], q);
      crosssigma = pow(sigmaSq[MeIndx],0.5)*pow(sigmaSq[ParentIndx],0.5);
      crossnu = 0.5*nu[MeIndx] + 0.5*nu[ParentIndx];
      
#ifdef _OPENMP
#pragma omp parallel for private(k, l, info, threadID, e)
#endif
      for(i = 0; i < n; i++){
#ifdef _OPENMP
        threadID = omp_get_thread_num();
#endif
        for(k = 0; k < nnIndxLUAll[n+i]; k++){
          if (k < m+1) {
            e = dist2(coords[i], coords[n+i],
                      coords[nnIndxParent[nnIndxLUParent[i]+k]],coords[n+nnIndxParent[nnIndxLUParent[i]+k]]);
            c[twomp1*threadID+k] = rho[MeIndx-1]*crosssigma*
              spCor(e, crossphi[MeIndx-1], crossnu, covModel, &bk[threadID*nb]); // change: rho[MeIndx-1]
          } else {
            k2 = k-m-1;
            e = dist2(coords[i], coords[n+i],
                      coords[nnIndx[nnIndxLU[i]+k2]], coords[n+nnIndx[nnIndxLU[i]+k2]]);
            c[twomp1*threadID+k] = sigmaSq[MeIndx]*
              spCor(e, phi[MeIndx], nu[MeIndx], covModel, &bk[threadID*nb]);
          }
          
          for(l = 0; l <= k; l++){
            if (k < m+1) {
              e = dist2(coords[nnIndxParent[nnIndxLUParent[i]+k]],coords[n+nnIndxParent[nnIndxLUParent[i]+k]],
                        coords[nnIndxParent[nnIndxLUParent[i]+l]],coords[n+nnIndxParent[nnIndxLUParent[i]+l]]);
              C[twomp1sq*threadID+l*nnIndxLUAll[n+i]+k] = sigmaSq[ParentIndx]*
                spCor(e, phi[ParentIndx], nu[ParentIndx], covModel, &bk[threadID*nb]);
            } else {
              k2 = k-m-1;
              if (l < m+1) {
                e = dist2(coords[nnIndx[nnIndxLU[i]+k2]],coords[n+nnIndx[nnIndxLU[i]+k2]],
                          coords[nnIndxParent[nnIndxLUParent[i]+l]],coords[n+nnIndxParent[nnIndxLUParent[i]+l]]);
                C[twomp1sq*threadID+l*nnIndxLUAll[n+i]+k] = rho[MeIndx-1]*crosssigma*
                  spCor(e, crossphi[MeIndx-1], crossnu, covModel, &bk[threadID*nb]); // change: rho[MeIndx-1]
              } else {
                l2 = l-m-1;
                e = dist2(coords[nnIndx[nnIndxLU[i]+k2]],coords[n+nnIndx[nnIndxLU[i]+k2]],
                          coords[nnIndx[nnIndxLU[i]+l2]],coords[n+nnIndx[nnIndxLU[i]+l2]]);
                C[twomp1sq*threadID+l*nnIndxLUAll[n+i]+k] = sigmaSq[MeIndx]*
                  spCor(e, phi[MeIndx], nu[MeIndx], covModel, &bk[threadID*nb]);
              }
            }
          }
        }
        F77_NAME(dpotrf)(&lower, &nnIndxLUAll[n+i], &C[twomp1sq*threadID], &nnIndxLUAll[n+i], &info FCONE);
        if(info != 0){error("c++ error: dpotrf failed\n");}
        F77_NAME(dpotri)(&lower, &nnIndxLUAll[n+i], &C[twomp1sq*threadID], &nnIndxLUAll[n+i], &info FCONE);
        F77_NAME(dsymv)(&lower, &nnIndxLUAll[n+i], &one, &C[twomp1sq*threadID], &nnIndxLUAll[n+i], &c[twomp1*threadID], &inc, &zero, &B[nIndx*MeIndx+nnIndxLUAll[i]], &inc FCONE); // change: B index
        F[n*MeIndx+i] = sigmaSq[MeIndx] - F77_NAME(ddot)(&nnIndxLUAll[n+i], &B[nIndx*MeIndx+nnIndxLUAll[i]], &inc, &c[twomp1*threadID], &inc); // change: B and F index                     
      }
    }
  }
}

extern "C" {
  
  SEXP sSTNNGP_NS(SEXP y_r, SEXP X_r, 
                  SEXP q_r,                                                     // BJ
                  SEXP p_r, SEXP n_r, SEXP m_r, SEXP coords_r, 
                  SEXP covModel_r, SEXP nnIndx_r, SEXP nnIndxLU_r, 
                  SEXP nnIndxParent_r, SEXP nnIndxLUParent_r, SEXP nnIndxLUAll_r, // BJ
                  SEXP uIndx_r, SEXP uIndxLU_r, SEXP uiIndx_r,     
                  SEXP uuIndx_r, SEXP uuIndxLU_r, SEXP uuiIndx_r,               // BJ
                  SEXP cIndx_r, SEXP cIndxLU_r,                                 // BJ
                  SEXP sigmaSqIGa_r, SEXP sigmaSqIGb_r,                         // BJ
                  SEXP tauSqIGa_r, SEXP tauSqIGb_r,                             // BJ
                  SEXP phiUnifa_r, SEXP phiUnifb_r,                             // BJ
                  SEXP nuUnifa_r, SEXP nuUnifb_r,                               // BJ
                  SEXP betaStarting_r, SEXP sigmaSqStarting_r, SEXP tauSqStarting_r, 
                  SEXP phiStarting_r, 
                  SEXP crossphiStarting_r,                                      // BJ
                  SEXP nuStarting_r, 
                  SEXP rhoStarting_r, SEXP adjmatStarting_r,                    // BJ
                  SEXP sigmaSqTuning_r,                                         // BJ
                  SEXP phiTuning_r, 
                  SEXP crossphiTuning_r,                                        // BJ
                  SEXP nuTuning_r, 
                  SEXP rhoTuning_r,                                             // BJ
                  SEXP nSamples_r, SEXP nThreads_r, SEXP verbose_r, SEXP nReport_r){ 
    
    int h, i, j, j2, k, k2, l, o, s, info, nProtect=0;
    const int inc = 1;
    const double one = 1.0;
    const double negOne = -1.0;
    const double zero = 0.0;
    char const *lower = "L";
    char const *ntran = "N";
    char const *ytran = "T";
    
    //get args
    double *y = REAL(y_r);
    double *X = REAL(X_r);
    int q = INTEGER(q_r)[0];                                                    // BJ
    int p = INTEGER(p_r)[0];
    int n = INTEGER(n_r)[0];
    int m = INTEGER(m_r)[0];
    double *coords = REAL(coords_r);
    int *nnIndx = INTEGER(nnIndx_r);
    int *nnIndxLU = INTEGER(nnIndxLU_r);
    int *nnIndxParent = INTEGER(nnIndxParent_r);                                // BJ
    int *nnIndxLUParent = INTEGER(nnIndxLUParent_r);                            // BJ
    int *nnIndxLUAll = INTEGER(nnIndxLUAll_r);                                  // BJ
    int *uIndx = INTEGER(uIndx_r);
    int *uIndxLU = INTEGER(uIndxLU_r);
    int *uiIndx = INTEGER(uiIndx_r);
    int *uuIndx = INTEGER(uuIndx_r);                                            // BJ
    int *uuIndxLU = INTEGER(uuIndxLU_r);                                        // BJ
    int *uuiIndx = INTEGER(uuiIndx_r);                                          // BJ
    int *cIndx = INTEGER(cIndx_r);                                              // BJ
    int *cIndxLU = INTEGER(cIndxLU_r);                                          // BJ
    int covModel = INTEGER(covModel_r)[0];
    std::string corName = getCorName(covModel);
    
    //priors
    double *sigmaSqIGa = REAL(sigmaSqIGa_r);                                    // BJ
    double *sigmaSqIGb = REAL(sigmaSqIGb_r);                                    // BJ
    double *phiUnifa = REAL(phiUnifa_r);                                        // BJ
    double *phiUnifb = REAL(phiUnifb_r);                                        // BJ                                                  // BJ
    double *nuUnifa = REAL(nuUnifa_r);                                          // BJ              
    double *nuUnifb = REAL(nuUnifb_r);                                          // BJ
    double *tauSqIGa = REAL(tauSqIGa_r);                                        // BJ 
    double *tauSqIGb = REAL(tauSqIGb_r);                                        // BJ 
    double crossphiUnifa, crossphiUnifb, rhoUnifa = -1, rhoUnifb = 1;           // BJ
    
    int nSamples = INTEGER(nSamples_r)[0];
    int nThreads = INTEGER(nThreads_r)[0];
    int verbose = INTEGER(verbose_r)[0];
    int nReport = INTEGER(nReport_r)[0];
    
#ifdef _OPENMP
    omp_set_num_threads(nThreads);
#else
    if(nThreads > 1){
      warning("n.omp.threads > %i, but source not compiled with OpenMP support.", nThreads);
      nThreads = 1;
    }
#endif
    
    if(verbose){
      Rprintf("----------------------------------------\n");
      Rprintf("\tModel description\n");
      Rprintf("----------------------------------------\n");
      Rprintf("Spanning Tree-Based NNGP Non-separable Latent model fit with %i observations.\n\n", n);
      Rprintf("Number of covariates %i (including intercept if specified).\n\n", p);
      Rprintf("Using the %s spatial correlation model.\n\n", corName.c_str());
      Rprintf("Using %i nearest neighbors.\n\n", m);
      Rprintf("Number of MCMC samples %i.\n\n", nSamples);
      Rprintf("Priors and hyperpriors:\n");
      Rprintf("\tbeta flat.\n");
      for (int r = 0; r < q; r++) {                                             // BJ
        Rprintf("\tsigma.sq[%i] IG hyperpriors shape=%.5f and scale=%.5f\n",    // BJ
                r+1, sigmaSqIGa[r], sigmaSqIGb[r]);                             // BJ
        Rprintf("\ttau.sq[%i] IG hyperpriors shape=%.5f and scale=%.5f\n",      // BJ
                r+1, tauSqIGa[r], tauSqIGb[r]);                                 // BJ
        Rprintf("\tphi[%i] Unif hyperpriors a=%.5f and b=%.5f\n",               // BJ
                r+1, phiUnifa[r], phiUnifb[r]);                                 // BJ
        if(corName == "matern"){                                                // BJ
          Rprintf("\tnu[%i] Unif hyperpriors a=%.5f and b=%.5f\n",              // BJ
                  r+1, nuUnifa[r], nuUnifb[r]);	                                // BJ
        }                                                                       // BJ
      }                                                                         // BJ
      Rprintf("\tphi(j->k) Unif hyperpriors a=min(phi[j],phi[k]) and b=max(phi[j],phi[k])\n"); // BJ
      Rprintf("\trho(j->k) Unif hyperpriors as a function of nu[j], nu[k], nu(j->k), phi[j], phi[k], phi(j->k)\n"); // BJ
      Rprintf("\tnu(j->k) = 0.5nu[j] + 0.5nu[k]\n");                            // BJ
      Rprintf("\tbased on validity of the full bivariate Matern model (Thm3 (c) of Gneiting, Kleiber, Schlather (2010)).\n"); // BJ

      #ifdef _OPENMP
      Rprintf("\nSource compiled with OpenMP support and model fit using %i thread(s).\n", nThreads);
#else
      Rprintf("\n\nSource not compiled with OpenMP support.\n");
#endif
    } 
    
    //starting	
    int nq = n*q;                                                               // BJ
    int pq = p*q;                                                               // BJ
    int qm1 = q-1;                                                              // BJ
    double *beta = (double *) R_alloc(pq, sizeof(double));                      // BJ
    double *betaj = (double *) R_alloc(p, sizeof(double));                      // BJ
    double *sigmaSq = (double *) R_alloc(q, sizeof(double));                    // BJ
    double *phi = (double *) R_alloc(q, sizeof(double));                        // BJ
    double *crossphi = (double *) R_alloc(qm1, sizeof(double));                 // BJ
    double *nu = (double *) R_alloc(q, sizeof(double));                         // BJ
    double *tauSq = (double *) R_alloc(q, sizeof(double));                      // BJ
    double *rho = (double *) R_alloc(qm1, sizeof(double));                      // BJ
    int *adjvec = INTEGER(adjmatStarting_r);                                    // BJ
    
    F77_NAME(dcopy)(&pq, REAL(betaStarting_r), &inc, beta, &inc);               // BJ
    F77_NAME(dcopy)(&q, REAL(sigmaSqStarting_r), &inc, sigmaSq, &inc);          // BJ
    F77_NAME(dcopy)(&q, REAL(phiStarting_r), &inc, phi, &inc);                  // BJ
    F77_NAME(dcopy)(&qm1, REAL(crossphiStarting_r), &inc, crossphi, &inc);      // BJ
    if(corName == "matern"){
      F77_NAME(dcopy)(&q, REAL(nuStarting_r), &inc, nu, &inc);                  // BJ
    } else if (corName == "gaussian") {                                         // BJ
      for(int o = 0; o < q; o++) { nu[o] = 1e+17; }                             // BJ: for rholim, 1e+17 role of infinity
    } else if (corName == "exponential") {                                      // BJ
      for(int o = 0; o < q; o++) { nu[o] = 0.5; }                               // BJ: for rholim
    } else { zeros(nu, q); }                                                    // BJ
    F77_NAME(dcopy)(&q, REAL(tauSqStarting_r), &inc, tauSq, &inc);              // BJ
    F77_NAME(dcopy)(&qm1, REAL(rhoStarting_r), &inc, rho, &inc);                // BJ
    
    //tuning and fixed
    double *sigmaSqTuning = REAL(sigmaSqTuning_r);                              // BJ
    double *phiTuning = REAL(phiTuning_r);                                      // BJ
    double *nuTuning = (double *) R_alloc(q, sizeof(double));                   // BJ
    if(corName == "matern"){                                                    // BJ
      F77_NAME(dcopy)(&q, REAL(nuTuning_r), &inc, nuTuning, &inc);              // BJ
    } else { zeros(nuTuning, q); }                                              // BJ
    double *rhoTuning = REAL(rhoTuning_r);                                      // BJ
    double *crossphiTuning = REAL(crossphiTuning_r);                            // BJ
    
    //return stuff  
    SEXP wSamples_r, betaSamples_r, sigmaSqSamples_r, phiSamples_r,             // BJ
    crossphiSamples_r, tauSqSamples_r, rhoSamples_r, nuSamples_r;               // BJ
    PROTECT(wSamples_r = allocMatrix(REALSXP, nq, nSamples)); nProtect++;       // BJ
    PROTECT(betaSamples_r = allocMatrix(REALSXP, pq, nSamples)); nProtect++;    // BJ
    PROTECT(sigmaSqSamples_r = allocMatrix(REALSXP, q, nSamples)); nProtect++;  // BJ
    PROTECT(phiSamples_r = allocMatrix(REALSXP, q, nSamples)); nProtect++;      // BJ
    PROTECT(crossphiSamples_r = allocMatrix(REALSXP, qm1, nSamples)); nProtect++;// BJ
    if (corName == "matern") {
      PROTECT(nuSamples_r = allocMatrix(REALSXP, q, nSamples)); nProtect++;     // BJ
    }
    PROTECT(tauSqSamples_r = allocMatrix(REALSXP, q, nSamples)); nProtect++;    // BJ
    PROTECT(rhoSamples_r = allocMatrix(REALSXP, qm1, nSamples)); nProtect++;    // BJ
    
    //miscellaneous
    int nIndx = static_cast<int>(static_cast<double>(m*(m-1)/2+(n-m)*m+(m+1)*n)); // BJ
    int twomp1 = 2*m+1;                                                         // BJ
    int twomp1sq = pow(twomp1,2);                                               // BJ
    double *B = (double *) R_alloc(nIndx*q, sizeof(double));                    // BJ, change: q times longer
    double *F = (double *) R_alloc(nq, sizeof(double));                         // BJ, change: q times longer
    double *BCand = (double *) R_alloc(nIndx*q, sizeof(double));                // BJ, change: q times longer
    double *FCand = (double *) R_alloc(nq, sizeof(double));                     // BJ, change: q times longer
    double *c =(double *) R_alloc(twomp1*nThreads, sizeof(double));             // BJ
    double *C = (double *) R_alloc(twomp1sq*nThreads, sizeof(double));          // BJ

    double *sigmaSqCand = (double *) R_alloc(q, sizeof(double));                // BJ
    double *phiCand = (double *) R_alloc(q, sizeof(double));                    // BJ
    double *crossphiCand = (double *) R_alloc(qm1, sizeof(double));             // BJ
    double *nuCand = (double *) R_alloc(q, sizeof(double));                     // BJ
    F77_NAME(dcopy)(&q, nu, &inc, nuCand, &inc);                                // BJ
    double *rhoCand = (double *) R_alloc(qm1, sizeof(double));                  // BJ
    double crossphiCandUnifa, crossphiCandUnifb,
    rhoCandUnifa = -1, rhoCandUnifb = 1;                                        // BJ
    double crossnu, crossnuCand, infimum;                                       // BJ

    double logPostCand, logPostCurrent, logDet;
    double logPriorJacobianCand, logPriorJacobianCurrent;                       // BJ
    int accept = 0, batchAccept = 0, status = 0;
    int jj, kk, ll, pp = p*p;
    double *tmp_pp = (double *) R_alloc(pp, sizeof(double));
    double *tmp_p = (double *) R_alloc(p, sizeof(double));
    double *tmp_p2 = (double *) R_alloc(p, sizeof(double));
    double *tmp_n = (double *) R_alloc(n, sizeof(double)); zeros(tmp_n, n);
    double *tmp_n_parent = (double *) R_alloc(n, sizeof(double));               // BJ
    double *tmp_zero = (double *) R_alloc(n, sizeof(double));                   // BJ
    zeros(tmp_zero, n);                                                         // BJ
    double *XtX = (double *) R_alloc(pp, sizeof(double));
    F77_NAME(dgemm)(ytran, ntran, &p, &p, &n, &one, X, &n, X, &n, &zero, XtX, &p FCONE FCONE);
    double *w = (double *) R_alloc(nq, sizeof(double)); zeros(w, nq);           // BJ
    double a, v, b, e, mu, var, aij;
    double ac, a_ss, a_th, vc;                                                  // BJ

    int MeIndx, ParentIndx, ChildIndx;                                          // BJ

    if(verbose){
      Rprintf("----------------------------------------\n");
      Rprintf("\t\tSampling\n");
      Rprintf("----------------------------------------\n");
#ifdef Win32
      R_FlushConsole();
#endif
    }

    GetRNGstate();
    
    updateBF_parent_NS(B, F, c, C, coords,
                       nnIndx, nnIndxLU,
                       nnIndxParent, nnIndxLUParent, nnIndxLUAll,
                       n, m, nIndx, twomp1, twomp1sq, q,
                       sigmaSq, phi, nu,
                       rho, crossphi, adjvec,
                       covModel, nThreads, nuUnifb);
    
    for(s = 0; s < nSamples; s++){
      ParentIndx = 0;                                                           // BJ: for MeIndx = 0, does not affect results
      for (MeIndx = 0; MeIndx < q; MeIndx++) {                                  // BJ
        if (MeIndx > 0) {                                                       // BJ
          ParentIndx = which(1, &adjvec[q*MeIndx], q);                          // BJ
        }

        ///////////////
        //update w
        ///////////////
        if (MeIndx == 0) {
          for(i = 0; i < n; i++){
            // BJ: [case 1] conditional distribution of w_MeIdx(s_i) 
            e = 0;
            for(j = 0; j < nnIndxLU[n+i]; j++){                                 // BJ
              e += B[nIndx*MeIndx+nnIndxLU[i]+j]*w[n*MeIndx+nnIndx[nnIndxLU[i]+j]]; // BJ
            } 
            
            // BJ: [case 2] w_MeIndx(s_i) is a parent of w_MeIndx(s_jj)
            a = 0;
            v = 0;
            if(uIndxLU[n+i] > 0){//is i a neighbor for anybody
              for(j = 0; j < uIndxLU[n+i]; j++){//how many location have i as a neighbor
                b = 0;
                //now the neighbors for the jth location who has i as a neighbor
                jj = uIndx[uIndxLU[i]+j]; //jj is the index of the jth location who has i as a neighbor
                for(k = 0; k < nnIndxLU[n+jj]; k++){// these are the neighbors of the jjth location // BJ
                  kk = nnIndx[nnIndxLU[jj]+k];// kk is the index for the jth locations neighbors // BJ
                  if(kk != i) {//if the neighbor of jj is not i                 // BJ
                    b += B[nIndx*MeIndx+nnIndxLU[jj]+k]*w[n*MeIndx+kk];//covariance between jj and kk and the random effect of kk // BJ
                  }
                }
                aij = w[n*MeIndx+jj] - b;
                a += B[nIndx*MeIndx+nnIndxLU[jj]+uiIndx[uIndxLU[i]+j]]*aij/F[n*MeIndx+jj]; // BJ
                v += pow(B[nIndx*MeIndx+nnIndxLU[jj]+uiIndx[uIndxLU[i]+j]],2)/F[n*MeIndx+jj]; // BJ
              }
            }
            
            // BJ: [case 3] w_MeIndx(s_i) is a parent of w_ChildIndx(s_jj) (s_jj may be s_i)
            ac = 0;                                                             // BJ
            vc = 0;                                                             // BJ
            if(cIndxLU[q+MeIndx] > 0){                                          // BJ: does MeIndx have children
              for(j = 0; j < cIndxLU[q+MeIndx]; j++){                           // BJ: how many children does MeIndx have
                ChildIndx = cIndx[cIndxLU[MeIndx]+j];                           // BJ: ChildIndx is the child variable of MeIndx
                
                if(uuIndxLU[n+i] > 0){//is i a neighbor for anybody
                  for(l = 0; l < uuIndxLU[n+i]; l++){//how many location have i as a neighbor
                    b = 0;
                    ll = uuIndx[uuIndxLU[i]+l]; //ll is the index of the lth location who has i as a neighbor (BJ: i*)
                    for(k = 0; k < nnIndxLUAll[n+ll]; k++){// these are the neighbors of the llth location
                      if (k < m+1) {
                        kk = nnIndxParent[nnIndxLUParent[ll]+k];// kk is the index for the llth locations neighbors
                        if (kk != i) {
                          b += B[nIndx*ChildIndx+nnIndxLUAll[ll]+k]*w[n*MeIndx+kk];
                        } 
                      } else {
                        k2 = k-m-1;
                        b += B[nIndx*ChildIndx+nnIndxLUAll[ll]+k]*w[n*ChildIndx+nnIndx[nnIndxLU[ll]+k2]];
                      }
                    }
                    aij = w[n*ChildIndx+ll] - b;
                    ac += B[nIndx*ChildIndx+nnIndxLUAll[ll]+uuiIndx[uuIndxLU[i]+l]]*aij/F[n*ChildIndx+ll];
                    vc += pow(B[nIndx*ChildIndx+nnIndxLUAll[ll]+uuiIndx[uuIndxLU[i]+l]],2)/F[n*ChildIndx+ll];
                  }
                }
              }
            }
            
            mu = (y[n*MeIndx+i] - F77_NAME(ddot)(&p, &X[i], &n, &beta[p*MeIndx], &inc))/tauSq[MeIndx] +
              e/F[n*MeIndx+i] + a + ac;
            var = 1.0/(1.0/tauSq[MeIndx] + 1.0/F[n*MeIndx+i] + v + vc);         // BJ
            
            w[n*MeIndx+i] = rnorm(mu*var, sqrt(var));                           // BJ
            // // BJ: check
            // Rprintf("mean(w_%i[%i])=%f\n", MeIndx+1, i+1, mu*var);
            // Rprintf("var(w_%i[%i])=%f\n", MeIndx+1, i+1, var);
          }
        } else {
          for(i = 0; i < n; i++){
            // BJ: [case 1] conditional distribution of w_MeIdx(s_i) 
            e = 0;
            for(j = 0; j < nnIndxLUAll[n+i]; j++){                              // BJ
              if (j < m+1) {                                                    // BJ
                e += B[nIndx*MeIndx+nnIndxLUAll[i]+j]*w[n*ParentIndx+nnIndxParent[nnIndxLUParent[i]+j]]; // BJ
              } else {                                                          // BJ
                j2 = j-m-1;                                                     // BJ
                e += B[nIndx*MeIndx+nnIndxLUAll[i]+j]*w[n*MeIndx+nnIndx[nnIndxLU[i]+j2]]; // BJ
              }                                                                 // BJ
            }                                                                   // BJ
            
            // BJ: [case 2] w_MeIndx(s_i) is a parent of w_MeIndx(s_jj)
            a = 0;
            v = 0;
            if(uIndxLU[n+i] > 0){//is i a neighbor for anybody
              for(j = 0; j < uIndxLU[n+i]; j++){//how many location have i as a neighbor
                b = 0;
                //now the neighbors for the jth location who has i as a neighbor
                jj = uIndx[uIndxLU[i]+j]; //jj is the index of the jth location who has i as a neighbor
                for(k = 0; k < nnIndxLUAll[n+jj]; k++){// these are the neighbors of the jjth location // BJ
                  if (k < m+1) {                                                // BJ
                    b += B[nIndx*MeIndx+nnIndxLUAll[jj]+k]*w[n*ParentIndx+nnIndxParent[nnIndxLUParent[jj]+k]]; // BJ
                  } else {                                                      // BJ
                    k2 = k-m-1;                                                 // BJ
                    kk = nnIndx[nnIndxLU[jj]+k2];// kk is the index for the jth locations neighbors // BJ
                    if(kk != i) {//if the neighbor of jj is not i               // BJ
                      b += B[nIndx*MeIndx+nnIndxLUAll[jj]+k]*w[n*MeIndx+kk];//covariance between jj and kk and the random effect of kk // BJ
                    }
                  }
                }
                aij = w[n*MeIndx+jj] - b;
                a += B[nIndx*MeIndx+nnIndxLUAll[jj]+uiIndx[uIndxLU[i]+j]+m+1]*aij/F[n*MeIndx+jj]; // BJ: m+1 added as location i should come from the same variable
                v += pow(B[nIndx*MeIndx+nnIndxLUAll[jj]+uiIndx[uIndxLU[i]+j]+m+1],2)/F[n*MeIndx+jj]; // BJ
              }
            }
            
            // BJ: [case 3] w_MeIndx(s_i) is a parent of w_ChildIndx(s_jj) (s_jj may be s_i)
            ac = 0;                                                             // BJ
            vc = 0;                                                             // BJ
            if(cIndxLU[q+MeIndx] > 0){                                          // BJ: does MeIndx have children
              for(j = 0; j < cIndxLU[q+MeIndx]; j++){                           // BJ: how many children does MeIndx have
                ChildIndx = cIndx[cIndxLU[MeIndx]+j];                           // BJ: ChildIndx is the child variable of MeIndx
                
                if(uuIndxLU[n+i] > 0){//is i a neighbor for anybody
                  for(l = 0; l < uuIndxLU[n+i]; l++){//how many location have i as a neighbor
                    b = 0;
                    ll = uuIndx[uuIndxLU[i]+l]; //ll is the index of the lth location who has i as a neighbor (BJ: i*)
                    for(k = 0; k < nnIndxLUAll[n+ll]; k++){// these are the neighbors of the llth location
                      if (k < m+1) {
                        kk = nnIndxParent[nnIndxLUParent[ll]+k];// kk is the index for the llth locations neighbors
                        if (kk != i) {
                          b += B[nIndx*ChildIndx+nnIndxLUAll[ll]+k]*w[n*MeIndx+kk];
                        } 
                      } else {
                        k2 = k-m-1;
                        b += B[nIndx*ChildIndx+nnIndxLUAll[ll]+k]*w[n*ChildIndx+nnIndx[nnIndxLU[ll]+k2]];
                      }
                    }
                    aij = w[n*ChildIndx+ll] - b;
                    ac += B[nIndx*ChildIndx+nnIndxLUAll[ll]+uuiIndx[uuIndxLU[i]+l]]*aij/F[n*ChildIndx+ll];
                    vc += pow(B[nIndx*ChildIndx+nnIndxLUAll[ll]+uuiIndx[uuIndxLU[i]+l]],2)/F[n*ChildIndx+ll];
                  }
                }
              }
            }
            
            mu = (y[n*MeIndx+i] - F77_NAME(ddot)(&p, &X[i], &n, &beta[p*MeIndx], &inc))/tauSq[MeIndx] +
              e/F[n*MeIndx+i] + a + ac;
            var = 1.0/(1.0/tauSq[MeIndx] + 1.0/F[n*MeIndx+i] + v + vc);         // BJ
            
            w[n*MeIndx+i] = rnorm(mu*var, sqrt(var));                           // BJ
            // // BJ: check
            // Rprintf("mean(w_%i[%i])=%f\n", MeIndx+1, i+1, mu*var);
            // Rprintf("var(w_%i[%i])=%f\n", MeIndx+1, i+1, var);
          }
        }

        ///////////////
        //update beta
        ///////////////
        for(i = 0; i < n; i++){
          tmp_n[i] = (y[n*MeIndx+i] - w[n*MeIndx+i])/tauSq[MeIndx];             // BJ
        }
        F77_NAME(dgemv)(ytran, &n, &p, &one, X, &n, tmp_n, &inc, &zero, tmp_p, &inc FCONE);

        for(i = 0; i < pp; i++){
          tmp_pp[i] = XtX[i]/tauSq[MeIndx];                                     // BJ
        }
        F77_NAME(dpotrf)(lower, &p, tmp_pp, &p, &info FCONE);
        if(info != 0){error("c++ error: dpotrf failed\n");}
        F77_NAME(dpotri)(lower, &p, tmp_pp, &p, &info FCONE);
        if(info != 0){error("c++ error: dpotri failed\n");}
        F77_NAME(dsymv)(lower, &p, &one, tmp_pp, &p, tmp_p, &inc, &zero, tmp_p2, &inc FCONE);
        F77_NAME(dpotrf)(lower, &p, tmp_pp, &p, &info FCONE);
        if(info != 0){error("c++ error: dpotrf failed\n");}
        mvrnorm(betaj, tmp_p2, tmp_pp, p);                                      // BJ
        F77_NAME(dcopy)(&p, betaj, &inc, &beta[p*MeIndx], &inc);                // BJ
        // // BJ: check
        // Rprintf("mean(beta_%i)=(%f, %f)\n", MeIndx+1, tmp_p2[0], tmp_p2[1]);
        // Rprintf("var(beta_%i)=(%f, %f, %f, %f)\n", MeIndx+1, tmp_pp[0], tmp_pp[1], tmp_pp[2], tmp_pp[3]);

        /////////////////////
        //update tau^2
        /////////////////////
        for(i = 0; i < n; i++){
          tmp_n[i] = y[n*MeIndx+i] - w[n*MeIndx+i] -
            F77_NAME(ddot)(&p, &X[i], &n, betaj, &inc);                         // BJ
        }

        if (tauSq[MeIndx] != 0) {                                               // BJ: no nugget if tauSq starting = 0 for latent model due to no tuning value for tauSq
          tauSq[MeIndx] = 1.0/rgamma(tauSqIGa[MeIndx]+n/2.0, 1.0/(tauSqIGb[MeIndx]+0.5*F77_NAME(ddot)(&n, tmp_n, &inc, tmp_n, &inc))); // BJ
        }
      }

      ///////////////
      //update theta
      ///////////////
      //current
      updateBF_parent_NS(B, F, c, C, coords,
                         nnIndx, nnIndxLU,
                         nnIndxParent, nnIndxLUParent, nnIndxLUAll,
                         n, m, nIndx, twomp1, twomp1sq, q,
                         sigmaSq, phi, nu,
                         rho, crossphi, adjvec,
                         covModel, nThreads, nuUnifb);
                      
      ParentIndx = 0;                                                           // BJ: for MeIndx = 0, does not affect results
      a_th = 0;                                                                 // BJ
      logDet = 0;                                                               // BJ
      logPriorJacobianCurrent = 0;                                              // BJ
      logPriorJacobianCand = 0;                                                 // BJ
      for (MeIndx = 0; MeIndx < q; MeIndx++) {                                  // BJ
        // BJ: propose candidate and compute both Jacobian at once
        // BJ: sigmaSq
        logPriorJacobianCurrent += -1.0*(1.0+sigmaSqIGa[MeIndx])*log(sigmaSq[MeIndx]) - 
          sigmaSqIGb[MeIndx]/sigmaSq[MeIndx]+log(sigmaSq[MeIndx]);              // BJ
        
        sigmaSqCand[MeIndx] = exp(rnorm(log(sigmaSq[MeIndx]), sigmaSqTuning[MeIndx])); // BJ
        logPriorJacobianCand += -1.0*(1.0+sigmaSqIGa[MeIndx])*log(sigmaSqCand[MeIndx])-
          sigmaSqIGb[MeIndx]/sigmaSqCand[MeIndx]+log(sigmaSqCand[MeIndx]);      // BJ
        
        // BJ: phi
        logPriorJacobianCurrent += log(phi[MeIndx] - phiUnifa[MeIndx]) +
          log(phiUnifb[MeIndx] - phi[MeIndx]);                                  // BJ

        phiCand[MeIndx] = logitInv(rnorm(logit(phi[MeIndx], phiUnifa[MeIndx], phiUnifb[MeIndx]),
                                         phiTuning[MeIndx]), phiUnifa[MeIndx], phiUnifb[MeIndx]); // BJ
        logPriorJacobianCand += log(phiCand[MeIndx] - phiUnifa[MeIndx]) +
          log(phiUnifb[MeIndx] - phiCand[MeIndx]);                              // BJ

        // BJ: nu
        if (corName == "matern"){
          logPriorJacobianCurrent += log(nu[MeIndx] - nuUnifa[MeIndx]) +
            log(nuUnifb[MeIndx] - nu[MeIndx]);                                  // BJ

          nuCand[MeIndx] = logitInv(rnorm(logit(nu[MeIndx], nuUnifa[MeIndx], nuUnifb[MeIndx]),
                                          nuTuning[MeIndx]), nuUnifa[MeIndx], nuUnifb[MeIndx]); // BJ
          logPriorJacobianCand += log(nuCand[MeIndx] - nuUnifa[MeIndx]) +
            log(nuUnifb[MeIndx] - nuCand[MeIndx]);                              // BJ
        }

        if (MeIndx == 0) {
#ifdef _OPENMP
#pragma omp parallel for private (e, j, b) reduction(+:a_th, logDet)
#endif

          for(i = 0; i < n; i++){
            e = 0;
            for(j = 0; j < nnIndxLU[n+i]; j++){                                 // BJ
              e += B[nIndx*MeIndx+nnIndxLU[i]+j]*w[n*MeIndx+nnIndx[nnIndxLU[i]+j]]; // BJ
            }                                                                   // BJ
            b = w[n*MeIndx+i] - e;                                              // BJ
            a_th += b*b/F[n*MeIndx+i];                                          // BJ: a_th sum over MeIndx and i
            logDet += log(F[n*MeIndx+i]);                                       // BJ: logDet sum over MeIndx and i
          }
        } else {
          ParentIndx = which(1, &adjvec[q*MeIndx], q);                          // BJ
          // BJ: propose candidate and compute both Jacobian at once
          // BJ: crossphi
          crossphiUnifa = std::min(phi[MeIndx], phi[ParentIndx]);               // BJ
          crossphiUnifb = std::max(phi[MeIndx], phi[ParentIndx]);               // BJ
          logPriorJacobianCurrent += log(crossphi[MeIndx-1] - crossphiUnifa) +         // BJ
            log(crossphiUnifb - crossphi[MeIndx-1]);                            // BJ

          crossphiCandUnifa = std::min(phiCand[MeIndx], phiCand[ParentIndx]);   // BJ
          crossphiCandUnifb = std::max(phiCand[MeIndx], phiCand[ParentIndx]);   // BJ
          crossphiCand[MeIndx-1] = logitInv(rnorm(logit(crossphi[MeIndx-1], crossphiUnifa, crossphiUnifb),
                                                  crossphiTuning[MeIndx-1]), crossphiCandUnifa, crossphiCandUnifb); // BJ
          logPriorJacobianCand += log(crossphiCand[MeIndx-1] - crossphiCandUnifa) +// BJ
            log(crossphiCandUnifb - crossphiCand[MeIndx-1]);                    // BJ

          // BJ: rho
          crossnu = 0.5*nu[MeIndx] + 0.5*nu[ParentIndx];                        // BJ
          rhoUnifb = rholim(crossphi[MeIndx-1], phi[MeIndx], phi[ParentIndx],
                            nu[MeIndx], nu[ParentIndx], crossnu);               // BJ
          rhoUnifa = -1.0*rhoUnifb;                                             // BJ
          logPriorJacobianCurrent += log(rho[MeIndx-1] - rhoUnifa) +                   // BJ
            log(rhoUnifb - rho[MeIndx-1]);                                      // BJ

          crossnuCand = 0.5*nuCand[MeIndx] + 0.5*nuCand[ParentIndx];            // BJ
          rhoCandUnifb = rholim(crossphiCand[MeIndx-1], phiCand[MeIndx],
                                phiCand[ParentIndx], nuCand[MeIndx],
                                nuCand[ParentIndx], crossnuCand);               // BJ
          rhoCandUnifa = -1.0*rhoCandUnifb;                                     // BJ
          rhoCand[MeIndx-1] = logitInv(rnorm(logit(rho[MeIndx-1], rhoUnifa, rhoUnifb),
                                             rhoTuning[MeIndx-1]), rhoCandUnifa, rhoCandUnifb); // BJ
          logPriorJacobianCand += log(rhoCand[MeIndx-1] - rhoCandUnifa) +       // BJ
            log(rhoCandUnifb - rhoCand[MeIndx-1]);                              // BJ

#ifdef _OPENMP
#pragma omp parallel for private (e, j, b) reduction(+:a_th, logDet)
#endif

          for(i = 0; i < n; i++){
            e = 0;
            for(j = 0; j < nnIndxLUAll[n+i]; j++){                              // BJ
              if (j < m+1) {                                                    // BJ
                e += B[nIndx*MeIndx+nnIndxLUAll[i]+j]*w[n*ParentIndx+nnIndxParent[nnIndxLUParent[i]+j]]; // BJ
              } else {                                                          // BJ
                j2 = j-m-1;                                                     // BJ
                e += B[nIndx*MeIndx+nnIndxLUAll[i]+j]*w[n*MeIndx+nnIndx[nnIndxLU[i]+j2]]; // BJ
              }                                                                 // BJ
            }                                                                   // BJ
            b = w[n*MeIndx+i] - e;                                              // BJ
            a_th += b*b/F[n*MeIndx+i];                                          // BJ: a_th sum over MeIndx and i
            logDet += log(F[n*MeIndx+i]);                                       // BJ: logDet sum over MeIndx and i
          }
        }
      }

      logPostCurrent = -0.5*logDet - 0.5*a_th + logPriorJacobianCurrent;        // BJ

      //candidate
      updateBF_parent_NS(BCand, FCand, c, C, coords,
                         nnIndx, nnIndxLU,
                         nnIndxParent, nnIndxLUParent, nnIndxLUAll,
                         n, m, nIndx, twomp1, twomp1sq, q,
                         sigmaSqCand, phiCand, nuCand,
                         rhoCand, crossphiCand, adjvec,
                         covModel, nThreads, nuUnifb);

      ParentIndx = 0;                                                           // BJ: for MeIndx = 0, does not affect results
      a_th = 0;                                                                 // BJ
      logDet = 0;                                                               // BJ
      for (MeIndx = 0; MeIndx < q; MeIndx++) {                                  // BJ
        if (MeIndx == 0) {
#ifdef _OPENMP
#pragma omp parallel for private (e, j, b) reduction(+:a_th, logDet)
#endif

          for(i = 0; i < n; i++){
            e = 0;
            for(j = 0; j < nnIndxLU[n+i]; j++){                                 // BJ
              e += BCand[nIndx*MeIndx+nnIndxLU[i]+j]*w[n*MeIndx+nnIndx[nnIndxLU[i]+j]]; // BJ
            }                                                                   // BJ
            b = w[n*MeIndx+i] - e;                                              // BJ
            a_th += b*b/FCand[n*MeIndx+i];                                      // BJ: a_th sum over MeIndx and i
            logDet += log(FCand[n*MeIndx+i]);                                   // BJ: logDet sum over MeIndx and i
          }
        } else {
          ParentIndx = which(1, &adjvec[q*MeIndx], q);                          // BJ
#ifdef _OPENMP
#pragma omp parallel for private (e, j, b) reduction(+:a_th, logDet)
#endif

          for(i = 0; i < n; i++){
            e = 0;
            for(j = 0; j < nnIndxLUAll[n+i]; j++){                              // BJ
              if (j < m+1) {                                                    // BJ
                e += BCand[nIndx*MeIndx+nnIndxLUAll[i]+j]*w[n*ParentIndx+nnIndxParent[nnIndxLUParent[i]+j]]; // BJ
              } else {                                                          // BJ
                j2 = j-m-1;                                                     // BJ
                e += BCand[nIndx*MeIndx+nnIndxLUAll[i]+j]*w[n*MeIndx+nnIndx[nnIndxLU[i]+j2]]; // BJ
              }                                                                 // BJ
            }                                                                   // BJ
            b = w[n*MeIndx+i] - e;                                              // BJ
            a_th += b*b/FCand[n*MeIndx+i];                                      // BJ: a_th sum over MeIndx and i
            logDet += log(FCand[n*MeIndx+i]);                                   // BJ: logDet sum over MeIndx and i
          }
        }
      }

      logPostCand = -0.5*logDet - 0.5*a_th + logPriorJacobianCand;              // BJ

      if(runif(0.0,1.0) <= exp(logPostCand - logPostCurrent)){

        std::swap(sigmaSqCand, sigmaSq);                                        // BJ
        std::swap(phiCand, phi);                                                // BJ
        if(corName == "matern"){
          std::swap(nuCand, nu);                                                // BJ
        }
        std::swap(rhoCand, rho);                                                // BJ
        std::swap(crossphiCand, crossphi);                                      // BJ

        std::swap(BCand, B);
        std::swap(FCand, F);

        accept++;
        batchAccept++;
      }

      //save samples
      F77_NAME(dcopy)(&pq, beta, &inc, &REAL(betaSamples_r)[s*pq], &inc);       // BJ
      F77_NAME(dcopy)(&q, sigmaSq, &inc, &REAL(sigmaSqSamples_r)[s*q], &inc);   // BJ
      F77_NAME(dcopy)(&q, phi, &inc, &REAL(phiSamples_r)[s*q], &inc);           // BJ
      F77_NAME(dcopy)(&qm1, crossphi, &inc, &REAL(crossphiSamples_r)[s*qm1], &inc); // BJ
      if (corName == "matern") {                                                // BJ
        F77_NAME(dcopy)(&q, nu, &inc, &REAL(nuSamples_r)[s*q], &inc);           // BJ
      }                                                                         // BJ
      F77_NAME(dcopy)(&qm1, rho, &inc, &REAL(rhoSamples_r)[s*qm1], &inc);       // BJ
      F77_NAME(dcopy)(&q, tauSq, &inc, &REAL(tauSqSamples_r)[s*q], &inc);       // BJ
      F77_NAME(dcopy)(&nq, w, &inc, &REAL(wSamples_r)[s*nq], &inc);             // BJ

      //report
      if(status == nReport){
        if(verbose){
          Rprintf("Sampled: %i of %i, %3.2f%%\n", s, nSamples, 100.0*s/nSamples);
          Rprintf("Report interval Metrop. Acceptance rate: %3.2f%%\n", 100.0*batchAccept/nReport);
          Rprintf("Overall Metrop. Acceptance rate: %3.2f%%\n", 100.0*accept/s);
          Rprintf("-------------------------------------------------\n");
#ifdef Win32
          R_FlushConsole();
#endif
        }
        batchAccept = 0;
        status = 0;
      }

      status++;

      R_CheckUserInterrupt();
    }
    
    if(verbose){
      Rprintf("Sampled: %i of %i, %3.2f%%\n", s, nSamples, 100.0);
      Rprintf("Report interval Metrop. Acceptance rate: %3.2f%%\n", 100.0*batchAccept/nReport);
      Rprintf("Overall Metrop. Acceptance rate: %3.2f%%\n", 100.0*accept/nSamples);
      Rprintf("-------------------------------------------------\n");
#ifdef Win32
      R_FlushConsole();
#endif
    }
    
    PutRNGstate();
    
    //make return object
    SEXP result_r, resultName_r;
    int nResultListObjs = 7;
    if (corName == "matern") { nResultListObjs += 1;}
    
    PROTECT(result_r = allocVector(VECSXP, nResultListObjs)); nProtect++;
    PROTECT(resultName_r = allocVector(VECSXP, nResultListObjs)); nProtect++;
    
    SET_VECTOR_ELT(result_r, 0, betaSamples_r);
    SET_VECTOR_ELT(resultName_r, 0, mkChar("p.beta.samples"));
    
    SET_VECTOR_ELT(result_r, 1, tauSqSamples_r);
    SET_VECTOR_ELT(resultName_r, 1, mkChar("p.tausq.samples"));
    
    SET_VECTOR_ELT(result_r, 2, sigmaSqSamples_r);
    SET_VECTOR_ELT(resultName_r, 2, mkChar("p.sigmasq.samples"));
    
    SET_VECTOR_ELT(result_r, 3, phiSamples_r);
    SET_VECTOR_ELT(resultName_r, 3, mkChar("p.phi.samples"));
    
    SET_VECTOR_ELT(result_r, 4, crossphiSamples_r);
    SET_VECTOR_ELT(resultName_r, 4, mkChar("p.crossphi.samples"));
    
    SET_VECTOR_ELT(result_r, 5, rhoSamples_r);
    SET_VECTOR_ELT(resultName_r, 5, mkChar("p.rho.samples"));
    
    SET_VECTOR_ELT(result_r, 6, wSamples_r);
    SET_VECTOR_ELT(resultName_r, 6, mkChar("p.w.samples"));
    
    if (corName == "matern") {
      SET_VECTOR_ELT(result_r, 7, nuSamples_r);
      SET_VECTOR_ELT(resultName_r, 7, mkChar("p.nu.samples"));
    }
    namesgets(result_r, resultName_r);
    
    //unprotect
    UNPROTECT(nProtect);
    
    return(result_r);
  }
}