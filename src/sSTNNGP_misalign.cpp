#define USE_FC_LEN_T
#include <string>
#include <iostream>
#include "util.h"
#include "STNNGP.h"                                                             // BJ: for updateBF_parent which is used in sSTNNGP as well.

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

extern "C" {
  
  SEXP sSTNNGP_misalign(SEXP y_r, SEXP X_r, 
                        SEXP XtX_r,                                             // BJ
                        SEXP q_r,                                               // BJ
                        SEXP p_r, SEXP n_r, SEXP m_r, SEXP coords_r, 
                        SEXP nj_r, SEXP Y_missing_r,                            // BJ
                        SEXP covModel_r, SEXP nnIndx_r, SEXP nnIndxLU_r, 
                        SEXP nnIndxParent_r, SEXP nnIndxLUParent_r, SEXP nnIndxLUAll_r,  // BJ
                        SEXP uIndx_r, SEXP uIndxLU_r, SEXP uiIndx_r,   
                        SEXP uuIndx_r, SEXP uuIndxLU_r, SEXP uuiIndx_r,         // BJ
                        SEXP cIndx_r, SEXP cIndxLU_r,                           // BJ
                        SEXP sigmaSqIG_r, 
                        SEXP tauSqIGa_r, SEXP tauSqIGb_r,                       // BJ
                        SEXP phiUnif_r, SEXP nuUnif_r, 
                        SEXP rhoUnif_r,                                         // BJ
                        SEXP betaStarting_r, SEXP sigmaSqStarting_r, SEXP tauSqStarting_r, 
                        SEXP phiStarting_r, SEXP nuStarting_r, 
                        SEXP rhoStarting_r, SEXP adjmatStarting_r,              // BJ
                        SEXP wStarting_r,                                       // BJ
                        SEXP phiTuning_r, SEXP nuTuning_r, 
                        SEXP rhoTuning_r,                                       // BJ
                        SEXP nSamples_r, SEXP nThreads_r, SEXP verbose_r, SEXP nReport_r, 
                        SEXP savew_r){                                          // BJ
    
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
    double *XtX = REAL(XtX_r);                                                  // BJ
    int q = INTEGER(q_r)[0];                                                    // BJ
    int p = INTEGER(p_r)[0];
    int n = INTEGER(n_r)[0];
    int m = INTEGER(m_r)[0];
    double *coords = REAL(coords_r);
    int *nj = INTEGER(nj_r);                                                    // BJ
    int *missing = INTEGER(Y_missing_r);                                        // BJ
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
    double sigmaSqIGa = REAL(sigmaSqIG_r)[0]; double sigmaSqIGb = REAL(sigmaSqIG_r)[1];
    double phiUnifa = REAL(phiUnif_r)[0]; double phiUnifb = REAL(phiUnif_r)[1];
    double nuUnifa = 0, nuUnifb = 0;
    if(corName == "matern"){
      nuUnifa = REAL(nuUnif_r)[0]; nuUnifb = REAL(nuUnif_r)[1]; 
    }
    double *tauSqIGa = REAL(tauSqIGa_r);                                        // BJ 
    double *tauSqIGb = REAL(tauSqIGb_r);                                        // BJ 
    double rhoUnifa = REAL(rhoUnif_r)[0]; double rhoUnifb = REAL(rhoUnif_r)[1]; // BJ
    
    // for iterations
    int nSamples = INTEGER(nSamples_r)[0];
    int nThreads = INTEGER(nThreads_r)[0];
    int verbose = INTEGER(verbose_r)[0];
    int nReport = INTEGER(nReport_r)[0];
    int savew = INTEGER(savew_r)[0];
    
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
      Rprintf("Spanning Tree-Based NNGP Separable Latent model fit with %i observations.\n\n", n);
      Rprintf("Number of covariates %i (including intercept if specified).\n\n", p);
      Rprintf("Using the %s spatial correlation model.\n\n", corName.c_str());
      Rprintf("Using %i nearest neighbors.\n\n", m);
      Rprintf("Number of MCMC samples %i.\n\n", nSamples);
      Rprintf("Priors and hyperpriors:\n");
      Rprintf("\tbeta flat.\n");
      Rprintf("\tsigma.sq IG hyperpriors shape=%.5f and scale=%.5f\n", sigmaSqIGa, sigmaSqIGb);
      for (int r = 0; r < q; r++) {                                             // BJ
        Rprintf("\ttau.sq[%i] IG hyperpriors shape=%.5f and scale=%.5f\n",      // BJ
                r+1, tauSqIGa[r], tauSqIGb[r]);                                 // BJ
      }                                                                         // BJ
      Rprintf("\tphi Unif hyperpriors a=%.5f and b=%.5f\n", phiUnifa, phiUnifb);
      if(corName == "matern"){
        Rprintf("\tnu Unif hyperpriors a=%.5f and b=%.5f\n", nuUnifa, nuUnifb);	  
      }
      Rprintf("\trho Unif hyperpriors a=%.5f and b=%.5f\n", rhoUnifa, rhoUnifb);// BJ
#ifdef _OPENMP
      Rprintf("\nSource compiled with OpenMP support and model fit using %i thread(s).\n", nThreads);
#else
      Rprintf("\n\nSource not compiled with OpenMP support.\n");
#endif
    } 
    
    //parameters                                                                // BJ: remove tauSqIndx
    int nTheta, sigmaSqIndx, phiIndx, nuIndx;                                   // BJ
    
    if(corName != "matern"){                                                    // BJ
      nTheta = 2;//sigma^2, phi                                                 // BJ
      sigmaSqIndx = 0; phiIndx = 1;                                             // BJ
    }else{                                                                      // BJ
      nTheta = 3;//sigma^2, phi, nu                                             // BJ
      sigmaSqIndx = 0; phiIndx = 1; nuIndx = 2;                                 // BJ 
    }                                                                           // BJ
    
    //starting	
    int nq = n*q;                                                               // BJ
    int pq = p*q;                                                               // BJ
    int qm1 = q-1;                                                              // BJ
    double *beta = (double *) R_alloc(pq, sizeof(double));                      // BJ
    double *betaj = (double *) R_alloc(p, sizeof(double));                      // BJ
    double *theta = (double *) R_alloc(nTheta, sizeof(double));
    double *tauSq = (double *) R_alloc(q, sizeof(double));                      // BJ
    double *rho = (double *) R_alloc(qm1, sizeof(double));                      // BJ
    int *adjvec = INTEGER(adjmatStarting_r);                                    // BJ
    double *w = (double *) R_alloc(nq, sizeof(double));                         // BJ
    
    F77_NAME(dcopy)(&pq, REAL(betaStarting_r), &inc, beta, &inc);               // BJ
    theta[sigmaSqIndx] = REAL(sigmaSqStarting_r)[0];
    theta[phiIndx] = REAL(phiStarting_r)[0];
    if(corName == "matern"){
      theta[nuIndx] = REAL(nuStarting_r)[0];
    }
    F77_NAME(dcopy)(&q, REAL(tauSqStarting_r), &inc, tauSq, &inc);              // BJ
    F77_NAME(dcopy)(&qm1, REAL(rhoStarting_r), &inc, rho, &inc);                // BJ
    F77_NAME(dcopy)(&nq, REAL(wStarting_r), &inc, w, &inc);                     // BJ
    
    //tuning and fixed
    double *tuning = (double *) R_alloc(nTheta, sizeof(double));
    
    tuning[sigmaSqIndx] = 0; //not accessed
    tuning[phiIndx] = REAL(phiTuning_r)[0];
    if(corName == "matern"){
      tuning[nuIndx] = REAL(nuTuning_r)[0];
    }
    double *rhoTuning = REAL(rhoTuning_r);                                      // BJ
    
    //return stuff  
    SEXP betaSamples_r, thetaSamples_r, wSamples_r, tauSqSamples_r, rhoSamples_r; // BJ
    PROTECT(betaSamples_r = allocMatrix(REALSXP, pq, nSamples)); nProtect++;    // BJ
    PROTECT(thetaSamples_r = allocMatrix(REALSXP, nTheta, nSamples)); nProtect++; 
    if (savew) {
      PROTECT(wSamples_r = allocMatrix(REALSXP, nq, nSamples)); nProtect++;       // BJ
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
    
    double logPostCand, logPostCurrent, logDet;     
    int accept = 0, batchAccept = 0, status = 0;
    int jj, kk, ll, pp = p*p;
    double *tmp_pp = (double *) R_alloc(pp, sizeof(double));
    double *tmp_p = (double *) R_alloc(p, sizeof(double));
    double *tmp_p2 = (double *) R_alloc(p, sizeof(double));
    double *tmp_n = (double *) R_alloc(n, sizeof(double)); zeros(tmp_n, n);
    double *tmp_n_parent = (double *) R_alloc(n, sizeof(double));               // BJ
    double *bk = (double *) R_alloc(nThreads*(1.0+static_cast<int>(floor(nuUnifb))), sizeof(double));
    double *tmp_zero = (double *) R_alloc(n, sizeof(double));                   // BJ
    zeros(tmp_zero, n);                                                         // BJ
    double a, v, b, e, mu, var, aij, phiCand, nuCand = 0, nu = 0; 
    if (corName == "matern") { nu = theta[nuIndx]; }
    double ac, a_ss, a_th, vc;                                                  // BJ
    double *rhoCand = (double *) R_alloc(qm1, sizeof(double));                  // BJ
    
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
    
    updateBF_parent(B, F, c, C, coords,
                    nnIndx, nnIndxLU,
                    nnIndxParent, nnIndxLUParent, nnIndxLUAll,
                    n, m, nIndx, twomp1, twomp1sq, q,
                    theta[sigmaSqIndx], theta[phiIndx], nu,
                    rho,
                    covModel, bk, nuUnifb);
    
    for(s = 0; s < nSamples; s++){
      ParentIndx = 0;                                                           // BJ: for MeIndx = 0, does not affect results
      a_ss = 0;                                                                 // BJ
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
            
            // BJ: if MeIndx is missing at location i, impute w from a prior
            if (missing[n*MeIndx+i] == 1.0) {
              mu = e/F[n*MeIndx+i] + a + ac;                                    // BJ
              var = 1.0/(1.0/F[n*MeIndx+i] + v + vc);                           // BJ
            } else {
              mu = (y[n*MeIndx+i] - F77_NAME(ddot)(&p, &X[i], &n, &beta[p*MeIndx], &inc))/tauSq[MeIndx] +
                e/F[n*MeIndx+i] + a + ac;
              var = 1.0/(1.0/tauSq[MeIndx] + 1.0/F[n*MeIndx+i] + v + vc);       // BJ
            }
            
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
            
            // BJ: if MeIndx is missing at location i, impute w from a prior
            if (missing[n*MeIndx+i] == 1.0) {
              mu = e/F[n*MeIndx+i] + a + ac;                                    // BJ
              var = 1.0/(1.0/F[n*MeIndx+i] + v + vc);                           // BJ
            } else {
              mu = (y[n*MeIndx+i] - F77_NAME(ddot)(&p, &X[i], &n, &beta[p*MeIndx], &inc))/tauSq[MeIndx] +
                e/F[n*MeIndx+i] + a + ac;
              var = 1.0/(1.0/tauSq[MeIndx] + 1.0/F[n*MeIndx+i] + v + vc);       // BJ
            }
            
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
          if (missing[n*MeIndx+i] == 1.0) {                                     // BJ
            tmp_n[i] = 0;                                                       // BJ
          } else {
            tmp_n[i] = (y[n*MeIndx+i] - w[n*MeIndx+i])/tauSq[MeIndx];           // BJ
          }
        }
        F77_NAME(dgemv)(ytran, &n, &p, &one, X, &n, tmp_n, &inc, &zero, tmp_p, &inc FCONE);
        
        for(i = 0; i < pp; i++){
          tmp_pp[i] = XtX[pp*MeIndx+i]/tauSq[MeIndx];                           // BJ: XtX based on Xj having zeros for locations where variable j is not observed
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
          if (missing[n*MeIndx+i] == 1.0) {                                     // BJ
            tmp_n[i] = 0;                                                       // BJ
          } else {
            tmp_n[i] = y[n*MeIndx+i] - w[n*MeIndx+i] -
              F77_NAME(ddot)(&p, &X[i], &n, betaj, &inc);                       // BJ
          }
        }

        if (tauSq[MeIndx] != 0) {                                               // BJ: no nugget if tauSq starting = 0 for latent model due to no tuning value for tauSq
          tauSq[MeIndx] = 1.0/rgamma(tauSqIGa[MeIndx]+nj[MeIndx]/2.0, 1.0/(tauSqIGb[MeIndx]+0.5*F77_NAME(ddot)(&n, tmp_n, &inc, tmp_n, &inc))); // BJ
        }

        /////////////////////
        //update sigma^2
        /////////////////////
        if (MeIndx == 0) {
#ifdef _OPENMP
#pragma omp parallel for private (e, j, b) reduction(+:a_ss, logDet)
#endif

          for(i = 0; i < n; i++){
            e = 0;
            for(j = 0; j < nnIndxLU[n+i]; j++){                                 // BJ
              e += B[nIndx*MeIndx+nnIndxLU[i]+j]*w[n*MeIndx+nnIndx[nnIndxLU[i]+j]]; // BJ
            }                                                                   // BJ
            b = w[n*MeIndx+i] - e;                                              // BJ
            a_ss += b*b/F[n*MeIndx+i];                                          // BJ: a_ss sum over MeIndx and i
          }
        } else {
#ifdef _OPENMP
#pragma omp parallel for private (e, j, b) reduction(+:a_ss, logDet)
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
            a_ss += b*b/F[n*MeIndx+i];                                          // BJ: a_ss sum over MeIndx and i
          }
        }
      }
      theta[sigmaSqIndx] = 1.0/rgamma(sigmaSqIGa+nq/2.0, 1.0/(sigmaSqIGb+0.5*a_ss*theta[sigmaSqIndx])); // BJ
      
      updateBF_parent(B, F, c, C, coords,
                      nnIndx, nnIndxLU,
                      nnIndxParent, nnIndxLUParent, nnIndxLUAll,
                      n, m, nIndx, twomp1, twomp1sq, q,
                      theta[sigmaSqIndx], theta[phiIndx], nu,
                      rho,
                      covModel, bk, nuUnifb);
      
      ///////////////
      //update theta
      ///////////////
      //current
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
              e += B[nIndx*MeIndx+nnIndxLU[i]+j]*w[n*MeIndx+nnIndx[nnIndxLU[i]+j]]; // BJ
            }                                                                   // BJ
            b = w[n*MeIndx+i] - e;                                              // BJ
            a_th += b*b/F[n*MeIndx+i];                                          // BJ: a_th sum over MeIndx and i
            logDet += log(F[n*MeIndx+i]);                                       // BJ: logDet sum over MeIndx and i
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

      logPostCurrent = -0.5*logDet - 0.5*a_th;
      logPostCurrent += log(theta[phiIndx] - phiUnifa) + log(phiUnifb - theta[phiIndx]);
      if(corName == "matern"){
        logPostCurrent += log(theta[nuIndx] - nuUnifa) + log(nuUnifb - theta[nuIndx]);
      }
      // BJ
      for (o = 1; o < q; o++) {
        logPostCurrent += log(rho[o-1] - rhoUnifa) + log(rhoUnifb - rho[o-1]);
      }

      //candidate
      phiCand = logitInv(rnorm(logit(theta[phiIndx], phiUnifa, phiUnifb), tuning[phiIndx]), phiUnifa, phiUnifb);
      if(corName == "matern"){
        nuCand = logitInv(rnorm(logit(theta[nuIndx], nuUnifa, nuUnifb), tuning[nuIndx]), nuUnifa, nuUnifb);
      }
      // BJ
      for (o = 1; o < q; o++) {
        rhoCand[o-1] = logitInv(rnorm(logit(rho[o-1], rhoUnifa, rhoUnifb), rhoTuning[o-1]), rhoUnifa, rhoUnifb);
      }
      
      updateBF_parent(BCand, FCand, c, C, coords,
                      nnIndx, nnIndxLU,
                      nnIndxParent, nnIndxLUParent, nnIndxLUAll,
                      n, m, nIndx, twomp1, twomp1sq, q,
                      theta[sigmaSqIndx], phiCand, nuCand,
                      rhoCand,
                      covModel, bk, nuUnifb);

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

      logPostCand = -0.5*logDet - 0.5*a_th;
      logPostCand += log(phiCand - phiUnifa) + log(phiUnifb - phiCand);
      if(corName == "matern"){
        logPostCand += log(nuCand - nuUnifa) + log(nuUnifb - nuCand);
      }
      // BJ
      for (o = 1; o < q; o++) {
        logPostCand += log(rhoCand[o-1] - rhoUnifa) + log(rhoUnifb - rhoCand[o-1]);
      }

      if(runif(0.0,1.0) <= exp(logPostCand - logPostCurrent)){

        theta[phiIndx] = phiCand;
        if(corName == "matern"){
          std::swap(nuCand, nu);                                                // BJ
          theta[nuIndx] = nuCand;
        }
        std::swap(rhoCand, rho);                                                // BJ
        
        std::swap(BCand, B);
        std::swap(FCand, F);

        accept++;
        batchAccept++;
      }

      //save samples
      F77_NAME(dcopy)(&pq, beta, &inc, &REAL(betaSamples_r)[s*pq], &inc);       // BJ
      F77_NAME(dcopy)(&nTheta, theta, &inc, &REAL(thetaSamples_r)[s*nTheta], &inc);
      F77_NAME(dcopy)(&qm1, rho, &inc, &REAL(rhoSamples_r)[s*qm1], &inc);       // BJ
      F77_NAME(dcopy)(&q, tauSq, &inc, &REAL(tauSqSamples_r)[s*q], &inc);       // BJ
      if (savew) {
        F77_NAME(dcopy)(&nq, w, &inc, &REAL(wSamples_r)[s*nq], &inc);             // BJ
      }

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
    int nResultListObjs = 4;
    if (savew) { nResultListObjs += 1; }
    
    PROTECT(result_r = allocVector(VECSXP, nResultListObjs)); nProtect++;
    PROTECT(resultName_r = allocVector(VECSXP, nResultListObjs)); nProtect++;
    
    SET_VECTOR_ELT(result_r, 0, betaSamples_r);
    SET_VECTOR_ELT(resultName_r, 0, mkChar("p.beta.samples"));
    
    SET_VECTOR_ELT(result_r, 1, thetaSamples_r);
    SET_VECTOR_ELT(resultName_r, 1, mkChar("p.theta.samples"));
    
    SET_VECTOR_ELT(result_r, 2, tauSqSamples_r);
    SET_VECTOR_ELT(resultName_r, 2, mkChar("p.tausq.samples"));
    
    SET_VECTOR_ELT(result_r, 3, rhoSamples_r);
    SET_VECTOR_ELT(resultName_r, 3, mkChar("p.rho.samples"));
    
    if (savew) {
      SET_VECTOR_ELT(result_r, 4, wSamples_r);
      SET_VECTOR_ELT(resultName_r, 4, mkChar("p.w.samples"));
    }
    
    namesgets(result_r, resultName_r);
    
    //unprotect
    UNPROTECT(nProtect);
    
    return(result_r);
  }
}