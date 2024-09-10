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

// BJ
double func_tsq(double tsq, double crossphi, double phi, double phiParent,
                double nu, double nuParent, double crossnu) {
  return pow(pow(crossphi, 2) + tsq, 2*crossnu+2)/pow(pow(phi, 2) + tsq, nu+1)/pow(pow(phiParent, 2) + tsq, nuParent+1);
}

// BJ
double rholim(double crossphi, double phi, double phiParent,
              double nu, double nuParent, double crossnu) {
  double tsq = ((2*nu+2)*pow(phiParent,2)*pow(crossphi,2) + (2*nuParent+2)*pow(phi,2)*pow(crossphi,2) - 2*(nu+nuParent+2)*pow(phi,2)*pow(phiParent,2))/
    ((2*nu+2)*pow(phi,2) + (2*nuParent+2)*pow(phiParent,2) - 2*(nu+nuParent+2)*pow(crossphi,2));
  
  // std::cout << "tsq=0: " << func_tsq(0, crossphi, phi, phiParent, nu, nuParent, crossnu) << "\n";
  // std::cout << "tsq=tsq: " << "(" << tsq << ", " << func_tsq(tsq, crossphi, phi, phiParent, nu, nuParent, crossnu) << ")" << "\n";
  // std::cout << "tsq=Infty: " << func_tsq(1e+17, crossphi, phi, phiParent, nu, nuParent, crossnu) << "\n";
  
  double infimum = std::min(std::min(func_tsq(0, crossphi, phi, phiParent, nu, nuParent, crossnu),
                                     func_tsq(tsq, crossphi, phi, phiParent, nu, nuParent, crossnu)), // BJ: tsq can be negative. when it is, func_tsq = -nan. excluded for finding the min
                                     func_tsq(1e+17, crossphi, phi, phiParent, nu, nuParent, crossnu)); // BJ: 1e+17 role of infinity
  // std::cout << "infimum: " << infimum << "\n";
  // std::cout << "infimum: " << infimum << "\n";
  double value = gammafn(nu+1)/gammafn(nu)*gammafn(nuParent+1)/gammafn(nuParent)*
    pow(gammafn(crossnu),2)/pow(gammafn(crossnu+1),2)*
    pow(phi,2*nu)*pow(phiParent,2*nuParent)/pow(crossphi,4*crossnu)*infimum;
  
  return pow(value, 0.5);
}

// BJ
double updateBF_parent(double *B, double *F, double *c, double *C, double *coords,
                       int *nnIndx, int *nnIndxLU,
                       int *nnIndxParent, int *nnIndxLUParent, int *nnIndxLUAll,
                       int n, int m, int twomp1, int twomp1sq,
                       double sigmaSq, double phi, double nu, double tauSq, 
                       double sigmaSqParent, double phiParent, double nuParent, 
                       double tauSqParent, double rho, double crossphi, double crossnu,
                       int covModel, double nThreads, double nuUnifb, bool root) {
  
  int i, k, l;
  int info = 0;
  int inc = 1;
  double one = 1.0;
  double zero = 0.0;
  char lower = 'L';
  double logDet = 0;
  
  //bk must be 1+(int)floor(alpha) * nthread
  double *bk = (double *) R_alloc(nThreads*(1.0+static_cast<int>(floor(nuUnifb))), sizeof(double));
  int nb = 1+static_cast<int>(floor(nuUnifb));
  int threadID = 0;
  double e; 
  double crosssigma = pow(sigmaSq,0.5)*pow(sigmaSqParent,0.5);
  
  if (root) {
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
          c[twomp1*threadID+k] = sigmaSq*spCor(e, phi, nu, covModel, &bk[threadID*nb]);
          for(l = 0; l <= k; l++){
            e = dist2(coords[nnIndx[nnIndxLU[i]+k]], coords[n+nnIndx[nnIndxLU[i]+k]], coords[nnIndx[nnIndxLU[i]+l]], coords[n+nnIndx[nnIndxLU[i]+l]]);  
            C[twomp1sq*threadID+l*nnIndxLU[n+i]+k] = sigmaSq*spCor(e, phi, nu, covModel, &bk[threadID*nb]); 
            if(l == k){
              C[twomp1sq*threadID+l*nnIndxLU[n+i]+k] += tauSq;
            }
          }
        }
        F77_NAME(dpotrf)(&lower, &nnIndxLU[n+i], &C[twomp1sq*threadID], &nnIndxLU[n+i], &info FCONE);
        if(info != 0){error("c++ error: dpotrf failed\n");}
        F77_NAME(dpotri)(&lower, &nnIndxLU[n+i], &C[twomp1sq*threadID], &nnIndxLU[n+i], &info FCONE);
        F77_NAME(dsymv)(&lower, &nnIndxLU[n+i], &one, &C[twomp1sq*threadID], &nnIndxLU[n+i], &c[twomp1*threadID], &inc, &zero, &B[nnIndxLU[i]], &inc FCONE);
        F[i] = sigmaSq - F77_NAME(ddot)(&nnIndxLU[n+i], &B[nnIndxLU[i]], &inc, &c[twomp1*threadID], &inc) + tauSq;                     // BJ
      } else {
        B[i] = 0;
        F[i] = sigmaSq + tauSq;
      }
    }
  } else {
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
          c[twomp1*threadID+k] = rho*crosssigma*
            spCor(e, crossphi, crossnu, covModel, &bk[threadID*nb]);
        } else {
          int k2 = k-m-1;
          e = dist2(coords[i], coords[n+i],
                    coords[nnIndx[nnIndxLU[i]+k2]], coords[n+nnIndx[nnIndxLU[i]+k2]]);
          c[twomp1*threadID+k] = sigmaSq*
            spCor(e, phi, nu, covModel, &bk[threadID*nb]);
        }
        
        for(l = 0; l <= k; l++){
          if (k < m+1) {
            e = dist2(coords[nnIndxParent[nnIndxLUParent[i]+k]],coords[n+nnIndxParent[nnIndxLUParent[i]+k]],
                      coords[nnIndxParent[nnIndxLUParent[i]+l]],coords[n+nnIndxParent[nnIndxLUParent[i]+l]]);
            C[twomp1sq*threadID+l*nnIndxLUAll[n+i]+k] = sigmaSqParent*
              spCor(e, phiParent, nuParent, covModel, &bk[threadID*nb]);
            if (l == k) {
              C[twomp1sq*threadID+l*nnIndxLUAll[n+i]+k] += tauSqParent;
            }
          } else {
            int k2 = k-m-1;
            if (l < m+1) {
              e = dist2(coords[nnIndx[nnIndxLU[i]+k2]],coords[n+nnIndx[nnIndxLU[i]+k2]],
                        coords[nnIndxParent[nnIndxLUParent[i]+l]],coords[n+nnIndxParent[nnIndxLUParent[i]+l]]);
              C[twomp1sq*threadID+l*nnIndxLUAll[n+i]+k] = rho*crosssigma*
                spCor(e, crossphi, crossnu, covModel, &bk[threadID*nb]);
            } else {
              int l2 = l-m-1;
              e = dist2(coords[nnIndx[nnIndxLU[i]+k2]],coords[n+nnIndx[nnIndxLU[i]+k2]],
                        coords[nnIndx[nnIndxLU[i]+l2]],coords[n+nnIndx[nnIndxLU[i]+l2]]);
              C[twomp1sq*threadID+l*nnIndxLUAll[n+i]+k] = sigmaSq*
                spCor(e, phi, nu, covModel, &bk[threadID*nb]);
              if (l == k) {
                C[twomp1sq*threadID+l*nnIndxLUAll[n+i]+k] += tauSq;
              }
            }
          }
        }
      }
      F77_NAME(dpotrf)(&lower, &nnIndxLUAll[n+i], &C[twomp1sq*threadID], &nnIndxLUAll[n+i], &info FCONE);
      if(info != 0){error("c++ error: dpotrf failed\n");}
      F77_NAME(dpotri)(&lower, &nnIndxLUAll[n+i], &C[twomp1sq*threadID], &nnIndxLUAll[n+i], &info FCONE);
      F77_NAME(dsymv)(&lower, &nnIndxLUAll[n+i], &one, &C[twomp1sq*threadID], &nnIndxLUAll[n+i], &c[twomp1*threadID], &inc, &zero, &B[nnIndxLUAll[i]], &inc FCONE);
      F[i] = sigmaSq - F77_NAME(ddot)(&nnIndxLUAll[n+i], &B[nnIndxLUAll[i]], &inc, &c[twomp1*threadID], &inc) + tauSq;                     // BJ
    }
  }
  
  for(i = 0; i < n; i++){
    logDet += log(F[i]);
  }
  
  return(logDet);
}

extern "C" {
  
  SEXP rSTNNGP_NS(SEXP y_r, SEXP X_r, 
               SEXP q_r,                                                     // BJ
               SEXP p_r, SEXP n_r, SEXP m_r, SEXP coords_r, 
               SEXP covModel_r, SEXP nnIndx_r, SEXP nnIndxLU_r, 
               SEXP nnIndxParent_r, SEXP nnIndxLUParent_r, SEXP nnIndxLUAll_r, // BJ
               SEXP sigmaSqIGa_r, SEXP sigmaSqIGb_r,                         // BJ 
               SEXP tauSqIGa_r, SEXP tauSqIGb_r,                             // BJ
               SEXP phiUnifa_r, SEXP phiUnifb_r,                             // BJ
               SEXP nuUnifa_r, SEXP nuUnifb_r,                               // BJ
               SEXP betaStarting_r, SEXP sigmaSqStarting_r, SEXP tauSqStarting_r, 
               SEXP phiStarting_r, 
               SEXP crossphiStarting_r,                                      // BJ
               SEXP nuStarting_r, 
               SEXP rhoStarting_r, SEXP adjmatStarting_r,                    // BJ
               SEXP sigmaSqTuning_r, SEXP tauSqTuning_r, SEXP phiTuning_r, 
               SEXP crossphiTuning_r,                                        // BJ
               SEXP nuTuning_r, 
               SEXP rhoTuning_r,                                             // BJ
               SEXP nSamples_r, SEXP nThreads_r, SEXP verbose_r, SEXP nReport_r){
    
    int h, i, j, k, l, s, info, nProtect=0;
    const int inc = 1;
    const double one = 1.0;
    const double negOne = -1.0;
    const double zero = 0.0;
    char const *lower = "L";
    char const *ntran = "N";
    
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
    int covModel = INTEGER(covModel_r)[0];
    std::string corName = getCorName(covModel);
    
    //priors
    double *sigmaSqIGa = REAL(sigmaSqIGa_r);                                    // BJ
    double *sigmaSqIGb = REAL(sigmaSqIGb_r);                                    // BJ
    double *phiUnifa = REAL(phiUnifa_r);                                        // BJ
    double *phiUnifb = REAL(phiUnifb_r);                                        // BJ                                                  // BJ
    double *nuUnifa = (double *) R_alloc(q, sizeof(double));                    // BJ              
    double *nuUnifb = (double *) R_alloc(q, sizeof(double));                    // BJ
    if(corName == "matern"){                                                    // BJ
      F77_NAME(dcopy)(&q, REAL(nuUnifa_r), &inc, nuUnifa, &inc);                // BJ
      F77_NAME(dcopy)(&q, REAL(nuUnifb_r), &inc, nuUnifb, &inc);                // BJ
    }                                                                         // BJ
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
      Rprintf("Spanning Tree-Based NNGP Response model fit with %i observations.\n\n", n); // BJ
      Rprintf("Number of covariates %i (including intercept if specified).\n\n", p);
      Rprintf("Using the %s spatial correlation model.\n\n", corName.c_str());
      Rprintf("Using %i nearest neighbors.\n\n", m);
      Rprintf("Number of MCMC samples %i.\n\n", nSamples);
      Rprintf("Priors and hyperpriors:\n");
      Rprintf("\tbeta flat\n");
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
    int pq = p*q;                                                               // BJ
    int qm1 = q-1;                                                              // BJ
    double *beta = (double *) R_alloc(pq, sizeof(double));                      // BJ
    double *betaj = (double *) R_alloc(p, sizeof(double));                      // BJ
    double *beta2 = (double *) R_alloc(pq, sizeof(double));                     // BJ
    double *beta2j = (double *) R_alloc(p, sizeof(double));                     // BJ
    double *sigmaSq = (double *) R_alloc(q, sizeof(double));                    // BJ
    double *phi = (double *) R_alloc(q, sizeof(double));                        // BJ
    double *crossphi = (double *) R_alloc(qm1, sizeof(double));                 // BJ
    double *nu = (double *) R_alloc(q, sizeof(double));                         // BJ
    double *tauSq = (double *) R_alloc(q, sizeof(double));                      // BJ
    double *rho = (double *) R_alloc(qm1, sizeof(double));                      // BJ
    int *adjvec = INTEGER(adjmatStarting_r);                                    // BJ
    
    F77_NAME(dcopy)(&pq, REAL(betaStarting_r), &inc, beta, &inc);               // BJ
    F77_NAME(dcopy)(&pq, REAL(betaStarting_r), &inc, beta2, &inc);              // BJ
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
    double *tauSqTuning = REAL(tauSqTuning_r);                                  // BJ
    double *rhoTuning = REAL(rhoTuning_r);                                      // BJ
    double *crossphiTuning = REAL(crossphiTuning_r);                            // BJ
    
    //return stuff  
    SEXP betaSamples_r, sigmaSqSamples_r, phiSamples_r, crossphiSamples_r,      // BJ
    tauSqSamples_r, rhoSamples_r, nuSamples_r;                                  // BJ
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
    int nIndx = static_cast<int>(static_cast<double>(m*(m-1)/2+(n-m)*m+(m+1)*n));// BJ
    int twomp1 = 2*m+1;                                                         // BJ
    int twomp1sq = pow(twomp1,2);                                               // BJ
    double *sigmaSqCand = (double *) R_alloc(q, sizeof(double));                // BJ
    double *phiCand = (double *) R_alloc(q, sizeof(double));                    // BJ
    double *crossphiCand = (double *) R_alloc(qm1, sizeof(double));             // BJ
    F77_NAME(dcopy)(&q, crossphi, &inc, crossphiCand, &inc);                    // BJ
    double *nuCand = (double *) R_alloc(q, sizeof(double));                     // BJ
    F77_NAME(dcopy)(&q, nu, &inc, nuCand, &inc);                                // BJ
    double *tauSqCand = (double *) R_alloc(q, sizeof(double));                  // BJ
    F77_NAME(dcopy)(&q, tauSq, &inc, tauSqCand, &inc);                          // BJ
    double *rhoCand = (double *) R_alloc(qm1, sizeof(double));                  // BJ
    
    double crossphiCandUnifa, crossphiCandUnifb, rhoCandUnifa = -1, rhoCandUnifb = 1;    // BJ
    
    double *B = (double *) R_alloc(nIndx, sizeof(double));                      // BJ
    double *F = (double *) R_alloc(n, sizeof(double));                          // BJ
    double *c =(double *) R_alloc(twomp1*nThreads, sizeof(double));             // BJ
    double *C = (double *) R_alloc(twomp1sq*nThreads, sizeof(double));          // BJ
    
    double logPostCand, logPostCurrent, logDetCurrent, logDetCand, QCurrent, QCand;     
    double QCurrent2, logPriorJacobianCurrent, logPriorJacobianCand;            // BJ
    int accept = 0, batchAccept = 0, status = 0;
    int pp = p*p;
    double *tmp_pp = (double *) R_alloc(pp, sizeof(double));
    double *tmp_p = (double *) R_alloc(p, sizeof(double));
    double *tmp_p2 = (double *) R_alloc(p, sizeof(double));
    double *tmp_n = (double *) R_alloc(n, sizeof(double));
    double *tmp_n_parent = (double *) R_alloc(n, sizeof(double));               // BJ
    double *tmp_zero = (double *) R_alloc(n, sizeof(double));                   // BJ
    zeros(tmp_zero, n);                                                         // BJ
    
    bool thetaUpdate = true;
    bool root;                                                                  // BJ
    int ParentIndx;                                                             // BJ
    double crossnu, crossnuCand, infimum;                                       // BJ
    
    if(verbose){
      Rprintf("----------------------------------------\n");
      Rprintf("\t\tSampling\n");
      Rprintf("----------------------------------------\n");
#ifdef Win32
      R_FlushConsole();
#endif
    }
    
    GetRNGstate();
    
    for(s = 0; s < nSamples; s++){
      if(thetaUpdate){                                                          // BJ
        thetaUpdate = false;

        root = true;                                                            // BJ
        ParentIndx = 0;                                                         // BJ: for MeIndx = 0, does not affect results
        logDetCurrent = 0;
        QCurrent = 0;
        QCurrent2 = 0;                                                          // BJ
        logPriorJacobianCurrent = 0;                                            // BJ

        for (int MeIndx = 0; MeIndx < q; MeIndx++) {                            // BJ
          if (MeIndx > 0) {                                                     // BJ
            root = false;                                                       // BJ
            ParentIndx = which(1, &adjvec[q*MeIndx], q);                        // BJ

            crossphiUnifa = std::min(phi[MeIndx], phi[ParentIndx]);             // BJ
            crossphiUnifb = std::max(phi[MeIndx], phi[ParentIndx]);             // BJ
            logPriorJacobianCurrent += log(crossphi[MeIndx-1] - crossphiUnifa) +// BJ
              log(crossphiUnifb - crossphi[MeIndx-1]);                          // BJ
            
            crossnu = 0.5*nu[MeIndx] + 0.5*nu[ParentIndx];                      // BJ
            rhoUnifb = rholim(crossphi[MeIndx-1], phi[MeIndx], phi[ParentIndx],
                              nu[MeIndx], nu[ParentIndx], crossnu);             // BJ
            rhoUnifa = -1.0*rhoUnifb;                                           // BJ
            // std::cout << "(rhoUnifa, rhoUnifb): " << rhoUnifa << ", " << rhoUnifb << "\n";
            logPriorJacobianCurrent += log(rho[MeIndx-1] - rhoUnifa) +          // BJ
              log(rhoUnifb - rho[MeIndx-1]);                                    // BJ
            
          }                                                                     // BJ

          logPriorJacobianCurrent += log(phi[MeIndx] - phiUnifa[MeIndx]) + log(phiUnifb[MeIndx] - phi[MeIndx]);
          logPriorJacobianCurrent += -1.0*(1.0+sigmaSqIGa[MeIndx])*log(sigmaSq[MeIndx])-sigmaSqIGb[MeIndx]/sigmaSq[MeIndx]+log(sigmaSq[MeIndx]);
          if (tauSq[MeIndx] != 0) {
            logPriorJacobianCurrent += -1.0*(1.0+tauSqIGa[MeIndx])*log(tauSq[MeIndx])-tauSqIGb[MeIndx]/tauSq[MeIndx]+log(tauSq[MeIndx]);
          }
          if (corName == "matern"){
            logPriorJacobianCurrent += log(nu[MeIndx] - nuUnifa[MeIndx]) + log(nuUnifb[MeIndx] - nu[MeIndx]);
          }

          ///////////////////
          // BJ: update beta
          //////////////////
          logDetCurrent += updateBF_parent(B, F, c, C, coords,
                                           nnIndx, nnIndxLU,
                                           nnIndxParent, nnIndxLUParent, nnIndxLUAll,
                                           n, m, twomp1, twomp1sq,
                                           sigmaSq[MeIndx], phi[MeIndx],
                                           nu[MeIndx], tauSq[MeIndx],
                                           sigmaSq[ParentIndx], phi[ParentIndx],
                                           nu[ParentIndx], tauSq[ParentIndx],
                                           rho[MeIndx-1], crossphi[MeIndx-1], crossnu,
                                           covModel, nThreads, nuUnifb[MeIndx], root);

          for(i = 0; i < p; i++){
            tmp_p[i] = Q_parent(B, F, &X[n*i], &y[n*MeIndx], tmp_zero, tmp_zero,
                                n, nnIndx, nnIndxLU, nnIndxParent, nnIndxLUParent,
                                nnIndxLUAll, root);
            for(j = 0; j <= i; j++){
              tmp_pp[j*p+i] = Q_parent(B, F, &X[n*j], &X[n*i], tmp_zero, tmp_zero,
                                       n, nnIndx, nnIndxLU, nnIndxParent, nnIndxLUParent,
                                       nnIndxLUAll, root);
            }
          }

          F77_NAME(dpotrf)(lower, &p, tmp_pp, &p, &info FCONE);
          if(info != 0){
            error("c++ error: dpotrf failed\n");
          }
          F77_NAME(dpotri)(lower, &p, tmp_pp, &p, &info FCONE);
          if(info != 0){
            error("c++ error: dpotri failed\n");
          }
          F77_NAME(dsymv)(lower, &p, &one, tmp_pp, &p, tmp_p, &inc, &zero, tmp_p2, &inc FCONE);
          F77_NAME(dpotrf)(lower, &p, tmp_pp, &p, &info FCONE);
          if(info != 0){
            error("c++ error: dpotrf failed\n");
          }
          // std::cout << "beta posterior mean " << tmp_p2[0] << ", " << tmp_p2[1] << "\n";
          // std::cout << "beta posterior var " << tmp_pp[0] << ", " << tmp_pp[1] << ", " << tmp_pp[2] << ", " << tmp_pp[3] << "\n";

          mvrnorm(betaj, tmp_p2, tmp_pp, p); 
          F77_NAME(dcopy)(&p, betaj, &inc, &beta[p*MeIndx], &inc);
          // std::cout << "beta " << betaj[0] << ", " << betaj[1] << "\n";

          F77_NAME(dgemv)(ntran, &n, &p, &one, X, &n, betaj, &inc, &zero, tmp_n, &inc FCONE);
          F77_NAME(daxpy)(&n, &negOne, &y[n*MeIndx], &inc, tmp_n, &inc);
          F77_NAME(dgemv)(ntran, &n, &p, &one, X, &n, &beta[p*ParentIndx], &inc, &zero, tmp_n_parent, &inc FCONE);
          F77_NAME(daxpy)(&n, &negOne, &y[n*ParentIndx], &inc, tmp_n_parent, &inc);
          QCurrent += Q_parent(B, F, tmp_n, tmp_n, tmp_n_parent, tmp_n_parent,
                               n, nnIndx, nnIndxLU,
                               nnIndxParent, nnIndxLUParent, nnIndxLUAll, root);

          // BJ: backup in case theta is not updated
          mvrnorm(beta2j, tmp_p2, tmp_pp, p);
          F77_NAME(dcopy)(&p, beta2j, &inc, &beta2[p*MeIndx], &inc);

          F77_NAME(dgemv)(ntran, &n, &p, &one, X, &n, beta2j, &inc, &zero, tmp_n, &inc FCONE);
          F77_NAME(daxpy)(&n, &negOne, &y[n*MeIndx], &inc, tmp_n, &inc);
          F77_NAME(dgemv)(ntran, &n, &p, &one, X, &n, &beta2[p*ParentIndx], &inc, &zero, tmp_n_parent, &inc FCONE);
          F77_NAME(daxpy)(&n, &negOne, &y[n*ParentIndx], &inc, tmp_n_parent, &inc);
          QCurrent2 += Q_parent(B, F, tmp_n, tmp_n, tmp_n_parent, tmp_n_parent,
                                n, nnIndx, nnIndxLU,
                                nnIndxParent, nnIndxLUParent, nnIndxLUAll, root);
        }
      } else {
        QCurrent = QCurrent2;                                                   // BJ
        F77_NAME(dcopy)(&pq, beta2, &inc, beta, &inc);                          // BJ
      }
      // Rprintf("At %ith iteration,\n \tlogDet=%f\n", s, logDetCurrent);            // BJ
      // Rprintf("\tQ=%f\n", QCurrent);                                            // BJ
      // Rprintf("\tbeta[0]=%f\n", beta[0]);                                       // BJ

      /////////////////////
      // BJ: update theta
      ////////////////////
      logPostCurrent = -0.5*logDetCurrent - 0.5*QCurrent + logPriorJacobianCurrent; // BJ                                               // BJ
      // Rprintf("\tlogPostCurrent=%f\n", logPostCurrent);                         // BJ

      //update logDet and Q
      root = true;                                                              // BJ
      ParentIndx = 0;                                                           // BJ: for MeIndx = 0, does not affect results
      logDetCand = 0;
      QCand = 0;
      logPriorJacobianCand = 0;                                                 // BJ

      for (int MeIndx = 0; MeIndx < q; MeIndx++) {                              // BJ
        phiCand[MeIndx] = logitInv(rnorm(logit(phi[MeIndx], phiUnifa[MeIndx], phiUnifb[MeIndx]),
                                         phiTuning[MeIndx]), phiUnifa[MeIndx], phiUnifb[MeIndx]); // BJ
        sigmaSqCand[MeIndx] = exp(rnorm(log(sigmaSq[MeIndx]), sigmaSqTuning[MeIndx])); // BJ
        if (corName == "matern"){
          nuCand[MeIndx] = logitInv(rnorm(logit(nu[MeIndx], nuUnifa[MeIndx], nuUnifb[MeIndx]),
                                          nuTuning[MeIndx]), nuUnifa[MeIndx], nuUnifb[MeIndx]); // BJ
        }
        if (tauSqTuning[MeIndx] != 0) {                                         // BJ
          tauSqCand[MeIndx] = exp(rnorm(log(tauSq[MeIndx]), tauSqTuning[MeIndx])); // BJ
        }                                                                       // BJ
        if (MeIndx > 0) {                                                       // BJ
          root = false;                                                         // BJ
          ParentIndx = which(1, &adjvec[q*MeIndx], q);                          // BJ

          crossphiUnifa = std::min(phi[MeIndx], phi[ParentIndx]);               // BJ
          crossphiUnifb = std::max(phi[MeIndx], phi[ParentIndx]);               // BJ
          crossphiCandUnifa = std::min(phiCand[MeIndx], phiCand[ParentIndx]);   // BJ
          crossphiCandUnifb = std::max(phiCand[MeIndx], phiCand[ParentIndx]);   // BJ          
          
          crossphiCand[MeIndx-1] = logitInv(rnorm(logit(crossphi[MeIndx-1], crossphiUnifa, crossphiUnifb),
                                                  crossphiTuning[MeIndx-1]), crossphiCandUnifa, crossphiCandUnifb); // BJ
          
          crossnu = 0.5*nu[MeIndx] + 0.5*nu[ParentIndx];                      // BJ
          rhoUnifb = rholim(crossphi[MeIndx-1], phi[MeIndx], phi[ParentIndx],
                            nu[MeIndx], nu[ParentIndx], crossnu);             // BJ
          rhoUnifa = -1.0*rhoUnifb;                                           // BJ
          // std::cout << "(rhoUnifa, rhoUnifb): " << rhoUnifa << ", " << rhoUnifb << "\n";
          
          crossnuCand = 0.5*nuCand[MeIndx] + 0.5*nuCand[ParentIndx];            // BJ
          // std::cout << "crossphiCand[MeIndx-1]: " << crossphiCand[MeIndx-1] << "\n";
          // std::cout << "phiCand[MeIndx]: " << phiCand[MeIndx] << "\n";
          // std::cout << "phiCand[ParentIndx]: " << phiCand[ParentIndx] << "\n";
          // std::cout << "nuCand[MeIndx]: " << nuCand[MeIndx] << "\n";
          // std::cout << "nuCand[ParentIndx]: " << nuCand[ParentIndx] << "\n";
          // std::cout << "crossnuCand: " << crossnuCand << "\n";
          rhoCandUnifb = rholim(crossphiCand[MeIndx-1], phiCand[MeIndx], 
                                phiCand[ParentIndx], nuCand[MeIndx], 
                                nuCand[ParentIndx], crossnuCand);               // BJ
          rhoCandUnifa = -1.0*rhoCandUnifb;                                     // BJ
          // std::cout << "(rhoCandUnifa, rhoCandUnifb): " << rhoCandUnifa << ", " << rhoCandUnifb << "\n";

          rhoCand[MeIndx-1] = logitInv(rnorm(logit(rho[MeIndx-1], rhoUnifa, rhoUnifb), rhoTuning[MeIndx-1]), rhoCandUnifa, rhoCandUnifb);

          logPriorJacobianCand += log(crossphiCand[MeIndx-1] - crossphiCandUnifa) +// BJ
            log(crossphiCandUnifb - crossphiCand[MeIndx-1]);                    // BJ
          
          logPriorJacobianCand += log(rhoCand[MeIndx-1] - rhoCandUnifa) +       // BJ
            log(rhoCandUnifb - rhoCand[MeIndx-1]);                              // BJ
        }                                                                       // BJ

        logPriorJacobianCand += log(phiCand[MeIndx] - phiUnifa[MeIndx]) +
          log(phiUnifb[MeIndx] - phiCand[MeIndx]);
        logPriorJacobianCand += -1.0*(1.0+sigmaSqIGa[MeIndx])*log(sigmaSqCand[MeIndx])-
          sigmaSqIGb[MeIndx]/sigmaSqCand[MeIndx]+log(sigmaSqCand[MeIndx]);
        if (tauSqCand[MeIndx] != 0) {
          logPriorJacobianCand += -1.0*(1.0+tauSqIGa[MeIndx])*log(tauSqCand[MeIndx])-
            tauSqIGb[MeIndx]/tauSqCand[MeIndx]+log(tauSqCand[MeIndx]);
        }
        if (corName == "matern"){
          logPriorJacobianCand += log(nuCand[MeIndx] - nuUnifa[MeIndx]) +
            log(nuUnifb[MeIndx] - nuCand[MeIndx]);
        }

        logDetCand += updateBF_parent(B, F, c, C, coords,
                                      nnIndx, nnIndxLU,
                                      nnIndxParent, nnIndxLUParent, nnIndxLUAll,
                                      n, m, twomp1, twomp1sq,
                                      sigmaSqCand[MeIndx], phiCand[MeIndx],
                                      nuCand[MeIndx], tauSqCand[MeIndx],
                                      sigmaSqCand[ParentIndx], phiCand[ParentIndx],
                                      nuCand[ParentIndx], tauSqCand[ParentIndx],
                                      rhoCand[MeIndx-1], crossphiCand[MeIndx-1], 
                                      crossnuCand,
                                      covModel, nThreads, nuUnifb[MeIndx], root);

        F77_NAME(dgemv)(ntran, &n, &p, &one, X, &n, &beta[p*MeIndx], &inc, &zero, tmp_n, &inc FCONE);
        F77_NAME(daxpy)(&n, &negOne, &y[n*MeIndx], &inc, tmp_n, &inc);
        F77_NAME(dgemv)(ntran, &n, &p, &one, X, &n, &beta[p*ParentIndx], &inc, &zero, tmp_n_parent, &inc FCONE);
        F77_NAME(daxpy)(&n, &negOne, &y[n*ParentIndx], &inc, tmp_n_parent, &inc);
        QCand += Q_parent(B, F, tmp_n, tmp_n, tmp_n_parent, tmp_n_parent,
                          n, nnIndx, nnIndxLU,
                          nnIndxParent, nnIndxLUParent, nnIndxLUAll, root);
      }

      logPostCand = -0.5*logDetCand - 0.5*QCand + logPriorJacobianCand;         // BJ
      // Rprintf("\tlogDetCand=%f\n", logDetCand);                                 // BJ
      // Rprintf("\tQCand=%f\n", QCand);                                           // BJ
      // Rprintf("\tlogPostCand=%f\n", logPostCand);                               // BJ

      if (runif(0.0, 1.0) <= exp(logPostCand - logPostCurrent)) {
        thetaUpdate = true;
        dcopy_(&q, sigmaSqCand, &inc, sigmaSq, &inc);
        dcopy_(&q, phiCand, &inc, phi, &inc);
        dcopy_(&qm1, crossphiCand, &inc, crossphi, &inc);
        if (corName == "matern") {
          dcopy_(&q, nuCand, &inc, nu, &inc);
        }
        dcopy_(&qm1, rhoCand, &inc, rho, &inc);
        dcopy_(&q, tauSqCand, &inc, tauSq, &inc);
        accept++;
        batchAccept++;
      }
      //save samples
      F77_NAME(dcopy)(&pq, beta, &inc, &REAL(betaSamples_r)[s*pq], &inc);
      F77_NAME(dcopy)(&q, sigmaSq, &inc, &REAL(sigmaSqSamples_r)[s*q], &inc);
      F77_NAME(dcopy)(&q, phi, &inc, &REAL(phiSamples_r)[s*q], &inc);
      F77_NAME(dcopy)(&qm1, crossphi, &inc, &REAL(crossphiSamples_r)[s*qm1], &inc);
      if (corName == "matern") {
        F77_NAME(dcopy)(&q, nu, &inc, &REAL(nuSamples_r)[s*q], &inc);
      }
      F77_NAME(dcopy)(&qm1, rho, &inc, &REAL(rhoSamples_r)[s*qm1], &inc);
      F77_NAME(dcopy)(&q, tauSq, &inc, &REAL(tauSqSamples_r)[s*q], &inc);

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
    int nResultListObjs = 6;
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
    
    if (corName == "matern") {
      SET_VECTOR_ELT(result_r, 6, nuSamples_r);
      SET_VECTOR_ELT(resultName_r, 6, mkChar("p.nu.samples"));
    }
    namesgets(result_r, resultName_r);
    
    //unprotect
    UNPROTECT(nProtect);
    
    return(result_r);
  }
}