#define USE_FC_LEN_T
#include <string>
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

extern "C" {

  SEXP sSTNNGPPredict(SEXP X_r, SEXP y_r, SEXP coords_r, SEXP adjvec_r, 
                      SEXP n_r, SEXP q_r, SEXP p_r, SEXP m_r, 
                      SEXP X0_r, SEXP coords0_r, SEXP n0_r, SEXP nnIndx0_r, 
                      SEXP betaSamples_r, SEXP tauSqSamples_r, SEXP rhoSamples_r, 
                      SEXP thetaSamples_r, SEXP wSamples_r, 
                      SEXP nSamples_r, SEXP family_r, SEXP covModel_r, SEXP nThreads_r, 
                      SEXP verbose_r, SEXP nReport_r){

    int h, i, j, k, k2, l, l2, s, info, nProtect=0, MeIndx, ParentIndx;
    const int inc = 1;
    const double one = 1.0;
    const double negOne = -1.0;
    const double zero = 0.0;
    char const *lower = "L";
    
    //get args
    double *X = REAL(X_r);
    double *y = REAL(y_r);
    double *coords = REAL(coords_r);
    int *adjvec = INTEGER(adjvec_r);                                            // BJ
    int n = INTEGER(n_r)[0];
    int q = INTEGER(q_r)[0];                                                    // BJ
    int p = INTEGER(p_r)[0];
    int m = INTEGER(m_r)[0];
    int twom = 2*m;                                                             // BJ
    int twomsq = twom*twom;                                                     // BJ

    double *X0 = REAL(X0_r);
    double *coords0 = REAL(coords0_r);
    int n0 = INTEGER(n0_r)[0];                                                  // BJ

    int pq = p*q;                                                               // BJ
    int nq = n*q;                                                               // BJ
    int n0q = n0*q;                                                             // BJ
    
    int *nnIndx0 = INTEGER(nnIndx0_r);        
    double *beta = REAL(betaSamples_r);
    double *tauSq = REAL(tauSqSamples_r);                                       // BJ
    double *rho = REAL(rhoSamples_r);                                           // BJ
    double *theta = REAL(thetaSamples_r);
    double *w = REAL(wSamples_r);
    
    int nSamples = INTEGER(nSamples_r)[0];
    int family = INTEGER(family_r)[0];
    int covModel = INTEGER(covModel_r)[0];
    std::string corName = getCorName(covModel);
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
      Rprintf("\tPrediction description\n");
      Rprintf("----------------------------------------\n");
      if (family == 1) {
        Rprintf("Spanning Tree-Based NNGP Separable Latent model fit with %i observations.\n\n", n); // BJ
      } else {
        Rprintf("Spanning Tree-Based NNGP Separable Latent model fit with %i Bernoulli observations.\n\n", n); // BJ
      }
      Rprintf("Number of covariates %i (including intercept if specified).\n\n", p);
      Rprintf("Using the %s spatial correlation model.\n\n", corName.c_str());
      Rprintf("Using %i nearest neighbors.\n\n", m);
      Rprintf("Number of MCMC samples %i.\n\n", nSamples);
      Rprintf("Predicting at %i locations.\n\n", n0);  
#ifdef _OPENMP
      Rprintf("\nSource compiled with OpenMP support and model fit using %i threads.\n", nThreads);
#else
      Rprintf("\n\nSource not compiled with OpenMP support.\n");
#endif
    } 
    
    //parameters
    int nTheta, sigmaSqIndx, phiIndx, nuIndx;                                   // BJ

    if(corName != "matern"){
      nTheta = 2;//sigma^2, phi
      sigmaSqIndx = 0; phiIndx = 1;
    }else{
      nTheta = 3;//sigma^2, phi, nu
      sigmaSqIndx = 0; phiIndx = 1; nuIndx = 2;
    }
    
    //get max nu
    double nuMax = 0;
    int nb = 0;
    if(corName == "matern"){
      for(i = 0; i < nSamples; i++){
        if(theta[i*nTheta+nuIndx] > nuMax){
          nuMax = theta[i*nTheta+nuIndx];
        }
      }
      nb = 1+static_cast<int>(floor(nuMax));
    }
    double *bk = (double *) R_alloc(nThreads*nb, sizeof(double));
    
    double *C = (double *) R_alloc(nThreads*twomsq, sizeof(double)); zeros(C, nThreads*twomsq);
    double *c = (double *) R_alloc(nThreads*twom, sizeof(double)); zeros(c, nThreads*twom);
    double *tmp_B  = (double *) R_alloc(nThreads*twom, sizeof(double));         // BJ
    double phi, nu, sigmaSq, d; 
    double tauSqj, rhoj;                                                        // BJ
    int threadID = 0, status = 0;
    double psi;                                                                 // BJ

    SEXP y0_r, w0_r;
    PROTECT(y0_r = allocMatrix(REALSXP, n0q, nSamples)); nProtect++; 
    PROTECT(w0_r = allocMatrix(REALSXP, n0q, nSamples)); nProtect++;
    double *y0 = REAL(y0_r);
    double *w0 = REAL(w0_r);
 
    if(verbose){
      Rprintf("-------------------------------------------------\n");
      Rprintf("\t\tPredicting\n");
      Rprintf("-------------------------------------------------\n");
      #ifdef Win32
        R_FlushConsole();
      #endif
    }

    int zIndx = -1;
    double *wZ = (double *) R_alloc(n0q*nSamples, sizeof(double));              // BJ

    double *yZ = NULL;
    if(family == 1){
      yZ = (double *) R_alloc(n0q*nSamples, sizeof(double));                    // BJ
    }
    
    GetRNGstate();
    
    for(i = 0; i < n0q*nSamples; i++){                                          // BJ
      wZ[i] = rnorm(0.0,1.0);
    }
    
    if(family == 1){
      for(i = 0; i < n0q*nSamples; i++){                                        // BJ
        yZ[i] = rnorm(0.0,1.0);
      }
    }
    
    PutRNGstate();

    for (MeIndx = 0; MeIndx < q; MeIndx++) {                                    // BJ
      if(verbose){
        Rprintf("Variable %i\n", MeIndx+1);                                     // BJ: convert to R index
      }
      status = 0;                                                               // BJ
      if (MeIndx == 0) {                                                        // BJ: root
        for(i = 0; i < n0; i++){
#ifdef _OPENMP
#pragma omp parallel for private(threadID, phi, nu, sigmaSq, tauSqj, k, l, d, info) // BJ
#endif

          for(s = 0; s < nSamples; s++){
#ifdef _OPENMP
            threadID = omp_get_thread_num();
#endif

            phi = theta[s*nTheta+phiIndx];
            if(corName == "matern"){
              nu = theta[s*nTheta+nuIndx];
            }
            sigmaSq = theta[s*nTheta+sigmaSqIndx];
            tauSqj = tauSq[s*q+MeIndx];                                         // BJ

            for(k = 0; k < m; k++){
              d = dist2(coords[nnIndx0[i+n0*k]], coords[n+nnIndx0[i+n0*k]], coords0[i], coords0[n0+i]);
              c[threadID*twom+k] = sigmaSq*spCor(d, phi, nu, covModel, &bk[threadID*nb]);

              for(l = 0; l <= k; l++){
                d = dist2(coords[nnIndx0[i+n0*k]], coords[n+nnIndx0[i+n0*k]], coords[nnIndx0[i+n0*l]], coords[n+nnIndx0[i+n0*l]]);
                C[threadID*twomsq+l*m+k] = sigmaSq*spCor(d, phi, nu, covModel, &bk[threadID*nb]);
              }
            }

            F77_NAME(dpotrf)(lower, &m, &C[threadID*twomsq], &m, &info FCONE); if(info != 0){error("c++ error: dpotrf failed\n");}
            F77_NAME(dpotri)(lower, &m, &C[threadID*twomsq], &m, &info FCONE); if(info != 0){error("c++ error: dpotri failed\n");}
            F77_NAME(dsymv)(lower, &m, &one, &C[threadID*twomsq], &m, &c[threadID*twom], &inc, &zero, &tmp_B[threadID*twom], &inc FCONE);

            d = 0;
            for(k = 0; k < m; k++){
              d += tmp_B[threadID*twom+k]*w[s*nq+n*MeIndx+nnIndx0[i+n0*k]];     // BJ
            }

#ifdef _OPENMP
#pragma omp atomic
#endif

            zIndx++;

            w0[s*n0q+n0*MeIndx+i] = sqrt(sigmaSq - F77_NAME(ddot)(&m, &tmp_B[threadID*twom], &inc, &c[threadID*twom], &inc))*wZ[zIndx] + d; // BJ

            if(family == 1){
              y0[s*n0q+n0*MeIndx+i] = sqrt(tauSqj)*yZ[zIndx] + F77_NAME(ddot)(&p, &X0[i], &n0, &beta[s*pq+p*MeIndx], &inc) + w0[s*n0q+n0*MeIndx+i]; // BJ
            }else{//binomial
              psi = F77_NAME(ddot)(&p, &X0[i], &n0, &beta[s*pq+p*MeIndx], &inc) + w0[s*n0q+n0*MeIndx+i]; // BJ
              y0[s*n0q+n0*MeIndx+i] = 1/(1+exp(-psi));                          // BJ
            }
          }

          if(verbose){
            if(status == nReport){
              Rprintf("\tLocation %i of %i, %3.2f%%\n", i, n0, 100.0*(i+n0*MeIndx)/n0q); // BJ
#ifdef Win32
              R_FlushConsole();
#endif

              status = 0;
            }
          }
          status++;
          R_CheckUserInterrupt();
        }
      } else {
        ParentIndx = which(1, &adjvec[q*MeIndx], q);                            // BJ
        for(i = 0; i < n0; i++){
#ifdef _OPENMP
#pragma omp parallel for private(threadID, phi, nu, sigmaSq, tauSqj, rhoj, k, l, d, info) // BJ
#endif

          for(s = 0; s < nSamples; s++){
#ifdef _OPENMP
            threadID = omp_get_thread_num();
#endif

            phi = theta[s*nTheta+phiIndx];
            if(corName == "matern"){
              nu = theta[s*nTheta+nuIndx];
            }
            sigmaSq = theta[s*nTheta+sigmaSqIndx];
            tauSqj = tauSq[s*q+MeIndx];                                         // BJ
            rhoj = rho[s*(q-1)+(MeIndx-1)];                                     // BJ

            for(k = 0; k < twom; k++){
              if (k < m) {
                d = dist2(coords[nnIndx0[i+n0*k]], coords[n+nnIndx0[i+n0*k]], coords0[i], coords0[n0+i]);
                c[threadID*twom+k] = rhoj*sigmaSq*spCor(d, phi, nu, covModel, &bk[threadID*nb]);
              } else {
                k2 = k-m;
                d = dist2(coords[nnIndx0[i+n0*k2]], coords[n+nnIndx0[i+n0*k2]], coords0[i], coords0[n0+i]);
                c[threadID*twom+k] = sigmaSq*spCor(d, phi, nu, covModel, &bk[threadID*nb]);
              }

              for(l = 0; l <= k; l++){
                if (k < m) {
                  d = dist2(coords[nnIndx0[i+n0*k]], coords[n+nnIndx0[i+n0*k]], coords[nnIndx0[i+n0*l]], coords[n+nnIndx0[i+n0*l]]);
                  C[threadID*twomsq+l*twom+k] = sigmaSq*spCor(d, phi, nu, covModel, &bk[threadID*nb]);
                } else{
                  k2 = k-m;
                  if (l < m) {
                    d = dist2(coords[nnIndx0[i+n0*k2]], coords[n+nnIndx0[i+n0*k2]], coords[nnIndx0[i+n0*l]], coords[n+nnIndx0[i+n0*l]]);
                    C[threadID*twomsq+l*twom+k] = rhoj*sigmaSq*spCor(d, phi, nu, covModel, &bk[threadID*nb]);
                  } else{
                    l2 = l-m;
                    d = dist2(coords[nnIndx0[i+n0*k2]], coords[n+nnIndx0[i+n0*k2]], coords[nnIndx0[i+n0*l2]], coords[n+nnIndx0[i+n0*l2]]);
                    C[threadID*twomsq+l*twom+k] = sigmaSq*spCor(d, phi, nu, covModel, &bk[threadID*nb]);
                  }
                }
              }
            }

            F77_NAME(dpotrf)(lower, &twom, &C[threadID*twomsq], &twom, &info FCONE); if(info != 0){error("c++ error: dpotrf failed\n");}
            F77_NAME(dpotri)(lower, &twom, &C[threadID*twomsq], &twom, &info FCONE); if(info != 0){error("c++ error: dpotri failed\n");}
            F77_NAME(dsymv)(lower, &twom, &one, &C[threadID*twomsq], &twom, &c[threadID*twom], &inc, &zero, &tmp_B[threadID*twom], &inc FCONE);

            d = 0;
            for(k = 0; k < twom; k++){
              if (k < m) {
                d += tmp_B[threadID*twom+k]*w[s*nq+n*ParentIndx+nnIndx0[i+n0*k]]; // BJ
              } else {
                k2 = k-m;
                d += tmp_B[threadID*twom+k]*w[s*nq+n*MeIndx+nnIndx0[i+n0*k2]];  // BJ
              }
            }

#ifdef _OPENMP
#pragma omp atomic
#endif

            zIndx++;

            w0[s*n0q+n0*MeIndx+i] = sqrt(sigmaSq - F77_NAME(ddot)(&twom, &tmp_B[threadID*twom], &inc, &c[threadID*twom], &inc))*wZ[zIndx] + d; // BJ

            if(family == 1){
              y0[s*n0q+n0*MeIndx+i] = sqrt(tauSqj)*yZ[zIndx] + F77_NAME(ddot)(&p, &X0[i], &n0, &beta[s*pq+p*MeIndx], &inc) + w0[s*n0q+n0*MeIndx+i]; // BJ
            }else{//binomial
              psi = F77_NAME(ddot)(&p, &X0[i], &n0, &beta[s*pq+p*MeIndx], &inc) + w0[s*n0q+n0*MeIndx+i]; // BJ
              y0[s*n0q+n0*MeIndx+i] = 1/(1+exp(-psi));                          // BJ: report posterior predictive success probability
            }
          }

          if(verbose){
            if(status == nReport){
              Rprintf("\tLocation %i of %i, %3.2f%%\n", i, n0, 100.0*(i+n0*MeIndx)/n0q); // BJ
#ifdef Win32
              R_FlushConsole();
#endif

              status = 0;
            }
          }
          status++;
          R_CheckUserInterrupt();
        }
      }
      if(verbose){
        Rprintf("\tLocation %i of %i, %3.2f%%\n", i, n0, 100.0*(i+n0*MeIndx)/n0q); // BJ
#ifdef Win32
        R_FlushConsole();
#endif
      }
    }

    //make return object
    SEXP result_r, resultName_r;
    int nResultListObjs = 2;

    PROTECT(result_r = allocVector(VECSXP, nResultListObjs)); nProtect++;
    PROTECT(resultName_r = allocVector(VECSXP, nResultListObjs)); nProtect++;

    SET_VECTOR_ELT(result_r, 0, y0_r);
    SET_VECTOR_ELT(resultName_r, 0, mkChar("p.y.0")); 
    
    SET_VECTOR_ELT(result_r, 1, w0_r);
    SET_VECTOR_ELT(resultName_r, 1, mkChar("p.w.0"));

    namesgets(result_r, resultName_r);
    
    //unprotect
    UNPROTECT(nProtect);
    
    return(result_r);
  
  }
}

    
