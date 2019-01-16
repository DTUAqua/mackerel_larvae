#include <TMB.hpp>

template<class Type>
Type objective_function<Type>::operator() ()
{
  using namespace Eigen;
  using namespace density;

  DATA_VECTOR        (N);        /* Observations (counts) */
  DATA_FACTOR        (position); /* Observations */
  DATA_FACTOR        (time);     /* Observations */

  DATA_SPARSE_MATRIX (Q0);       /* For random field */
  DATA_SPARSE_MATRIX (I);        /* For random field */
  DATA_SPARSE_MATRIX (A);        /* Design matrix (standardization) */

  /* Fixed effects */
  PARAMETER          (logdelta); /* For random field (corr) */
  PARAMETER          (logkappa); /* For random field (scale) */
  PARAMETER          (logkappa_static); /* For random field (scale) */
  PARAMETER          (tphi_time);/* One-step time correlation */
  PARAMETER          (logsigma); /* Nugget error */
  PARAMETER_VECTOR   (beta);     /* For design matrix */

  /* Random effects */
  PARAMETER_ARRAY(eta);         // 2D: space x time
  PARAMETER_ARRAY(etanug);      // 1D: haul
  PARAMETER_VECTOR(eta_static); // 1D: space

  /* Parameter transforms */
  Type delta = exp(logdelta);
  Type kappa = exp(logkappa);
  Type kappa_static = exp(logkappa_static);
  Type sigma = exp(logsigma);
  Type phi_time = tphi_time / sqrt( 1.0 + tphi_time * tphi_time );

  /* Random fields */
  SparseMatrix<Type>       Q = (Q0 + delta * I);
  GMRF_t<Type>             nldens = GMRF(Q);
  SCALE_t<GMRF_t<Type> >   nldens_spatial = SCALE(nldens, kappa);
  AR1_t<N01<Type> >        nldens_time = AR1(phi_time);
  SCALE_t<GMRF_t<Type> >   nldens_spatial_static = SCALE(nldens, kappa_static);

  Type nll = 0; // Negative log likelhood

  // For report
  DATA_SCALAR(h);

  // Robustify
  DATA_SCALAR(huge);
  if (huge > 0) nll -= dnorm(beta, Type(0), huge, true).sum();

  // Process likelihood
  nll += SEPARABLE(nldens_time, nldens_spatial)(eta);

  // Static field
  nll += nldens_spatial_static(eta_static);

  // Nugget
  nll -= dnorm(vector<Type>(etanug), Type(0), sigma, true).sum();

  // Measurement likelihood
  vector<Type> predictor = A * beta;
  for(int i=0; i<N.size(); i++){
    Type loglambda =
      predictor(i) +
      eta_static(position(i)) +
      eta(position(i), time(i)) +
      etanug(i);
    // Robust dpois evaluation:
    nll += -N(i)*loglambda + exp(loglambda);
  }
  nll += lfactorial(N).sum();

  // REPORT
  vector<Type> sqrtindex(NLEVELS(time));
  sqrtindex.setZero();
  for(int j=0; j<eta.cols(); j++) {
    for(int i=0; i<eta.rows(); i++) {
      sqrtindex(j) += sqrt( exp( eta(i, j) + eta_static(i) + beta[j] ) );
    }
    sqrtindex(j) /= eta.rows();
  }
  ADREPORT(sqrtindex);
  REPORT(sqrtindex);

  Type T = -1./log(phi_time);
  ADREPORT(T);
  Type H = h / log(1 + (delta/2) + sqrt(delta + (delta*delta/4)));
  ADREPORT(H);

  // if(isDouble<Type>::value) {
  //   Type avgVariance = nldens_spatial.variance().mean();
  //   Type SpatialSD = sqrt(avgVariance);
  //   Type ForecastSD = sqrt(1. - phi*phi) * SpatialSD;
  //   REPORT(SpatialSD);
  //   REPORT(ForecastSD);
  //   Type sigma_square = nldens.variance().mean();
  // }

  return nll;

}

