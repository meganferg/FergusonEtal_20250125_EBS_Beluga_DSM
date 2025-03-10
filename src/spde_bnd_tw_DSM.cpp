//Script spde_bnd_tw_DSM.cpp...Megan C. Ferguson...31 January 2025

  //NOTES
  //
  //1. This script is identical to NSDL17dsm_tmb_spde_tw_bnd.cpp in MCF's 
  //original files.


#include <TMB.hpp>
using namespace Eigen;
using namespace tmbutils;

#include "barrierTMB.hpp"

// get sign of double, only for REPORT use
// Function copied from Thorson's tinyVAST.cpp.
template<class Type>
Type sign(Type x){
  return x / pow(pow(x,2),0.5);
}

// Deviance for the Tweedie
// https://en.wikipedia.org/wiki/Tweedie_distribution#Properties
// Function copied from Thorson's tinyVAST.cpp. 
//   y: data point. must be a scalar.
//   mu: Tweedie mean. must be a scalar.
//   p: Tweedie power parameter. must be a scalar. 
template<class Type>
Type devresid_tweedie( Type y,
                       Type mu,
                       Type p ){
  
  Type c1 = pow( y, 2.0-p ) / (1.0-p) / (2.0-p);
  Type c2 = y * pow( mu, 1.0-p ) / (1.0-p);
  Type c3 = pow( mu, 2.0-p ) / (2.0-p);
  Type deviance = 2 * (c1 - c2 + c3 );
  Type devresid = sign( y - mu ) * pow( deviance, 0.5 );
  return devresid;
  
}

template<class Type>
  Type objective_function<Type>::operator() ()
{
  using namespace R_inla; //includes SPDE-specific functions, e.g. Q_spde()
  using namespace Eigen;  //Needed for utilisation of sparse structures
  using namespace density;

  //Load data and parameters
  DATA_VECTOR(y);      //The response
  DATA_SPARSE_MATRIX(A);  //Matrix for interpolating observation points within triangles
  DATA_SPARSE_MATRIX(A_pred); //Matrix for interpolating prediction points within triangles
  DATA_STRUCT(spdeMatrices,spde_t); //Three matrices needed for representing the GMRF, see p. 8 in Lindgren et al. (2011)
  DATA_STRUCT(spdeMatricesBarrier,spde_barrier_t); //Structure needed for the barrier procedure
  DATA_VECTOR(c);      // Scaling of range
  DATA_INTEGER(barrier); // if barrier==1, barrier procedure is used
  DATA_MATRIX(X);         //Design matrix for fixed effects for observations
  DATA_MATRIX(X_pred);         //Design matrix for fixed effects for prediction
  DATA_VECTOR(offset);  //offset term for model-building
  DATA_VECTOR(offset_pred);  //offset term for predictions

  PARAMETER_VECTOR(beta); //Regression coefficients for intercept and any covariates
  PARAMETER(log_tau); //log precision for spatial effect
  PARAMETER(log_kappa); //log spatial scale parameter in spatial effect
  PARAMETER_VECTOR(x); //Spatial effect
  PARAMETER(log_phi);     // log(error SD)
  PARAMETER(finv_power);     // transformed power

  //Declare transformed parameters and derived variables
  Type tau = exp(log_tau);
  Type kappa = exp(log_kappa);
  Type phi = exp(log_phi);
  Type power = 1.0 + (exp(finv_power) / (1.0 + exp(finv_power)));
  Type Nhat = 0; // abundance derived from sampled segments in DATA
  vector<Type> devresid( y.size() ); // deviance residual vector
  
  // Spatial interpolation
  vector<Type> delta = (A*x)/tau;  //observations
  vector<Type> delta_pred = (A_pred*x)/tau;
  
  //Construct sparse precision matrix for latent field
  SparseMatrix<Type> Q;
  if(barrier==0){
    Q = Q_spde(spdeMatrices,kappa);
  }else if(barrier ==1){
    Q = Q_spde(spdeMatricesBarrier,kappa,c);
  }

  //Calculates nll
  Type nll = 0.0;
  nll += GMRF(Q)(x);

  vector<Type> mu(y.size());
  mu = X*beta + delta; // linear predictor for sampled segments;
                       // without offset; on link scale

  for(int i=0; i<y.size(); i++){
    nll -= dtweedie( y(i), offset(i)*exp(mu(i)), phi, power, true );
    Nhat += offset(i)*exp(mu(i));
    devresid(i) = devresid_tweedie( y(i), offset(i)*exp(mu(i)), power ); // deviance residual
  }
  
  // Prediction to grid. 
    vector<Type> mu_pred = X_pred*beta + delta_pred; // predicted density on log scale
    
    Type Nhat_pred = 0; // abundance derived from density surface
    for(int i=0; i<offset_pred.size(); i++){
      Nhat_pred += offset_pred(i)*exp(mu_pred(i));
    }

  //Report what we want to report
  Type range = sqrt(8)/kappa;   //Distance at which correlation has dropped to 0.1, see p. 4 in Lindgren et al. (2011)

  REPORT(nll);
  REPORT(range);
  REPORT(Nhat);
  REPORT(Nhat_pred);
  REPORT(devresid);
  
  ADREPORT(range);
  ADREPORT(Nhat);
  ADREPORT(Nhat_pred);
  //ADREPORT(x); //this takes forever because there are potentially lots of x's

  return nll;
}
