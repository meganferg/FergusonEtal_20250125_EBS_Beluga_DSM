//Script soap_tw_DSM.cpp...Megan C. Ferguson...25 January 2025

  //NOTES
  //
  //1. This script was based on Skaug's mgcv_soap.cpp, with modifications
  //   by Jim Thorson shown in mgcv_soap_JT.cpp.
  //
  //2  See https://kaskr.github.io/adcomp/tweedie_8cpp-example.html for example
  //   of a Tweedie distribution.
  //
  //3. This script is identical to NSDL17dsm_tmp_PDE.cpp in MCF's original files.

#include <TMB.hpp>
  
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
  using namespace atomic;

  // DATA declarations
  
    DATA_VECTOR(y);  // response variable
    
    // mgcv penalty matrices
      DATA_MATRIX(S1);  // boundary
      DATA_MATRIX(S2);  // internal
      
    DATA_MATRIX(X);  // design matrix for model-building
    DATA_MATRIX(X_pred); // design matrix for density surface predictions
    
    DATA_VECTOR(offset); //offset term for model-building
    DATA_VECTOR(offset_pred); //offset term for density surface predictions
  
    DATA_INTEGER(n_space_dim);  // dimension of null of smoothers (zero in this case)
  
  // PARAMETER declarations 
  
    PARAMETER(mu);     // intercept
    PARAMETER_VECTOR(beta);     // spline coefficients
    PARAMETER_VECTOR(log_lambda);     // log(smoothing parameters)
    PARAMETER(log_phi);     // log(error SD)
    PARAMETER(finv_power);     // transformed power
  
  // Derived variable declarations

    vector<Type> lambda = exp(log_lambda);
    Type phi = exp(log_phi);
    Type power = 1.0 + (exp(finv_power) / (1.0 + exp(finv_power)));
  
    Type nll=0; // negative log-likelihood
    
    matrix<Type> S = lambda(0)*S1+lambda(1)*S2;  // see Wood (2011)
    
    vector<Type> eta = mu + X*beta; // linear predictor for sampled segments;
                                    // without offset; on link scale
    vector<Type> eta_pred = mu + X_pred*beta; // linear predictor for predictions on density surface;
                                              //  without offset; on link scale
    
    Type Nhat = 0; // abundance derived from sampled segments in DATA
    Type Nhat_pred = 0; // abundance derived from density surface
    
    vector<Type> devresid( y.size() ); // deviance residual vector
    
  // Prior on the smoother induced by the penalty
  // From Wood (2011). Note that S has full rank for the soap film smoother.
    nll -= 0.5*logdet(S) - 0.5*(beta* vector<Type>(S*beta)).sum();  
  
  // Likelihood contribution from the observations
    
    //nll -= dnorm(y,eta,sigma,true).sum(); // Use this for a Gaussian model
    
    // Use the following for a Tweedie model
      for( int i=0; i<y.size(); i++ ){
        nll -= dtweedie( y(i), offset(i)*exp(eta(i)), phi, power, true );
        Nhat += offset(i)*exp(eta(i));
        devresid(i) = devresid_tweedie( y(i), offset(i)*exp(eta(i)), power ); // deviance residual
      }
      
  // Compute Nhat_pred, predicted abundance based on density surface
    for(int i=0; i<offset_pred.size(); i++){
      Nhat_pred += offset_pred(i)*exp(eta_pred(i));
    }
    
  // Reporting
    REPORT( nll );
    REPORT( Nhat );
    REPORT( Nhat_pred );
    REPORT( devresid );
    ADREPORT( Nhat );
    ADREPORT( Nhat_pred);
    
  return nll;
}
