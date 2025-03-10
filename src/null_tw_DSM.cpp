//Script null_tw_DSM.cpp...Megan C. Ferguson...25 January 2025
//
//  This script is identical to MCF's original NSDL1722_null_tw_PDE.cpp

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
    DATA_VECTOR(offset); //offset term for model-building
    
  // PARAMETER declarations 
  
    PARAMETER(mu);     // intercept
    PARAMETER(log_phi);     // log(error SD)
    PARAMETER(finv_power);     // transformed power
  
  // Derived variable declarations

    Type phi = exp(log_phi);
    Type power = 1.0 + (exp(finv_power) / (1.0 + exp(finv_power)));
  
    Type nll=0; // negative log-likelihood
    
    vector<Type> devresid( y.size() ); // deviance residual vector
    
  // Likelihood contribution from the observations
    
    //nll -= dnorm(y,eta,sigma,true).sum(); // Use this for a Gaussian model
    
    // Use the following for a Tweedie model
      for( int i=0; i<y.size(); i++ ){
        nll -= dtweedie( y(i), offset(i)*exp(mu), phi, power, true );
        devresid(i) = devresid_tweedie( y(i), offset(i)*exp(mu), power ); // deviance residual
      }
      
  // Reporting
    REPORT( nll );
    REPORT( devresid );

  return nll;
}
