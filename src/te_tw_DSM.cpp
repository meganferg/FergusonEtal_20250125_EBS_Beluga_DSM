//Script te_tw_DSM.cpp...Megan C. Ferguson...27 January 2025

  //NOTES
  //
  //1. This script was based on Skaug's pSplines.cpp and sample code for
  //   a gam with a bivariate tensor product smooth of x and y with shrinkage
  //   that I created using DLM's mgcvminusminus package
  //   (https://github.com/dill/mgcvminusminus).
  //
  //2  See https://kaskr.github.io/adcomp/tweedie_8cpp-example.html for example
  //   of a Tweedie distribution.
  //
  //3. This script is identical to NSDL17dsm_tmb_tw_te_PDE.cpp in MCF's original files.

#include <TMB.hpp>   //Links in the TMB libraries
  
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
Type objective_function<Type>::operator() (){
  using namespace density;
  using namespace Eigen; //Needed for utilisation of sparse structures
  
  //Read data from R------------------
  DATA_VECTOR(y);       //The response
  DATA_MATRIX(X);       //Design matrix for splines
  DATA_MATRIX(X_pred); // design matrix for density surface predictions
  DATA_VECTOR(offset);  //offset term for model-building
  DATA_VECTOR(offset_pred); //offset term for density surface predictions
  DATA_SPARSE_MATRIX(S1);//Penalization matrix, converted to a sparse matrix.
  DATA_SPARSE_MATRIX(S2);//Penalization matrix, converted to a sparse matrix.
  //----------------------------------
  
  //Read parameters from R------------
  PARAMETER(beta0);       //Intercept
  PARAMETER_VECTOR(beta); //Spline regression parameters
  PARAMETER_VECTOR(log_lambda);//Penalization parameters
  PARAMETER(log_phi);     // log(error SD)
  PARAMETER(finv_power);     // transformed power
  //----------------------------------
  
  //Derived variable declarations--
  vector<Type> lambda = exp(log_lambda);
  Type phi = exp(log_phi);
  Type power = 1.0 + (exp(finv_power) / (1.0 + exp(finv_power)));
  Type Nhat = 0; // abundance derived from sampled segments in DATA
  Type Nhat_pred = 0; // abundance derived from density surface
  vector<Type> devresid( y.size() ); // deviance residual vector
  //----------------------------------
  
  // define the real penalties
  Eigen::SparseMatrix<Type> K1 = lambda(0)*S1 + lambda(1)*S2;

  // initialize log-likelihood
  Type nll=0;

  using namespace atomic; // for logdet

  // calculate REML penalty
  nll = -0.5*(logdet(matrix<Type>(K1))) +
         0.5*(GMRF(K1).Quadform(beta));
  
  vector<Type> mu(y.size());
  mu = beta0 + X*beta; // linear predictor for sampled segments;
                       // without offset; on link scale
  
  for(int i=0; i<y.size(); i++){
    nll -= dtweedie( y(i), offset(i)*exp(mu(i)), phi, power, true );
    Nhat += offset(i)*exp(mu(i));
    devresid(i) = devresid_tweedie( y(i), offset(i)*exp(mu(i)), power ); // deviance residual
  }

  // Compute Nhat_pred, predicted abundance based on density surface
    vector<Type> mu_pred = beta0 + X_pred*beta; // linear predictor for predictions on density surface;
                                              //  without offset; on link scale
  
    for(int i=0; i<offset_pred.size(); i++){
      Nhat_pred += offset_pred(i)*exp(mu_pred(i));
    }
  
  // Reporting
    REPORT( nll );
    REPORT( Nhat );
    REPORT( Nhat_pred );
    REPORT( devresid );
  
    // ADREPORT( beta );
    ADREPORT( Nhat );
    ADREPORT( Nhat_pred );

  return nll;
}
