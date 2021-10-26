#include <TMB.hpp>

template<class Type>
Type objective_function<Type>::operator() ()
{
  DATA_VECTOR(y);
  DATA_VECTOR_INDICATOR(keep, y);
  DATA_VECTOR(age);
  DATA_VECTOR(weights);

  PARAMETER(logLinf);
  PARAMETER(lograte);
  PARAMETER(t0);
  PARAMETER(logSigma);
  
  Type Linf = exp(logLinf);
  Type rate = exp(lograte);
  Type sigma = exp(logSigma); 
  
  vector<Type> mu(y.size());
  vector<Type> residual(y.size());  
  
  Type nll = 0;
  //Type x = 0;
  for(int i = 0; i < y.size();  i++){
    mu(i) = Linf * (Type(1.0) - exp(-rate * (age(i) - t0)));
    residual(i) = (y(i)-mu(i));
    nll  -= keep(i) * weights(i) * dnorm(residual(i), Type(0.0), sigma, true);
  }
  ADREPORT(Linf);
  ADREPORT(rate);
  ADREPORT(sigma);
  REPORT(mu);
  REPORT(residual);
  return nll;
}

