#include <TMB.hpp>


template<class Type>
Type objective_function<Type>::operator() ()
{
  DATA_VECTOR(y);
  DATA_VECTOR_INDICATOR(keep, y);
  DATA_VECTOR(age);
  DATA_VECTOR(N1);
  DATA_VECTOR(N2);
  DATA_IVECTOR(age_indx);
  DATA_VECTOR(NbyYear);
  DATA_VECTOR(weight);

  PARAMETER_VECTOR(Linf);
  PARAMETER_VECTOR(rate);
  PARAMETER(t0);
  PARAMETER_VECTOR(logSigma);

  vector<Type> mu(y.size());
  vector<Type> residual(y.size());  
  vector<Type> SD(y.size());  
  vector<Type> sigma = exp(logSigma);
  
  Type nll = 0;
  for(int i = 0; i < y.size();  i++){
      mu(i) = (Linf(0) + Linf(1) * N1(i)) * (Type(1.0) - exp(-(rate(0) + rate(1) * N2(i)) * (age(i) - t0)));
      residual(i) = y(i)-mu(i);
      SD(i) = exp(logSigma(age_indx(i)));
      nll  -= keep(i) * weight(i) * dnorm(residual(i), Type(0.0), SD(i), true);
  }
  
  // SIMULATE{
  //   
  //   for(int i = 0; i < gamma.size(); i++){
  //     gamma(i) = rnorm(Type(0.0), exp(logSigmaGam(0)));
  //   }
  //   
  //   for(int i = 0; i < y.size(); i++){
  //     if(yisFishLength == 1){
  //       mu(i) = (Linf(0) + Linf(1) * N(i)) * (Type(1.0) - exp(-(rate(0) + rate(1) * N(i)) * (age(i) - (t0(0) + t0(1) * N(i))) ));
  //     }else if(yisFishLength == 0){
  //       mu(i) = (Linf(0) + Linf(1) * N(i)) * (Type(1.0) - exp( -(rate(0) + rate(1) * N(i)) * (age(i) - (t0(0) + t0(1) * N(i))) ))*
  //         (Type(1.0) - exp( -(rate(0) + rate(1) * N(i)) * (age(i) - (t0(0) + t0(1) * N(i))) ))*
  //         (Type(1.0) - exp( -(rate(0) + rate(1) * N(i)) * (age(i) - (t0(0) + t0(1) * N(i))) ));
  //     }else{
  //       mu(i) = Linf(0);
  //     }
  //     if(dist == 0){
  //       y(i) = rnorm(mu(i) + gamma(serial_indx(i)) , exp(logSigma(age_indx(i)))); 
  //     }
  //     if(dist == 1){
  //       y(i) = mu(i) + gamma(serial_indx(i))+rt(Type(2.0) + exp(logdegf(age_indx(i))));
  //     }
  //   }
  //   REPORT(y);
  //   
  // }
  vector<Type> predLinf(NbyYear.size());  
  for(int i = 0; i < NbyYear.size(); i++) {
    predLinf(i) = Linf(0) + Linf(1) * NbyYear(i);
  }
  ADREPORT(predLinf);
  ADREPORT(sigma);
  REPORT(mu);
  REPORT(residual);
  REPORT(SD);
  return nll;
  }

