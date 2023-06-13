// 2d.cpp

#include <TMB.hpp>

template <class Type>
Type objective_function<Type>::operator()()
{
  PARAMETER(theta1);
  PARAMETER(theta2);

  // Initialise negative log-likelihood
  Type nll;
  nll = Type(0.0);

  return(nll);
}
