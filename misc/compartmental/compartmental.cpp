// compartmental.cpp

#include <TMB.hpp>

template<class Type>
vector<Type> dz_dt(vector<Type> Z, Type t, Type b, Type g) {
  Type S = Z[0];
  Type I = Z[1];
  Type R = Z[2];

  Type prev = I / (Z.sum());
  Type dS_dt = -b * S * prev;
  Type dI_dt = b * S * prev - g * I;
  Type dR_dt = g * I;

  vector<Type> out(3);
  out << dS_dt, dI_dt, dR_dt;

  return out;
}

// Objective function (TMB requirement)
template<class Type>
Type objective_function<Type>::operator() ()
{
  DATA_INTEGER(N_t);
  DATA_SCALAR(dt);

  PARAMETER(b);
  PARAMETER(g);
  PARAMETER(I_0);

  Type nll = 0.0;

  matrix<Type> Z_out;
  Z_out.setZero(N_t, 3);

  vector<Type> Z(3);
  Z << 1.0 - I_0, I_0, 0.0;

  vector<Type> dZ(3);
  for (size_t i = 0; i < N_t; i++) {
    dZ = dz_dt(Z, Type(i), b, g);
    Z = Z + dt * dZ;
    Z_out.row(i) = Z;
  }

  REPORT(Z_out);

  return nll;
}
