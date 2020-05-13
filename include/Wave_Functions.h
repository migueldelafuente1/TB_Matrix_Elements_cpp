#ifndef WAVE_FUNCTIONS_H_INCLUDED
#define WAVE_FUNCTIONS_H_INCLUDED

double SHO_normalization(int n, int l, double b_param);
double SHO_Radial(int n, int l, double r, double b_param);
double Legendre_Polynomial (int n, double  x);

#endif // WAVE_FUNCTIONS_H_INCLUDED
