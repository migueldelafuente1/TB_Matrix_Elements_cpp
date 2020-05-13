#ifndef NUCLEAR_MATRIX_ELEMENTS_M_H_INCLUDED
#define NUCLEAR_MATRIX_ELEMENTS_M_H_INCLUDED

#include "Index_Coefficients.h"

double R_M_E_Potential(int n, int l, int n_q, int l_q,int option, double b_param);

double M_E_Central(int n1, int l1, int n2, int l2, int lambda, int n1_q, int l1_q, int n2_q, int l2_q,
								 int lambda_q,int option, double b_param);

double M_E_Spin_Orbit(int n1, int l1, int n2, int l2, int lambda, int n1_q, int l1_q, int n2_q, int l2_q,
								 int lambda_q,int option, double b_param);
double M_E_Spin_Orbit_jj(const QN_2body_jj_Coupling & BRA,const QN_2body_jj_Coupling & KET,int option, double b_param);

double M_E_Tensor(int n1, int l1, int n2, int l2, int lambda, int n1_q, int l1_q, int n2_q, int l2_q, 
								int lambda_q,int option, double b_param);
double M_E_Tensor_jj(const QN_2body_jj_Coupling & BRA,const QN_2body_jj_Coupling & KET,int option, double b_param);

#endif // NUCLEAR_MATRIX_ELEMENTS_M_H_INCLUDED
