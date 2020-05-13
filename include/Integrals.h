#ifndef INTEGRALS_H_INCLUDED
#define INTEGRALS_H_INCLUDED

#include "Index_Coefficients.h"
extern double oscillator_lenght;

double Talmi_Integral(int p, int option ,double b_param);

//struct QN_2body_radial;

double Multipolar_Rosenfeld(int lambda, int GL_ORDER, double r_1, double r_2, double mu_param);

double Radial_4_order_Subdivided_Quadrature(const QN_2body_radial &Q_Numbers_left,
				       const QN_2body_radial &Q_Numbers_right, int ORDER, 
				       int lambda, double mu_param, double b_param);

double Radial_4_order_Gauss_Quadrature(const QN_2body_radial &Q_N_radial_left,
				       const QN_2body_radial &Q_N_radial_right, int ORDER, 
					int lambda, double mu_param, double b_param);

double Radial_4_order(const QN_2body_radial &left, const QN_2body_radial &right, int lambda, double b_param);

double Rosenfeld_quadratures(int n_a,int l_a,int n_b,int l_b,int L,int n_c,int l_c,int n_d,int l_d,
			      int Lq,int GL_ORDER,double mu_param ,double b_param);

double Radial_SDI(const QN_2body_radial &WF);
#endif // INTEGRALS_H_INCLUDED
