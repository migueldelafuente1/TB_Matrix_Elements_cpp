#ifndef NUCLEAR_MATRIX_ELEMENTS_BB_H_INCLUDED
#define NUCLEAR_MATRIX_ELEMENTS_BB_H_INCLUDED
#include <complex>
#pragma once
/*
extern double mu_BB[2];

extern double A_BB [2];
extern double B_BB [2];
extern double C_BB [2];
extern double D_BB [2];
//*/
extern double LAMBDA;
extern bool Moshinsky_Method;

double M_E_Brink_Boeker(const QN_2body_jj_Coupling & Q_Numbers_left,
                        const QN_2body_jj_Coupling & Q_Numbers_right, int option,  double b_param);

 // ** Three ways of coupling for the average of the operators 
double Non_antisimetriced_BB_M_E(const QN_2body_jj_Coupling & Q_Numbers_left,
			  	 const QN_2body_jj_Coupling & Q_Numbers_right, double b_param);
 // and for exchange operators	
double Non_antisimetriced_BBME_by_exchange(const QN_2body_jj_Coupling & Q_Numbers_left, 
					   const QN_2body_jj_Coupling & Q_Numbers_right, double b_param);

double Non_antisimetriced_BBME_by_exchange_LSjj(const QN_2body_jj_Coupling & Q_Numbers_left,
					   const QN_2body_jj_Coupling & Q_Numbers_right, double b_param);
 // and the numerical approach
double Non_antisimetriced_BB_Multi(const QN_2body_jj_Coupling & BRA, const QN_2body_jj_Coupling & KET,int lambda_max, double b_param);

std::complex<double> conversion_NLM_XYZ(const QN_1body_radial & spherical, int n_x,int n_y,int n_z);

double Coeff_Talman(int n_a, int n_b, int n_mu);

double J_integral(int n_b, int n_d, int n_mu, double mu, double b);

double Radial_BB(const QN_1body_radial & a, const QN_1body_radial & b,
                  const QN_1body_radial & c, const QN_1body_radial & d , double mu, double b_param);

double Radial_BB_Moshinsky(const QN_1body_radial & a, const QN_1body_radial & b,
                  const QN_1body_radial & c, const QN_1body_radial & d , double mu, double b_param);

#endif // NUCLEAR_MATRIX_ELEMENTS_BB_H_INCLUDED
