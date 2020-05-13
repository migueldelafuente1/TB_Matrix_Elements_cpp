#ifndef NUCLEAR_MATRIX_ELEMENTS_ART_H_INCLUDED
#define NUCLEAR_MATRIX_ELEMENTS_ART_H_INCLUDED

double M_E_Central_jj( const QN_2body_jj_Coupling & Q_Numbers_jj_Coupling_left,
                         const QN_2body_jj_Coupling & Q_Numbers_jj_Coupling_right,int option, double b_param);

// Modification of the coefficients just only as an extension of the functions ,
// Used just here, for the moment I keep them in this module.

double Coeff_CLS(int n1,int l1,int  n2, int l2, int lambda, int n1_q, int l1_q, int n2_q,int l2_q,int lambda_q,int option, double b_param);
double M_E_Spin_Orbit(const QN_2body_jj_Coupling & Q_Numbers_left, const QN_2body_jj_Coupling & Q_Numbers_right,int option, double b_param);

double Coeff_CT(int n1,int l1,int  n2, int l2, int lambda, int n1_q, int l1_q, int n2_q,int l2_q,int lambda_q,int option, double b_param);
double M_E_Tensor( const QN_2body_jj_Coupling & Q_Numbers_left, const QN_2body_jj_Coupling & Q_Numbers_right,int option, double b_param);

#endif // NUCLEAR_MATRIX_ELEMENTS_ART_H_INCLUDED
