#ifndef NUCLEAR_MATRIX_ELEMENTS_SUH_H_INCLUDED
#define NUCLEAR_MATRIX_ELEMENTS_SUH_H_INCLUDED

#include "Index_Coefficients.h"

double M_E_Spin_Isospin_Multipolar(const QN_2body_jj_Coupling & BRA, const QN_2body_jj_Coupling &KET,
					int lambda_max, double b_param);

double M_E_SDI(const QN_2body_jj_Coupling & BRA, const QN_2body_jj_Coupling & KET, bool isospin_exchange);

double M_E_MSDI(const QN_2body_jj_Coupling & BRA, const QN_2body_jj_Coupling & KET);

#endif // NUCLEAR_MATRIX_ELEMENTS_SUH_H_INCLUDED
