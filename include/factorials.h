#include "Fractions.h"

#ifndef FACTORIALS_H_INCLUDED
#define FACTORIALS_H_INCLUDED
#pragma once

extern int Factorial_Maximum_Dimension;
extern double * Factorial_int_list;
extern double * Factorial_hfi_list;
extern double * Double_Factorial_list;

void run_Factotial_Base();
double factorial(Fraction n);
double factorial(int n);

double factorial_function(int * numerator, int * denominator,int N_num,int N_den);

double factorial_function(Fraction * num, Fraction * den,int N_num,int N_den);

double double_factorial(int n);

double double_factorial_function(int * numerator, int * denominator,int N_num,int N_den);

double gamma_function(Fraction x, int &sign);

#endif // FACTORIALS_H_INCLUDED
