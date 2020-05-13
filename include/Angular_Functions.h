#include "Fractions.h"

#ifndef ANGULAR_FUNCTIONS_H_INCLUDED
#define ANGULAR_FUNCTIONS_H_INCLUDED

bool triangular_condition(int a, int b, int c);
bool triangular_condition(Fraction a, Fraction b, Fraction c);

double Delta(int a,int b, int c);
double Delta(Fraction a, Fraction b, Fraction c);

double Racah_Coefficient(int a,int b, int c, int d, int e, int f);
double Racah_Coefficient(Fraction a,Fraction b, Fraction c, Fraction d, Fraction e, Fraction f);

double Racah_Coeficient_Analitic(int a_int,int b_int, int c_int, int d_int, int e_int, int f_int);

double Six_j_Coefficient(int j1,int j2, int j3, int J1, int J2, int J3);
double Six_j_Coefficient(Fraction j1,Fraction j2, Fraction j3, Fraction J1, Fraction J2, Fraction J3);

double U_Coefficient(int a,int b, int c, int d, int e, int f);
double U_Coefficient(Fraction a,Fraction b, Fraction c, Fraction d, Fraction e, Fraction f);

double Clebsh_Gordan(Fraction j1, Fraction j2, Fraction j, Fraction m1, Fraction m2, Fraction m);
double Clebsh_Gordan(int j1, int j2, int j, int m1, int m2, int m);

double Three_j_Symbol(int j1, int j2, int j, int m1, int m2, int m);
double Three_j_Symbol(Fraction j1, Fraction j2, Fraction j, Fraction m1, Fraction m2, Fraction m);

double Nine_j_Symbol(int l1, int s1, int j1, int l2, int s2, int j2, int L, int S, int J);
double Nine_j_Symbol(Fraction l1, Fraction s1, Fraction j1, Fraction l2, Fraction s2, Fraction j2, Fraction L, Fraction S, Fraction J);

double LS_jj_coupling_Coeff(Fraction l1, Fraction s1, Fraction j1, Fraction l2, Fraction s2, Fraction j2, Fraction L, Fraction S, Fraction J);
double LS_jj_coupling_Coeff(int l1, int s1, int j1, int l2, int s2, int j2, int L, int S, int J);

#endif // ANGULAR_FUNCTIONS_H_INCLUDED
