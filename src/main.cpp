#include <iostream>

#include "../include/factorials.h"
#include "../include/Index_Coefficients.h"
#include "../include/Angular_Functions.h"
#include "../include/Wave_Functions.h"
#include "../include/BM_Brackets.h"
#include "../include/Fractions.h"
#include "../include/Integrals.h"

#include "../include/checkers.h"
#include "../include/output.h"

#include "../include/Nuclear_Matrix_Elements_M.h"
#include "../include/Nuclear_Matrix_Elements_Art.h"
#include "../include/Nuclear_Matrix_Elements_Suh.h"
#include "../include/Nuclear_Matrix_Elements_BB.h"

#define _USE_MATH_DEFINES
#include <iostream>
#include <fstream>
#include <math.h>
#include <iomanip>
#include <stdlib.h>
#include <time.h>
#include <string>
#include <sstream>
//#pragma once
using namespace std;


int main()
{
    run_Factotial_Base();

    //double Bracket=BMB(n, l, N, L, lambda, n1, l1, n2, l2);
    //cout<<"Monshinsky ( )="<<BMB(1,3,6,1,4,2,4,2,6)<<endl;
    //cout<<"Monshinsky ( )="<<BMB(0,3,5,2,3,2,4,2,4)<<endl;
    //cout<<"Monshinsky ( )="<<BMB(1,3,4,2,3,2,4,2,4)<<endl;
    //cout<<"Monshinsky ( )="<<BMB(3,2,3,4,4,2,4,2,6)<<endl;
    //check_BMB_article_values();
    //bool continue_bool;
    //cin >> continue_bool;
    
    /*
    QN_2body_jj_Coupling BRA, KET;
    BRA = {0,2,Fraction(3,2), 0,2 ,Fraction(3,2), 0,0, 1,0};
    KET = BRA;
    //KET = {0,1,Fraction(1,2), 0,1,Fraction(1,2), 0,0,1,0};
    print_QN_2body_jj(BRA);
    print_QN_2body_jj(KET);
    bool Monsinsky_Method = true;
    cout << "B.B matrix element = " << M_E_Brink_Boeker(BRA, KET, 1.32) <<endl;
    //*/

    //std::cout<<"(3%2)="<<(3%2)<<"     (2%2)="<<(2%2)<<std::endl;
    //QN_1body_radial bra = {1,2,1};
    //bool Correct = check_Chasman_Conversion_Spherical2Cartesian(bra);
    
    
    
    cout << Clebsh_Gordan(Fraction(1), Fraction(3,2), Fraction(1,2),Fraction(-1), Fraction(1,2), Fraction(-1,2) )
         <<"   / alternado=" << Clebsh_Gordan(Fraction(3,2),Fraction(1),  Fraction(1,2), Fraction(1,2),Fraction(-1), Fraction(-1,2) ) <<endl;
    
    
    QN_1body_radial a = {0, 0, 0};
    QN_1body_radial b = {0, 1, 0};
    QN_1body_radial c = a;
    QN_1body_radial d = b;
    
    QN_2body_radial Q_N_radial_BRA;
    Q_N_radial_BRA = {0, 1, 0, 2, 0, 0};
    QN_2body_jj_Coupling BRA = {0,1,Fraction(3,2), 0,1,Fraction(3,2), Fraction(3),Fraction(0),0,0};

    QN_2body_radial Q_N_radial_KET;
    Q_N_radial_KET =  Q_N_radial_BRA; // {KET.n1, KET.l1, KET.n2, KET.l2, 0, 0}; //
    QN_2body_jj_Coupling KET = BRA; // {0,1,Fraction(3,2), 0,2,Fraction(3,2), Fraction(0),Fraction(0),1,0};  //
    
    double mu_param = 1.;
    double b_param  = 1.;
    
    
    int M = 0;
    BRA = {0,1,Fraction(1,2), 0,1,Fraction(1,2), Fraction(1),Fraction(M),0,0};
    Q_N_radial_BRA = {0, 0, 0, 1, 0, 0};
    KET = {0,1,Fraction(1,2), 0,1,Fraction(1,2), Fraction(1),Fraction(M),0,0};
    Q_N_radial_KET = {0, 0, 0, 1, 0, 0};
    
    
    check_Coupling_BBME_Methods(BRA, KET,mu_param,b_param,50, 10, true);
    
    /*Fraction J = Fraction(2);
    for(M = -int(J); M <= int(J); M++ ){
	BRA = {0,1,Fraction(3,2), 0,1,Fraction(3,2), Fraction(J),Fraction(M),0,0};
	check_LSjj_coupling_coeffs( BRA,  BRA);
    }*/
    //cout << "Radial BB=" << Radial_BB(a,b,c,d,mu_param,b_param) << endl;
    //cout << "M.E.Multipolar=" << M_E_Brink_Boeker(BRA,KET,3,b_param) << endl;
    //BRA.J = 2;
    //KET.J = 2;
    //cout << "M.E.Multipolar=" << M_E_Brink_Boeker(BRA,KET,3,b_param) << endl;
    
    //cout << "Radial BB=" << Non_antisimetriced_BB_Multi(BRA,KET,8,b_param) << endl;
    
    //cout << Nine_j_Symbol( Fraction(2), Fraction(1,2),Fraction(3,2),
    //                        Fraction(0),Fraction(1,2),Fraction(1,2), Fraction(2),  Fraction(0),  Fraction(2))<< endl;

    
    /*
    double gl_basic,gl_upgrade, Talman, Moshinsky ;
    Talman = Radial_BB(a, b, c,  d ,  mu_param,  b_param);
    Moshinsky = Radial_BB_Moshinsky(a, b, c,  d ,  mu_param,  b_param * sqrt(2));    
    cout << "Talman= " << Talman << " , Moshinsky= " << Moshinsky<< endl;
    
    
    for (int i = 2; i<100;i++){
	cout << "("<<i<<") "endl;
	run_Gauss_Legendre(i);
	for(int lambda = 0; lambda <= 10; lambda++){
	    gl_basic = Multipolar_Rosenfeld(lambda, i, 0.0400272, 0.0400272,  mu_param);
	    //cout << "=" << Multipolar_Rosenfeld(0, 0, 0.000400272, 6.31614,  mu_param)<< endl;
	    //Radial_4_order_Gauss_Quadrature(Q_N_radial_BRA, Q_N_radial_KET, i ,1,mu_param,b_param);
	    cout << gl_basic << " ";
	}
	cout << endl;
    }
    //*/
    
    
    //QN_1body_jj WF = {0,5,Fraction(11,2), Fraction(0)};
    //check_ME_multipolar(WF, 3);
    
    //WF = {1,1,Fraction(1,2),Fraction(-1,2)};
    //check_SHO_WF(WF , 15, b_lenght(40));
    
    //check_Cuadrature_Convergence(Q_N_radial_BRA,Q_N_radial_KET, 0,mu_param,b_param,300,8);
    
    /*
    double rad_aux, ang_aux;
    for (int lambda = 0; lambda <= 10; lambda++){
	rad_aux = Radial_4_order_Gauss_Quadrature(Q_N_radial_BRA, Q_N_radial_KET, lambda, mu_param,b_param);
	ang_aux = U_Coeff(BRA,KET, lambda);
	cout << "rad_aux (lambda = " << lambda <<") = "<< rad_aux << "/ * / ang_aux = "<< ang_aux << " /=/ "<< ang_aux*rad_aux << endl;}
    //*/
    
    /*double dr = 20./40;
    for( int i = 0; i <= 40; i++){
	cout << "SHO = "<< SHO_Radial(11, 1, dr*i ,  b_param) <<endl; 
    }*/
    
    
    //cout <<"Spin Orbit M.E. ="<< M_E_Spin_Isospin_Multipolar(BRA,KET,0)<< endl;
    //cout <<"   // SDI  M.E. =" << M_E_SDI(BRA, KET,true) << endl;
    //cout <<endl;
    //cout << "M E Spin Orbit_jj =" << M_E_Spin_Orbit_jj(BRA, KET)<< endl;
    //cout << "numeric M Tensor =" << M_E_Tensor_jj(BRA, KET) << endl;
    //cout <<endl;
    //cout << "numeric Art Central =" << M_E_Central_jj(BRA, KET)<< endl;
    //cout << "numeric Art LS =" << M_E_Spin_Orbit(BRA, KET)<< endl;
    //cout << "numeric Art Tensor =" << M_E_Tensor(BRA, KET) << endl;

    //QN_1body_jj WF = {1,3,Fraction(7,2)};
    //check_ME_multipolar(WF,6);
    //check_ME_Art(3, 2);
    //check_ME_SDI();
    
    /*
    cout << "B_coefficient(1, 0, 0, 2, (1)) = " << B_coefficient(1, 0, 0, 2, Fraction(1), 1) << endl;    
    cout << "B_coefficient(0, 1, 0, 1, (0)) (non)= " << B_coefficient(0, 1, 0, 1, Fraction(0), 1) << endl;
    cout << "B_coefficient(0, 1, 0, 1, (1)) = " << B_coefficient(0, 1, 0, 1, Fraction(1), 1) << endl;
    cout << "B_coefficient(1, 1, 1, 1, (1)) = " << B_coefficient(1, 1, 1, 1, Fraction(1), 1) << endl;
    cout << "B_coefficient(1, 1, 1, 1, (2)) = " << B_coefficient(1, 1, 1, 1, Fraction(2), 1) << endl;
    */
    
    /*
    cout<<"<<BMB(0,0,0,0,1,0,1,0,2)= "<< BMB(0,0,0,0,1,0,1,0,2)<<endl;
    cout<<"<<BMB(1,0,0,1,1,0,1,0,2)= "<< BMB(1,0,0,1,1,0,1,0,2)<<endl;
    cout<<"Monshinsky ( )="<<BMB(1,3,0,6,6,0,5,0,6)<<endl;
    
    QN_2body_radial bra = {0,1,0,2,1, 0};
    QN_2body_radial ket = bra;
    
    QN_1body_radial a = {0,1,0};
    QN_1body_radial b = {0,2,0};
    QN_1body_radial c = a;
    QN_1body_radial d = b;
    
    double b_param = 1;
    double mu_param = sqrt(2);
    int lambda = 1;
    
    Set_Manual_Lambda = true;
    //LAMBDA = 1./(sqrt(nu) * mu_param);
    LAMBDA = (sqrt(2) * b_param) / mu_param;
    
    cout << "M_E_Central =" << M_E_Central(a.n, a.l, b.n, b.l, lambda, c.n, c.l, d.n, d.l, lambda, b_param*sqrt(1)) <<endl;
    
    //Monsinsky_Method = true;
    //cout << " Mosinsky Radial Central = " << Radial_BB_Moshinsky (a , b, c, d, mu_param, b_param) << endl;
    //Monsinsky_Method = false;
    //cout << " Talman   Radial Central = " << Radial_BB(a , b, c, d, mu_param, b_param) << endl;
    //*/
    
    //antoine_output(1);

    return 0;
}
