#include "../include/Fractions.h"
#include "../include/Angular_Functions.h"
#include "../include/Index_Coefficients.h"
#include "../include/BM_Brackets.h"
#include "../include/Integrals.h"
#include <cmath>

// Functions and indexes Suhonen.
void inline lee_numeros_cuanticos(QN_2body_radial A){
    std::cout << "QN radial(n1 l1 n2 l2)=" << A.n1 << A.l1 <<A.n2 << A.l2 << std::endl;
}

double M_E_Spin_Isospin_Multipolar(const QN_2body_jj_Coupling & BRA, const QN_2body_jj_Coupling & KET,int lambda_max,  double b_param){
    
    // Constants for exchange operators and,
    double A = 1; // V_0
    double B = 0; // V_sigma
    double C = 0;//0.1; // V_tau
    double D = 0;//0.233; // V_tau_sigma
    
    double *Parameters;
    Parameters = new double[4];
    Parameters[0] = A;
    Parameters[1] = B;
    Parameters[2] = C;
    Parameters[3] = D;
    
    // Change from basic operators to the exchange.
    
    double mu_param = 1.; //1.48; // fm
    double V_0 = -1; //-70.82;

    // Exchanged wave function
    QN_2body_jj_Coupling KET_EXCH;
    KET_EXCH = {KET.n2, KET.l2, KET.j2,   KET.n1, KET.l1, KET.j1};

    double Aux_direct = 0.;
    double Aux_exchange = 0.;
    double aux_U;
    double rad_aux;
    
    int GL_order = 300;
    run_Gauss_Legendre(GL_order);
    
    // Reduced struct argument for the integrals.
    QN_2body_radial Q_N_radial_BRA;
    Q_N_radial_BRA = {BRA.n1, BRA.l1, BRA.n2, BRA.l2, 0, 0};// Total angular momentum and 3rd component are irrelevant

    QN_2body_radial Q_N_radial_KET, Q_N_radial_KET_EXCH;
    Q_N_radial_KET      = {KET.n1, KET.l1, KET.n2, KET.l2, 0, 0}; 
    Q_N_radial_KET_EXCH = {KET.n2, KET.l2, KET.n1, KET.l1, 0, 0};
    
    //lambda_max = min((BRA.j1 + KET.j1),(KET.j1 + KET.j2));
    
    for (int lambda = 0; lambda <= lambda_max; lambda ++){
        std::cout << " lambda =" << lambda << std::endl;
        // Computation of the direct Term V_abcd(JT)
        aux_U = U_Coeff(BRA, KET ,lambda,Parameters);
        std::cout <<" // U dir = "<<aux_U<<std::endl;
        if(fabs(aux_U) > 1.0e-10){

             //rad_aux = Radial_4_order(Q_N_radial_BRA, Q_N_radial_KET, lambda,b_param);
	     rad_aux = V_0 * Radial_4_order_Gauss_Quadrature(Q_N_radial_BRA, Q_N_radial_KET, GL_order ,lambda,mu_param,b_param);
             
	     Aux_direct +=  rad_aux * aux_U;
             std::cout <<"Direct: R= "<< rad_aux << std::endl;
             
	     //if(fabs(rad_aux) > 1e2){
             //   lee_numeros_cuanticos(Q_N_radial_BRA);
             //   lee_numeros_cuanticos(Q_N_radial_KET);
             //   std::cout <<"Direct: R= "<< rad_aux << std::endl;
             //}
        }

        //std::cout <<"Direct: R= "<< Radial_4_order(Q_N_radial_BRA, Q_N_radial_KET, lambda) << " // U = "<<U_Coeff(BRA, KET ,lambda)<<std::endl;
        // Computation of the exchange therm V_abdc(JT)
        aux_U = U_Coeff(BRA, KET_EXCH ,lambda,Parameters);
        std::cout <<" // U ex = "<< aux_U <<std::endl;
        if(fabs(aux_U) > 1.0e-10){
             
	     //rad_aux = Radial_4_order(Q_N_radial_BRA, Q_N_radial_KET_EXCH, lambda,b_param);
	     rad_aux = V_0 * Radial_4_order_Gauss_Quadrature(Q_N_radial_BRA, Q_N_radial_KET_EXCH, GL_order ,lambda, mu_param,b_param);
             Aux_exchange += rad_aux * aux_U;
             std::cout <<"Exchange: R= "<< rad_aux << std::endl;
             
	     //if(fabs(rad_aux) > 1e2){
             //   lee_numeros_cuanticos(Q_N_radial_BRA);
             //   lee_numeros_cuanticos(Q_N_radial_KET_EXCH);
             //   std::cout <<"Exchange: R= "<< rad_aux << std::endl;
             //}

        }
        //std::cout <<"Exchange: R= "<< Radial_4_order(Q_N_radial_BRA, Q_N_radial_KET_EXCH, lambda) << " // U = "<<U_Coeff(BRA, KET_EXCH ,lambda)<<std::endl;

    }

    // return Normalizated and antisimetrized.
    return Normalization_JT(BRA) * Normalization_JT(KET) *
	    ( Aux_direct + (pow(-1,double(KET.j1 + KET.j2  + KET.J) + KET.T) * Aux_exchange)); 

}

double M_E_SDI(const QN_2body_jj_Coupling & Q_Numbers_left,const QN_2body_jj_Coupling & Q_Numbers_right, bool isospin_exchange){

    // Exchange Forces constants
    double A = 1.;
    double C;
    if (isospin_exchange){
        C = 1.;
    }
    else{
        C = 0.;
    }

    //Passing of parameters

    int l_a = Q_Numbers_left.l1;
    Fraction j_a = Q_Numbers_left.j1;
    int l_b = Q_Numbers_left.l2;
    Fraction j_b = Q_Numbers_left.j2;

    Fraction J = Q_Numbers_left.J;
    int T = Q_Numbers_left.T;

    int l_c = Q_Numbers_right.l1;
    Fraction j_c = Q_Numbers_right.j1;
    int l_d = Q_Numbers_right.l2;
    Fraction j_d = Q_Numbers_right.j2;

    if ((Q_Numbers_left.M != Q_Numbers_right.M) || (Q_Numbers_left.M_T != Q_Numbers_right.M_T)){
        return 0;
    }
    if ((Q_Numbers_left.J != Q_Numbers_right.J) || (Q_Numbers_left.T != Q_Numbers_right.T)){
        return 0;
    }

    return  (A + C*(delta(T,1) - 3*delta(T,0)))*
            K_Coeff(Q_Numbers_left, Q_Numbers_right)* (1 + pow(-1,(l_a + l_b + l_c + l_d))) *
            Normalization_JT(Q_Numbers_left) * Normalization_JT(Q_Numbers_right) *
            sqrt((2*j_a + 1)*(2*j_b + 1)*(2*j_c + 1)*(2*j_d + 1)) *
            (
            ((1 + pow(-1,T))* Three_j_Symbol(j_a,j_b,J, Fraction(1,2), Fraction(1,2), Fraction(-1))*
                              Three_j_Symbol(j_c,j_d,J, Fraction(1,2), Fraction(1,2), Fraction(-1)))
             -
            ( ((1 - pow(-1, l_c + l_d + int((J + T)) ))* pow(-1,l_a + l_c + int((j_b + j_d)))) *
                              Three_j_Symbol(j_a,j_b,J, Fraction(1,2), Fraction(-1,2), Fraction(0))*
                              Three_j_Symbol(j_c,j_d,J, Fraction(1,2), Fraction(-1,2), Fraction(0)))
             );


}

double M_E_MSDI(const QN_2body_jj_Coupling & BRA, const QN_2body_jj_Coupling & KET){

    // Exchange Forces constants
    double B0 = -4.4;
    double B1 =  1.2;
    /*
     * Brussaard & Glaudemans book (1977)
     * 	 Parameters avaliable in Table 6.3 (pg 116)
     * 	  **  Considerations:
     *	 Change of variables with Suhonen Notation:
     * 	    A_T(suh) = 2 A_T(Bruss)
     * 	    B_1(suh) = B(Bruss) +   C(Bruss) 
     * 	    B_0(suh) = C(Bruss) - 3*B(Bruss)
     */

    // Building the j coupled functions
    QN_1body_jj  a,b,c,d;

    a = {BRA.n1, BRA.l1, BRA.j1};
    b = {BRA.n2, BRA.l2, BRA.j2};

    c = {KET.n1, KET.l1, KET.j1};
    d = {KET.n2, KET.l2, KET.j2};


    if ((BRA.M != KET.M) || (BRA.M_T != KET.M_T)){
        return 0;
    }
    if ((BRA.J != KET.J) || (BRA.T != KET.T)){
        return 0;
    }

    return (B1*delta(BRA.T,1) + B0*delta(BRA.T,0))*
            ( ((delta_2_1body_j(a,c)*delta_2_1body_j(b,d)) +
                (pow(-1,int((c.j + d.j + BRA.J) + BRA.T))*delta_2_1body_j(a,d)*delta_2_1body_j(b,c))
	      )/(1 + delta_2_1body_j(a,b)) )
             + M_E_SDI(BRA, KET, false);
            

}
