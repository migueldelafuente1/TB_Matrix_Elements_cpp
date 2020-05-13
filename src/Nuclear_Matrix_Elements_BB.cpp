#include "../include/Fractions.h"
#include "../include/factorials.h"
#include "../include/Angular_Functions.h"
#include "../include/Index_Coefficients.h"
#include "../include/BM_Brackets.h"
#include "../include/Integrals.h"
#include "../include/checkers.h"

#define _USE_MATH_DEFINES
#include <cmath>
#include <complex>

#include "../include/Nuclear_Matrix_Elements_M.h"
#include "../include/Nuclear_Matrix_Elements_BB.h"
#include "../include/Nuclear_Matrix_Elements_Suh.h"

/**     Brink Brooker Matrix elements, I have create an additional space
 *      for all the functions involved in the process:
 *          *  All Changes of Basis, to uncouple the |a b (JT)>  w.f
 *          *  Specific change of Basis for n l m to n_x n_y n_z
 *          *  Talman coefficients and J integral to perform the radial integral
 *      for the gaussian interaction
 *
 *      BB interaction or Gogny force employ a pair of gaussians and their corresponding
 *      set of parameters for each one. Must be defined
 *
 */

/*
double mu_BB[0] =  1;
double mu_BB[1] =  1;

double A_BB[0] = 1;
double A_BB[1] = 1;

double B_BB[0] = 1;
double B_BB[1] = 1;

double C_BB[0] = 1;
double C_BB[1] = 1;

double D_BB[0] = 1;
double D_BB[1] = 1;
*/

//*

/** Constant specifications, there are two forms of describing the spin-isospin
 *  exchange, which are by Permutation operators (with the constants V_W, V_M,
 *  V_B and V_H) and by the dot product of the two spin or isospin operators (A B C & D).
 *  The code is performed with the latter ones, but the constants are usually for the former ones.
 *      These are related by a linear transformation.
 *      Uncomment to perform the computation for specific A,B,C & D constants
 */

double mu_BB[2] = { 0.7 , 1.4 };
//double mu_BB[2] = { 1.48, 1.48 };
//double mu_BB[2] = { 1., 0 };

/////////    Permutation Constants conversion to ABCD   ////////////////
// Parameters for Rosenfeld
/*
double V_W_BB[2] = { -1., 0. };    //{ 595.55 , -72.21 }; //  // Wigner, Central
double V_M_BB[2] = { 0. , 0. };         // Majorana, Spatial Permutation
double V_B_BB[2] = { 0. , 0. };         // Barlett, Spin Permutation
double V_H_BB[2] = { 0. , 0. };         // Heisenberg, Isospin Permutation
//*/
// Parameters for hw =23.8

double V_W_BB[2] = { 714.085  , -72.6119 };    //{ 595.55 , -72.21 }; //  // Wigner, Central
double V_M_BB[2] = { -443.085 , -44.7881 };         // Majorana, Spatial Permutation
double V_B_BB[2] = { 0. , 0. };         // Barlett, Spin Permutation
double V_H_BB[2] = { 0. , 0. };         // Heisenberg, Isospin Permutation
//*/
/*
// Parameters for hw =21
double V_W_BB[2] = { 595.5455  , -72.2122 };    //{ 595.55 , -72.21 }; //  // Wigner, Central
double V_M_BB[2] = { -206.0455 , -68.3878 };         // Majorana, Spatial Permutation
double V_B_BB[2] = { 0. , 0. };         // Barlett, Spin Permutation
double V_H_BB[2] = { 0. , 0. };         // Heisenberg, Isospin Permutation
//*/
///*
// Inverted expression
double A_BB[2] = { V_W_BB[0] - (0.25*V_M_BB[0]) + (0.5*V_B_BB[0]) + (0.5*V_H_BB[0]),
                     V_W_BB[1] - (0.25*V_M_BB[1]) + (0.5*V_B_BB[1]) +(0.5*V_H_BB[1])};    // Central
double B_BB[2] = { (0.5*V_B_BB[0]) - (0.25*V_M_BB[0]) ,
                    (0.5*V_B_BB[1]) - (0.25*V_M_BB[1]) };            			   // Spin
double C_BB[2] = { (0.5*V_H_BB[0]) - (0.25*V_M_BB[0]),
                     (0.5*V_H_BB[1]) - (0.25*V_M_BB[1]) };            			   // Isospin
double D_BB[2] = { -0.25*V_M_BB[0] , -0.25*V_M_BB[1] };                			   // Spin * Isospin


//*/

//*/////////    ABCD Scheme conversion to Permutation Constants  //////////
/*double A_BB[2] = { 0.375 , 0.375 };     // Central
double B_BB[2] = { -.125 , -.125 };     // Spin
double C_BB[2] = { -.125 , -.125 };     // Isospin
double D_BB[2] = { -.125 , -.125 };     // Spin * Isospin
*/
/*
double A_BB[2] = { -1 , 0};     // Central
double B_BB[2] = { 0 , 0 };     // Spin
double C_BB[2] = { 0 , 0 };     // Isospin
double D_BB[2] = { 0 , 0 };     // Spin * Isospin
//double C_BB[2] = { -3.541 , -3.541};     // Isospin
//double D_BB[2] = { -8.25053 , -8.25053 };     // Spin * Isospin

// Inverted expression (sign change in V_H)
double V_W_BB[2] = { A_BB[0] - B_BB[0] - C_BB[0] + D_BB[0] ,
		      A_BB[1] - B_BB[1] - C_BB[1] + D_BB[1] };    // Wigner, Central
double V_M_BB[2] = { -4*D_BB[0] , -4*D_BB[1] };                     // Majorana, Spatial Permutation
double V_B_BB[2] = { (2*B_BB[0]) - (2*D_BB[0]) ,
		      (2*B_BB[1]) - (2*D_BB[1]) };                // Barlett, Spin Permutation
double V_H_BB[2] = { (2*C_BB[0]) - (2*D_BB[0]) ,
		    (2*C_BB[1]) - (2*D_BB[1]) };                // Heisenberg, Isospin Permutation
//*/

/////////////////////////////////////////////////////////////////////////////
///////////////                FUNCTIONS                /////////////////////
/////////////////////////////////////////////////////////////////////////////


std::complex<double> conversion_NLM_XYZ(const QN_1body_radial & spherical, int n_x,int n_y,int n_z){

      // All definitions of the conversion coefficients are with n>=1 in the spherical QN
    int n = spherical.n + 1;
    int l = spherical.l;
    int m = spherical.m_l;

    // Conservation of the energy
    int n_axial = n_x + n_y;

    // If we impose the condition that the combinatorials of the transformations have integer arguments,
    // the following condition has to be fulfilled (necessary condition)
    std::complex <double> Result (0., 0.);
    if( ((n_axial + m) % 2 != 0) || ((n_axial - abs(m)) % 2 != 0) || (n_axial < abs(m)) ){ //
	return Result;
    }
    // The Conservation law imposes the correspondence of the sum between reference frames
    if( (2*(n-1) + l) != (n_axial + n_z)){
	std::cout << "		Not preserves E" << std::endl;
	return Result;
    }

    ///////////// Auxiliary Real Coefficients ///////////////////

    double C_coeff, A_coeff;

    A_coeff = exp( 0.5 * ( factorial(n-1) + double_factorial(2*(n + l) - 1) - ((n + l + 1)* log(2)) ));

    int *num_1;
    num_1 = new int[3];
    int *den_1;
    den_1 = new int[2];

    num_1[0] = l + m;
    num_1[1] = l - m;
    num_1[2] = n_z;

    den_1[0] = l;
    den_1[1] = l;

    C_coeff = pow(-1, (n_axial + abs(m)) / 2) *
		exp(0.5 * (log(2*l + 1) + factorial_function(num_1,den_1,3,2) - (n_z*log(2)) ))  / (2 * A_coeff);

     // Note:  pow(-1) will always be real and the B_coeff desappears with when the two conversions are used


    ///////////// Conversion Spherical to Axial /////////////////

    double Aux_axial2spheric = 0.;

    //   The shorthand using of the rest  of  the division (%) avoid the multiple
    // conditions to maintain the limits and combinatorial arguments as integers.
    int minimum = max( 0, ( (n_axial - l) + ((n_axial - l)%2) )/2 );
    int maximum = min(min(min((n - 1) ,((n_axial - m)/2)), ((n_axial + m ) / 2)),
			   ( (n_axial + (n_axial%2))/2) );

    int *num_2;
    num_2 = new int[3];
    int *den_2;
    den_2 = new int[6];

    // The restriction of integer argument in the factorial imply alpha to be integer
    //std::cout <<"Loop1 start (min,max)="<<minimum<<" , "<<maximum <<std::endl;
    for( int alpha = minimum; alpha <= maximum; alpha++){

        num_2[0] = n - 1;
        num_2[1] = l;
        num_2[2] = n_axial - 2*alpha;

        den_2[0] = alpha;
        den_2[1] = num_2[0] - den_2[0];
        den_2[2] = n_axial - 2*alpha;
        den_2[3] = num_2[1] - den_2[2];
        den_2[4] = (n_axial + m - 2*alpha) / 2;   // its integer since we have not returned 0. (n_axial + m is even)
        den_2[5] = num_2[2] - den_2[4];

        Aux_axial2spheric += pow(-1,(alpha + n_axial - 2*alpha)) *
                                exp(factorial_function(num_2,den_2,3,6))
                                 / pow(2,(n_axial - 2*alpha));
    }
    Aux_axial2spheric *= C_coeff;

    ////////////// Conversion Axial to Cartesian ////////////////////////////

    std::complex <double> Aux_cartesian2axial (0., 0.);
    double Aux = 0.;

    int *num_3;
    num_3 = new int[2];
    int *den_3;
    den_3 = new int[4];

    std::complex< double> unit_i (0,1);
    std::complex< double> fase;

    // If m is negative, the Conversion coeffient need to be conjugated.
    bool Conjugate = false;
    if( m < 0){
	m = (-1) * m;
	Conjugate = true;
    }

    minimum = max(m + n_x - n_axial , 0);
    maximum = min(m ,n_x);

    //std::cout <<"Loop2 start (min,max)="<<minimum<<" , "<<maximum <<std::endl;
    for( int beta = minimum; beta <= maximum; beta++){
	// n_axial - |m| is even due to the previous condition
	if((n_x - beta) % 2 == 0){

	    fase = pow(unit_i, m - beta);

	    num_3[0] = m;
	    num_3[1] = int(n_axial - m) /2 ;

	    den_3[0] = beta;
	    den_3[1] = num_3[0] - den_3[0];
	    den_3[2] = (n_x - beta) /2 ;
	    den_3[3] = num_3[1] - den_3[2];

	    Aux = exp(factorial_function(num_3,den_3,2,4));
	    Aux_cartesian2axial += Aux * fase ;
	}
    }

    if(Conjugate){
	Aux_cartesian2axial = conj(Aux_cartesian2axial);
    }

        // The constant factor
    num_1[0] = n_x;
    num_1[1] = n_y;

    Aux_cartesian2axial *= exp( 0.5 * (factorial_function(num_1,den_1,2,0) -  (log(2)*n_axial) ));

    delete[] num_1;
    delete[] den_1;
    delete[] num_2;
    delete[] den_2;
    delete[] num_3;
    delete[] den_3;

    return  Aux_cartesian2axial * Aux_axial2spheric;
}



// /*
double Coeff_Talman(int n_a, int n_b, int n_mu){
	  // Robledo Talman
    // Geometric coefficient to express the product of two Tridimensional wave functions (TD-SHO)
    // as a linear combination of one TD-SHO.   (tau)
    
    // parity condition
    if( ((n_a + n_b)%2) != (n_mu%2)){return 0.;}
  
    int s = (n_a + n_b + n_mu);
    if( s%2 == 0 ){
        s = s/2;

        int *numerator;
        numerator = new int[3];
        int *denominator;
        denominator = new int[0];

        numerator[0] = n_a;
        numerator[1] = n_b;
        numerator[2] = n_mu;

        double aux;
        aux = 0.5 * factorial_function(numerator, denominator,3,0);  // log value

        numerator[0] = s - n_a;
        numerator[1] = s - n_b;
        numerator[2] = s - n_mu;

        aux -= factorial_function(numerator, denominator,3,0);    // log value

        delete[] numerator;
        delete[] denominator;
	//std::cout << " << Result  Talman = " << exp(aux) << std::endl;
        return exp(aux);
    }
    //std::cout << " Talman = 0  s= H.I" << std::endl;
    return 0;
}


double J_integral(int n_b, int n_d, int n_mu, double mu, double b_param){
    
    // 	Robledo's appendix method for the J integral using Fourier Transform

    int s = (n_d + n_b + n_mu);
    if( s%2 == 0 ){
        s = s/2;

        //double b = b_length(50);

        double G = 1 + (0.5 * pow((mu / b_param),2)) ;
        double Result = pow((-1), s + n_mu) * (mu / (b_param * sqrt(2 * M_PI) * pow( G, s+0.5)) );

        int *numerator;
        numerator= new int[3];
        int *denominator;
        denominator = new int[1];
        double aux;

        numerator[0] = n_d;
        numerator[1] = n_b;
        numerator[2] = 1;

        denominator[0] = n_mu;

        aux = exp(0.5 * factorial_function( numerator, denominator, 2, 1));

        double Sum = 0.;
        int sign = NULL;
        for (int p = 0; p <= min(n_b, n_d); p++){

            numerator[0] = n_d - p;
            numerator[1] = n_b - p;
            numerator[2] = p;

            denominator[0] = 1;

            Sum += exp(gamma_function( s - p + Fraction(1,2), sign) - factorial_function(numerator, denominator,3,0)  )
                        * pow( (-1)*G, p);
                        // no need for the sign, p is always lower than s.
        }

        Result *= aux * Sum;

        delete[] numerator;
        delete[] denominator;
	//std::cout << " << Result  J = " << Result << std::endl;
        return Result;
    }
    //std::cout << " J_integral = 0  s= H.I" << std::endl;
    
    return 0;
}
// */

 /*
double Coeff_Talman(int n_a, int n_b, int n_mu){
	  // Girod I_ac mu
  
    // Geometric coefficient to express the product of two Tridimensional wave functions (TD-SHO)
    // as a linear combination of one TD-SHO.   (tau)
    
    // parity condition
    if( ((n_a + n_b)%2) != (n_mu%2)){return 0.;}
  
    int s = (n_a + n_b + n_mu);
    if( s%2 == 0 ){
        s = s/2;

        int *numerator;
        numerator = new int[3];
        int *denominator;
        denominator = new int[1];

        numerator[0] = n_a;
        numerator[1] = n_b;
        numerator[2] = 0;

        double aux;
        aux = 0.5 * ((n_mu * log(2)) + factorial_function(numerator, denominator,2,0));  // log value
	
	denominator[0] = n_mu;
        numerator[0] = s - n_a;
        numerator[1] = s - n_b;
        numerator[2] = s - n_mu;

        aux -= factorial_function(numerator, denominator,3,1);    // log value

        delete[] numerator;
        delete[] denominator;
	//std::cout << " << Result  Talman = " << exp(aux) << std::endl;
        return exp(aux);
    }
    //std::cout << " Talman = 0  s= H.I" << std::endl;
    return 0;
}


double J_integral(int n_b, int n_d, int n_mu, double mu_param, double b_param){
     
    // PHYSICAL REVIEW C VOLUME 27, NUMBER 5 MAY 1983. (M. Girod & B.Grammaticos');
  
    double Result = 0.;
    
    double G = 2 + pow(mu_param / b_param,2);
    
    double Moshinsky_0, Integral;
    
    int *num;
    int *den;
    num = new int[1];
    den = new int[2];
    double fact_aux, Constant;
    
    int k;
    for (int n_nu = abs(n_b - n_d); n_nu <= (n_b + n_d);n_nu++){
	
	// Calculation of the integral, which is cero for odd n_mu + n_nu 
	if( (n_mu + n_nu)%2 == 0 ){
	      // Calculation of the integral,   
	    k = (n_mu + n_nu)/2;
	    
	    num[0] = 2*k;
	    den[0] = k;
	    den[1] = k; 
	    
	    fact_aux = exp(0.5 * factorial_function(num,den,1,2));
	    
	    Integral = mu_param * sqrt( sqrt(M_PI) / (b_param * G)) * fact_aux * pow((-1 / G) ,k);
	  
	      // Calculation of the talman coefficient times constants
	    Constant = Coeff_Talman(n_b,n_d,n_nu) / (sqrt( pow(2,n_nu) * exp(factorial(n_nu)) ));
	    
	      // Calculation of the 1D Moshinsky Coefficient M_0^{n_mu,n_nu}
	    num[0] = n_mu + n_nu;
	    den[0] = n_mu;
	    den[1] = n_nu;
	    
	    fact_aux = exp(factorial_function(num,den,1,2));
	    Moshinsky_0 = pow(-1,n_mu) * sqrt(fact_aux / pow(2,n_mu + n_nu));
	    
	    // Due the condition, only odd cases contribute, then
	    Result += Constant * Moshinsky_0 * Integral;
	}
    }
    // Due the differences with the Robledo's Appendix, it's necesary 
    // to add a factor for the 1D matrix element.

    return sqrt(1. / ( b_param * sqrt(M_PI) * pow(2, n_mu) * exp(factorial(n_mu)) )) * Result;
    
}
// */

double Radial_BB(const QN_1body_radial & a, const QN_1body_radial & b,
                  const QN_1body_radial & c, const QN_1body_radial & d , double mu_param, double b_param){
    
    //std::cout << "Talman begin: "<< std::endl;

    /** Two steps process:
    *       * Transform from radial n,l,m QN to n_x n_y n_z
    *       * Evaluate the product of the three sub-matrix elements in x,y & z
    *           - Sub_M_E are a sum of Talman coeff and the J integral
    */
    //std::cout << "RADIAL FUNCTION REACHED" <<  std::endl;
    //print_QN_1body_rad(a);

    // Conversion, conservation of the energy (top equals to the major quantum oscillator shell number N)
    int top_a = 2*a.n + a.l;
    int top_b = 2*b.n + b.l;
    int top_c = 2*c.n + c.l;
    int top_d = 2*d.n + d.l;


    double Aux_I_x = 0.;
    double Aux_I_y = 0.;
    double Aux_I_z = 0.;

    std::complex< double > Aux_Conversion, aux_coeff_a, aux_coeff_b, aux_coeff_c, aux_coeff_d;
    std::complex< double > Result (0. ,0.);

    // In the loop, all physical parameters are summarized in this quantity. It's Computed here to save time.
    // 		It consists in the product of the three (dimensional) factors asocied to
    // the sum (Talman-J_integral) multiplied by the Normalization factor of the J_integral.
    double Physical_Constant =  1.; //pow((b_param * sqrt(M_PI) ),(-3.));
    //std::cout << "Physical_Constant" << Physical_Constant << std::endl;

    //std::cout << "Top_a,b,c,d =  " <<top_a <<" , "<<top_b <<" , "<<top_c <<" , "<<top_d <<std::endl;
    int  n_a_z, n_b_z, n_c_z, n_d_z;

    for(int n_a_x = 0; n_a_x <= top_a; n_a_x++ ){
        for(int n_a_y = 0; n_a_y <= (top_a - n_a_x); n_a_y++ ){
	  //for(int n_a_z = n_a_y; n_a_z <= (top_a - n_a_x - n_a_y); n_a_z++ ){

	    n_a_z = top_a - n_a_y - n_a_x;
	    //std::cout << "(a) (n_x,n_y,n_z) = "<< n_a_x<<" , "<<n_a_y<<" , "<<n_a_z << std::endl;

	    aux_coeff_a.real(0.);
	    aux_coeff_a.imag(0.);

	    if( (top_a) == (n_a_x + n_a_y + n_a_z)){
		aux_coeff_a = conj(conversion_NLM_XYZ( a , n_a_x, n_a_y, n_a_z));
		//aux_coeff_a = conversion_NLM_XYZ( a , n_a_x, n_a_y, n_a_z);
	    }

	    if( abs(aux_coeff_a) > 1.e-10){

		for(int n_b_x = 0; n_b_x <= top_b; n_b_x++ ){
		    for(int n_b_y = 0; n_b_y <= (top_b - n_b_x); n_b_y++ ){
		    //for(int n_b_z = n_b_y; n_b_z <= (top_b - n_b_x - n_b_y); n_b_z++ ){

			n_b_z = top_b - n_b_y - n_b_x;
			//std::cout << "(b) (n_x,n_y,n_z) = "<< n_b_x<<" , "<<n_b_y<<" , "<<n_b_z << std::endl;

			aux_coeff_b.real(0.);
			aux_coeff_b.imag(0.);

			if( (top_b) == (n_b_x + n_b_y + n_b_z)){
			    aux_coeff_b = conj( conversion_NLM_XYZ( b , n_b_x, n_b_y, n_b_z));
			    //aux_coeff_b = conversion_NLM_XYZ( b , n_b_x, n_b_y, n_b_z);
			}

			if( abs(aux_coeff_b) > 1.e-10){

			    for(int n_c_x = 0; n_c_x <= top_c; n_c_x++ ){
				for(int n_c_y = 0; n_c_y <= (top_c - n_c_x); n_c_y++ ){
				//for(int n_c_z = n_c_y; n_c_z <= (top_c - n_c_x - n_c_y); n_c_z++ ){

				    n_c_z =  top_c - n_c_y - n_c_x;
				    //std::cout << "(c) (n_x,n_y,n_z) = "<< n_c_x<<" , "<<n_c_y<<" , "<<n_c_z << std::endl;

				    aux_coeff_c.real(0.);
				    aux_coeff_c.imag(0.);

				    if( (top_c) == (n_c_x + n_c_y + n_c_z)){
					aux_coeff_c = conversion_NLM_XYZ( c , n_c_x, n_c_y, n_c_z);
				    }

				    if( abs(aux_coeff_c) > 1.e-10){

					for(int n_d_x = 0; n_d_x <= top_d; n_d_x++ ){
					    for(int n_d_y = 0; n_d_y <= (top_d - n_d_x); n_d_y++ ){
					    //for(int n_d_z = n_d_y; n_d_z <= (top_d - n_d_x - n_d_y); n_d_z++ ){

						n_d_z = top_d - n_d_y - n_d_x;

						//std::cout << "(d) (n_x,n_y,n_z) = "<< n_d_x<<" , "<<n_d_y<<" , "<<n_d_z << std::endl;

						aux_coeff_d.real(0.);
						aux_coeff_d.imag(0.);

						if( (top_d) == (n_d_x + n_d_y + n_d_z)){
						    aux_coeff_d = conversion_NLM_XYZ( d , n_d_x, n_d_y, n_d_z);
						}

						if( abs(aux_coeff_d) > 1.e-10){

						    Aux_Conversion = ((aux_coeff_a * aux_coeff_b) * aux_coeff_c) * aux_coeff_d;
						    /*
						    std::cout << std::endl;
						    std::cout <<" *** J integral step reached, Aux_Conversion= "<< Aux_Conversion << std::endl;

						    std::cout << std::endl;
						    std::cout << "aux_coeff_a= "<< aux_coeff_a <<"  ( "<< n_a_x<<" , "<<n_a_y<<" , "<<n_a_z <<" )"<< std::endl;
						    std::cout << "aux_coeff_b= "<< aux_coeff_b <<"  ( "<< n_b_x<<" , "<<n_b_y<<" , "<<n_b_z <<" )"<< std::endl;
						    std::cout << "aux_coeff_c= "<< aux_coeff_c <<"  ( "<< n_c_x<<" , "<<n_c_y<<" , "<<n_c_z <<" )"<< std::endl;
						    std::cout << "aux_coeff_d= "<< aux_coeff_d <<"  ( "<< n_d_x<<" , "<<n_d_y<<" , "<<n_d_z <<" )"<< std::endl;
						    */

						    Aux_I_x = 0.;
						    Aux_I_y = 0.;
						    Aux_I_z = 0.;

						      // Product X,Y,Z [0,1,2] as the positions in the QN arrays ( A class QN is not necessary)

							// Sum for each Sub - M_E

						    for(int n_mu = abs(n_a_x - n_c_x); n_mu <= (n_a_x + n_c_x) ; n_mu++){
							//std::cout << "(X) n_mu " << n_mu << std::endl;
							Aux_I_x += Coeff_Talman(n_a_x, n_c_x, n_mu) * J_integral(n_b_x, n_d_x, n_mu, mu_param,  b_param);
						    }
						    for(int n_mu = abs(n_a_y - n_c_y); n_mu <= (n_a_y + n_c_y) ; n_mu++){
							//std::cout << "(Y) n_mu " << n_mu << std::endl;
							Aux_I_y += Coeff_Talman(n_a_y, n_c_y, n_mu) * J_integral(n_b_y, n_d_y, n_mu, mu_param,  b_param);
						    }
						    for(int n_mu = abs(n_a_z - n_c_z); n_mu <= (n_a_z + n_c_z) ; n_mu++){
							//std::cout << "(Z) n_mu " << n_mu << std::endl;
							Aux_I_z += Coeff_Talman(n_a_z, n_c_z, n_mu) * J_integral(n_b_z, n_d_z, n_mu, mu_param,  b_param);
						    }

						    //std::cout << "Talman * J integtal=" <<  Aux_I_x * Aux_I_y * Aux_I_z <<  std::endl;

						    Result = Result + ((Aux_I_x * Aux_I_y * Aux_I_z)* Physical_Constant * Aux_Conversion);

						    //std::cout << " Result = " << Result << std::endl;
						    //std::cout << std::endl;


						}//end if d

					    //}//z
					    }
					}//  end of d
				    }// end if c

				//} //z
				}
			    }//  end of c
			}// end if b

		    //}//z
		    }
		} //  end of b
	    }// end if a
	//}z
        }
    } //  end of a

    //std::cout << "Result (RADIAL)" << Result << std::endl;

    // Non-zero Imaginary part Warning
    bool continua;
    if(fabs(Result.imag()) > 1.e-10){
	std::cout << ""<<std::endl;
	std::cout << "Result (RADIAL)=" << Result << std::endl;
	std::cin>>continua;
    }
    //std::cout << "Talman end: "<< std::endl;
    return real(Result);  	// Conversion to double
}

double Radial_BB_Moshinsky(const QN_1body_radial & a, const QN_1body_radial & b,
                  const QN_1body_radial & c, const QN_1body_radial & d , double mu_param, double b_param){
    //std::cout << "Moshinsky begin: "<< std::endl;
    /** In this variant of the previous function, we evaluate the central ME by the BMB brackets:
    *       * Evaluate ME
    *       *
    */

    //double nu = M_PROTON * h_bar_omega(A) / pow( H_BAR_C,2);
    Set_Manual_Lambda = true;
    //LAMBDA = 1./(sqrt(nu) * mu_param);
    LAMBDA = b_param/ mu_param;  // b_param is already sqrt(2)*b_param

    double Aux_Conversion, Aux_I, aux_coeff;
    double Result = 0.0;

    int mu,mu_q;
    // lambda ket

    for(int lambda = max(abs(a.l - b.l),abs(c.l - d.l)); lambda <= min((a.l + b.l),(c.l + d.l)); lambda++){

	//for( mu = -lambda; mu <= lambda; mu++){
	    //if( (mu == (a.m_l + b.m_l)) && (mu == (c.m_l + d.m_l))){
	mu = a.m_l + b.m_l;
	if( mu == (c.m_l + d.m_l)){
	    Result += Clebsh_Gordan(a.l, b.l, lambda, a.m_l, b.m_l, mu)*
		      Clebsh_Gordan(c.l, d.l, lambda, c.m_l, d.m_l, mu)*
		      M_E_Central(a.n, a.l, b.n, b.l, lambda, c.n, c.l, d.n, d.l, lambda,0, b_param);

	    //}
	}

    }
    /* identico
    for(int lambda = abs(a.l - b.l); lambda <= (a.l + b.l); lambda++){
        mu = a.m_l + b.m_l;

	if( lambda >= abs(mu)){
	    for(int lambda_q = abs(c.l - d.l); lambda_q <= (c.l + d.l); lambda_q++ ){
		mu_q = c.m_l + d.m_l;

		if( lambda_q >= abs(mu_q)){
		    if((lambda == lambda_q) && (mu == mu_q)){
			Result += Clebsh_Gordan(a.l, b.l, lambda, a.m_l, b.m_l, mu)*
				    Clebsh_Gordan(c.l, d.l, lambda_q, c.m_l, d.m_l, mu_q)*
				      M_E_Central(a.n, a.l, b.n, b.l, lambda, c.n, c.l, d.n, d.l, lambda_q,0, b_param);

		    }
		}

	    }
	}

    }*/
    //std::cout << "Moshinsky end: "<< std::endl;
    return Result;
}

double Non_antisimetriced_BB_Multi(const QN_2body_jj_Coupling & BRA, const QN_2body_jj_Coupling & KET,int lambda_max, double b_param){
    
    // Constants for exchange operators and,
    double *Parameters;
    Parameters = new double[4];
    // Change from basic operators to the exchange.

            // Passing of the parameters
    int n_a = BRA.n1;
    int l_a = BRA.l1;
    Fraction j_a = BRA.j1;
    int n_b = BRA.n2;
    int l_b = BRA.l2;
    Fraction j_b = BRA.j2;
    Fraction J = BRA.J;
    Fraction M = BRA.M;
    int T = BRA.T;

    int n_c = KET.n1;
    int l_c = KET.l1;
    Fraction j_c = KET.j1;
    int n_d = KET.n2;
    int l_d = KET.l2;
    Fraction j_d = KET.j2;
    Fraction J_q = KET.J;
    Fraction M_q = KET.M;
    int T_q = KET.T;

    // Starting of computation
    if(M != M_q){return 0;}
    if(J != J_q){return 0;}
    if(T != T_q){return 0;}

        // Isospin and total coupled angular momentum must be odd it the states are in the same orbit:
    if(((n_a == n_b) && (l_a == l_b) && (j_a == j_b)) || ((n_c == n_d) && (l_c == l_d) && (j_c == j_d))){
	if((J + T)%2 != 1){return 0;}
    }
    
    //double mu_param = 1.; //1.48; // fm
    //double V_0 = -1; //-70.82;

    double Result = 0.;
    double aux_U;
    double rad_aux;
    
    int GL_order = 80;
    run_Gauss_Legendre(GL_order);
    
    // Reduced struct argument for the integrals.
    QN_2body_radial Q_N_radial_BRA;
    Q_N_radial_BRA = {BRA.n1, BRA.l1, BRA.n2, BRA.l2, 0, 0};// Total angular momentum and 3rd component are irrelevant
    QN_2body_radial Q_N_radial_KET;
    Q_N_radial_KET = {KET.n1, KET.l1, KET.n2, KET.l2, 0, 0}; 
    
    // One element for each gaussian. For pure Rosenfeld force set maximum i=0 or A_BB[1] = 0;
    for (int i=0; i<2; i++){ 
	// lambda descomposition
	for (int lambda = 0; lambda <= lambda_max; lambda ++){
	    //std::cout << "i, lambda =" <<i << " , "<< lambda << std::endl;
	    Parameters[0] = A_BB[i];
	    Parameters[1] = B_BB[i];
	    Parameters[2] = C_BB[i];
	    Parameters[3] = D_BB[i];
	    // Computation of the direct Term V_abcd(JT)
	    aux_U = U_Coeff(BRA, KET ,lambda, Parameters);
	    //std::cout <<" // U dir = "<<aux_U<<std::endl;
	    
	    if(fabs(aux_U) > 1.0e-10){

		//rad_aux = Radial_4_order(Q_N_radial_BRA, Q_N_radial_KET, lambda,b_param);
		rad_aux = Radial_4_order_Gauss_Quadrature(Q_N_radial_BRA, Q_N_radial_KET, GL_order ,lambda, mu_BB[i],b_param);
		
		Result +=  rad_aux * aux_U;
		//std::cout <<"Radial 4= "<< rad_aux << std::endl;
	    }

	    //std::cout <<"Direct: R= "<< Radial_4_order(Q_N_radial_BRA, Q_N_radial_KET, lambda) << " // U = "<<U_Coeff(BRA, KET ,lambda)<<std::endl;
	}
    }

    return Result; 
}

double Non_antisimetriced_BB_M_E(const QN_2body_jj_Coupling & Q_Numbers_left, const QN_2body_jj_Coupling & Q_Numbers_right, double b_param){

    /*	This function compute the matrix element of the BB interaction, without consider the antisimetrization.
     * 		0) evaluate isospin part and constants
     *  	1) Uncouple J and jj scheme
     * 		2) recouple to S and  evaluate spin part
     * 		3) evaluate radial part via Moshinsky method or Talman method and sum up
     */

    bool Moshinsky_Method = true;	// Option to use or not this method instead of the Tridimensional decomposition.

        // Passing of the parameters
    int n_a = Q_Numbers_left.n1;
    int l_a = Q_Numbers_left.l1;
    Fraction j_a = Q_Numbers_left.j1;
    int n_b = Q_Numbers_left.n2;
    int l_b = Q_Numbers_left.l2;
    Fraction j_b = Q_Numbers_left.j2;
    Fraction J = Q_Numbers_left.J;
    Fraction M = Q_Numbers_left.M;
    int T = Q_Numbers_left.T;

    int n_c = Q_Numbers_right.n1;
    int l_c = Q_Numbers_right.l1;
    Fraction j_c = Q_Numbers_right.j1;
    int n_d = Q_Numbers_right.n2;
    int l_d = Q_Numbers_right.l2;
    Fraction j_d = Q_Numbers_right.j2;
    Fraction J_q = Q_Numbers_right.J;
    Fraction M_q = Q_Numbers_right.M;
    int T_q = Q_Numbers_right.T;

    // Starting of computation
    if(M != M_q){return 0;}
    if(J != J_q){return 0;}
    if(T != T_q){return 0;}

        // Isospin and total coupled angular momentum must be odd it the states are in the same orbit:
    if(((n_a == n_b) && (l_a == l_b) && (j_a == j_b)) || ((n_c == n_d) && (l_c == l_d) && (j_c == j_d))){
	if((J + T)%2 != 1){return 0;}
    }

    // Despite the Isospin parts, the computation requires for both spin and radial part
    // an uncoupled LS expression. The uncoupling is carried out with the same loop.
    //
    // Evaluation of the Isospin part (already T coupled)
    double Isospin_Radial[2] = {A_BB[0] + (2*T*(T+1) - 3)*C_BB[0], A_BB[1] + (2*T*(T+1) - 3)*C_BB[1]} ;
    double Isospin_Spin[2]   = {B_BB[0] + (2*T*(T+1) - 3)*D_BB[0], B_BB[1] + (2*T*(T+1) - 3)*D_BB[1]} ;
    double radial_value[2]   = { 0. , 0. };

    std::cout << "Isospin_Radial value [0][1]= " << Isospin_Radial[0] << " , " << Isospin_Radial[1] << std::endl;
    std::cout << "Isospin_Spin   value [0][1]= " << Isospin_Spin[0] << " , " << Isospin_Spin[1] << std::endl;


    Fraction m_b, m_d;
    int m_l_a, m_l_b, m_l_c, m_l_d;
    int S;
    double unc_J_total, unc_j, coup_S, spin_value;   	//"unc" from uncoupling
    double Result = 0.;
    QN_1body_radial radial_a, radial_b, radial_c, radial_d;

    for( Fraction m_a = Fraction(-1*j_a.numerator, j_a.denominator); m_a <= j_a ; m_a += 1){
        m_b = M - m_a;
        for( Fraction m_c = Fraction(-1*j_c.numerator, j_c.denominator); m_c <= j_c ; m_c += 1){
            m_d = M - m_c;

            unc_J_total = Clebsh_Gordan(j_a,j_b,J, m_a, m_b, M) * Clebsh_Gordan(j_c,j_d,J, m_c, m_d, M);
	    //std::cout << "Uncoupling J total= " << unc_J_total << std::endl;

            if( 1.e-10 < fabs(unc_J_total)){
                for( Fraction m_s_a = Fraction(1,2); m_s_a >= Fraction(-1,2) ; m_s_a -= 1){
                    m_l_a = m_a - m_s_a;
                    for( Fraction m_s_b = Fraction(1,2); m_s_b >= Fraction(-1,2) ; m_s_b -= 1){
                        m_l_b = m_b - m_s_b;
                        for( Fraction m_s_c = Fraction(1,2); m_s_c >= Fraction(-1,2) ; m_s_c -= 1){
                            m_l_c = m_c - m_s_c;
                            for( Fraction m_s_d = Fraction(1,2); m_s_d >= Fraction(-1,2) ; m_s_d -= 1){
                                m_l_d = m_d - m_s_d;

                                unc_j = Clebsh_Gordan(Fraction(l_a),Fraction(1,2),j_a, Fraction(m_l_a),m_s_a, m_a )
                                        * Clebsh_Gordan(Fraction(l_b),Fraction(1,2),j_b, Fraction(m_l_b),m_s_b, m_b )
                                          * Clebsh_Gordan(Fraction(l_c),Fraction(1,2),j_c, Fraction(m_l_c),m_s_c, m_c )
                                            * Clebsh_Gordan(Fraction(l_d),Fraction(1,2),j_d, Fraction(m_l_d),m_s_d, m_d );
				std::cout << std::endl;
				std::cout << "Uncoupling j LS= " << unc_j << "  // Uncoupling J total= " << unc_J_total << std::endl;

                                if( 1.e-10 < fabs(unc_j)){
                                        // Evaluation of the Spin part (before because it's simpler than the radial one)
                                    spin_value = 0.;
                                    S = 0;
				    /*std::cout << "SPIN PART: "<< std::endl;
				    std::cout << Fraction(1,2)<<" "<<Fraction(1,2)<<" "<<Fraction(S)<<" "<<m_s_a<<" "<<m_s_b
					      <<" "<<Fraction(0) << " = " <<
					      Clebsh_Gordan(Fraction(1,2),Fraction(1,2),Fraction(S), m_s_a, m_s_b, Fraction(0) )<<std::endl;
				    std::cout << Fraction(1,2)<<" "<<Fraction(1,2)<<" "<<Fraction(S)<<" "<<m_s_c<<" "<<m_s_d
					      <<" "<<Fraction(0) << " = "<<
					      Clebsh_Gordan(Fraction(1,2),Fraction(1,2),Fraction(S), m_s_c, m_s_d, Fraction(0) ) <<std::endl;
				    */
                                    spin_value += Clebsh_Gordan(Fraction(1,2),Fraction(1,2),Fraction(S), m_s_a, m_s_b, Fraction(0) )
                                                  * Clebsh_Gordan(Fraction(1,2),Fraction(1,2),Fraction(S), m_s_c, m_s_d, Fraction(0) )
                                                    * (-3);  // <S> = (2*S*(S+1) - 3);

                                    S = 1;
                                    for(int M_S = -S; M_S <= S; M_S++){
					/*
					std::cout << Fraction(1,2)<<" "<<Fraction(1,2)<<" "<<Fraction(S)<<" "<<m_s_a<<" "<<m_s_b
						  <<" "<<Fraction(M_S) << " = " <<
						  Clebsh_Gordan(Fraction(1,2),Fraction(1,2),Fraction(S), m_s_a, m_s_b, Fraction(M_S) )<<std::endl;
					std::cout << Fraction(1,2)<<" "<<Fraction(1,2)<<" "<<Fraction(S)<<" "<<m_s_c<<" "<<m_s_d
						  <<" "<<Fraction(M_S) << " = "<<
						  Clebsh_Gordan(Fraction(1,2),Fraction(1,2),Fraction(S), m_s_c, m_s_d, Fraction(M_S) ) <<std::endl;
					*/
                                        spin_value += Clebsh_Gordan(Fraction(1,2),Fraction(1,2),Fraction(S), m_s_a, m_s_b, Fraction(M_S) )
                                                      * Clebsh_Gordan(Fraction(1,2),Fraction(1,2),Fraction(S), m_s_c, m_s_d, Fraction(M_S) )
                                                        * (1); // <S> = (2*S*(S+1) - 3);
                                    }

                                    //std::cout << std::endl;
                                    std::cout << "Spin value= " << spin_value << std::endl;
				    //std::cout << std::endl;


                                    //  MONSHINSKY PROCEDURE //
				    if(Moshinsky_Method){

					//if(1.e-10 < fabs(spin_value)){

						// If the spin part is non-negative, evaluate the lengthy radial part
						    // Evaluation of the Radial part,
					    radial_a = {n_a, l_a, m_l_a};
					    radial_b = {n_b, l_b, m_l_b};
					    radial_c = {n_c, l_c, m_l_c};
					    radial_d = {n_d, l_d, m_l_d};

					    radial_value[0] = Radial_BB_Moshinsky(radial_a, radial_b, radial_c, radial_d, mu_BB[0], sqrt(2)*b_param);
					    std::cout << " --------"<<std::endl;
					    radial_value[1] = Radial_BB_Moshinsky(radial_a, radial_b, radial_c, radial_d, mu_BB[1], sqrt(2)*b_param);

					//}
				    }
                                    //*/

                                    // PROCEDURE FOR DECOMPOSITION IN nX nY nZ //
                                    else{

					//if(1.e-10 < fabs(spin_value)){
						// If the spin part is non-negative, evaluate the lengthy radial part
						    // Evaluation of the Radial part,
					    radial_a = {n_a, l_a, m_l_a};
					    radial_b = {n_b, l_b, m_l_b};
					    radial_c = {n_c, l_c, m_l_c};
					    radial_d = {n_d, l_d, m_l_d};

					    radial_value[0] = Radial_BB(radial_a, radial_b, radial_c, radial_d, mu_BB[0],b_param);
					    radial_value[1] = Radial_BB(radial_a, radial_b, radial_c, radial_d, mu_BB[1],b_param);

					//}
				    }
                                    std::cout << "Radial value [0][1]= " << radial_value[0] << " , " << radial_value[1] << std::endl;

                                        // Add everything
                                    Result += (unc_J_total * unc_j) *
                                                  ((radial_value[0] * (Isospin_Radial[0] + (spin_value * Isospin_Spin[0])))
                                                   + (radial_value[1] * (Isospin_Radial[1] + (spin_value * Isospin_Spin[1]))));

                                } // end if j

                            } // end for m_s
                        }
                    }
                }
            }// end if J

        }
    }
    std::cout << "Result = "<< Result << std::endl;
    return Result; //Normalization ???
}
//*//

    // DEVANATHAN APPROACH WITH LS-jj COUPLING COEFFICENT
double Non_antisimetriced_BBME_by_exchange_LSjj(const QN_2body_jj_Coupling & Q_Numbers_left,
					   const QN_2body_jj_Coupling & Q_Numbers_right, double b_param){

    /**
     *  This variant of the function  evaluate the BB interaction by the exchange operators, and not by < sigma(1).sigma(2)>/< tau(1).tau(2)>
     * 	in order to simplify the computation and comparation.
     *		1) uncoupling  J and jj scheme
     * 		2) evaluate radial with Moshinsky Method
     * 		3) couple to LS and apply exhange operators
     */
    
    bool Moshinsky_Method = true;

        // Passing of the parameters
    int n_a = Q_Numbers_left.n1;
    int l_a = Q_Numbers_left.l1;
    Fraction j_a = Q_Numbers_left.j1;
    int n_b = Q_Numbers_left.n2;
    int l_b = Q_Numbers_left.l2;
    Fraction j_b = Q_Numbers_left.j2;
    Fraction J = Q_Numbers_left.J;
    Fraction M = Q_Numbers_left.M;
    int T = Q_Numbers_left.T;

    int n_c = Q_Numbers_right.n1;
    int l_c = Q_Numbers_right.l1;
    Fraction j_c = Q_Numbers_right.j1;
    int n_d = Q_Numbers_right.n2;
    int l_d = Q_Numbers_right.l2;
    Fraction j_d = Q_Numbers_right.j2;
    Fraction J_q = Q_Numbers_right.J;
    Fraction M_q = Q_Numbers_right.M;
    int T_q = Q_Numbers_right.T;

    // Starting of computation
    if(M != M_q){return 0;}
    if(J != J_q){return 0;}
    if(T != T_q){return 0;}

    // Isospin and total coupled angular momentum must be odd it the states are in the same orbit:
    if(((n_a == n_b) && (l_a == l_b) && (j_a == j_b)) || ((n_c == n_d) && (l_c == l_d) && (j_c == j_d))){
	if((J + T)%2 != 1){return 0;}
    }

    double radial_value[2]   = { 0. , 0. };
    double exchange_value[2] = { 0. , 0. };

    //double nu = M_PROTON * h_bar_omega(A) / pow( H_BAR_C,2);
    Set_Manual_Lambda = true;
    //LAMBDA = 1./(sqrt(nu) * mu_param);
    
    int GL_order = 150;
    run_Gauss_Legendre(GL_order);
    
    int m_l_a, m_l_b, m_l_c, m_l_d;
    double unc_J_total , unc_l, LSjj_1,LSjj_2;   	//"unc" from uncoupling
    double Result = 0.;
    QN_1body_radial radial_a, radial_b, radial_c, radial_d;

    for( int L = abs(l_a-l_b); L <= l_a + l_b; L++ ){
        for( int S = 0; S <= 1; S++){

            LSjj_1 = LS_jj_coupling_Coeff( Fraction(l_a), Fraction(1,2),j_a,
                            Fraction(l_b),Fraction(1,2),j_b, Fraction(L),  Fraction(S),  J);
            if(1.e-10 < fabs(LSjj_1) ){

                LSjj_2 = LS_jj_coupling_Coeff( Fraction(l_c), Fraction(1,2),j_c,
				Fraction(l_d),Fraction(1,2),j_d, Fraction(L),  Fraction(S),  J);
                unc_J_total = LSjj_1 * LSjj_2;

                //std::cout << "Uncoupling J total= " << unc_J_total << std::endl;

                if( 1.e-10 < fabs(unc_J_total)){

                      //std::cout << std::endl;
                      //std::cout << "Uncoupling j LS= " << unc_j << "  // Uncoupling J total= " << unc_J_total << std::endl;

                          // If the spin part is non-negative, evaluate the lengthy radial part
                          // Evaluation of the Radial part,
		      
		      if(Moshinsky_Method){
			    // suposing that the 
			  LAMBDA = (b_param*sqrt(2.)/ mu_BB[0]);
			  radial_value[0] = M_E_Central(n_a, l_a, n_b, l_b, L, n_c, l_c, n_d, l_d, L,0, b_param*sqrt(2)) * sqrt(2*L+1);
			  LAMBDA = (b_param*sqrt(2.)/ mu_BB[1]);
			  radial_value[1] = M_E_Central(n_a, l_a, n_b, l_b, L, n_c, l_c, n_d, l_d, L,0, b_param*sqrt(2)) * sqrt(2*L+1);
		      }  
		      else{
			  // obtaining a M component by Wigner Echkart, its arbitrary chosen to 0.
			  
			  // uncoupling in l's from L,M_L
			  for(m_l_a = ((-1)*l_a); m_l_a <= l_a; m_l_a++){
			      for(m_l_c = ((-1)*l_c); m_l_c <= l_c; m_l_c++){
				  
				  m_l_b = -m_l_a;
				  m_l_d = -m_l_c;
				  
				  unc_l = sqrt(2*L+1) * Clebsh_Gordan(l_a,l_b,L, m_l_a,m_l_b,0) * Clebsh_Gordan(l_c,l_d,L, m_l_c,m_l_d,0);
				  if(1.e-10 < unc_l ){
				  
				      radial_a = {n_a, l_a, m_l_a};
				      radial_b = {n_b, l_b, m_l_b};
				      radial_c = {n_c, l_c, m_l_c};
				      radial_d = {n_d, l_d, m_l_d};

				      radial_value[0] += unc_l * Radial_BB(radial_a, radial_b, radial_c, radial_d, mu_BB[0],b_param);
				      radial_value[1] += unc_l * Radial_BB(radial_a, radial_b, radial_c, radial_d, mu_BB[1],b_param);
				      
				  }
			      }
			  }
			  
		      }

		      /*
                      LAMBDA = (sqrt(2)*b_param/ mu_BB[0]);
                      radial_value[0] = M_E_Central(n_a, l_a, n_b, l_b, L, n_c, l_c, n_d, l_d, L,0, b_param*sqrt(2));
                      LAMBDA = (sqrt(2)*b_param/ mu_BB[1]);
                      radial_value[1] = M_E_Central(n_a, l_a, n_b, l_b, L, n_c, l_c, n_d, l_d, L,0, b_param*sqrt(2));
                      
                      

                      exchange_value[0] = (V_W_BB[0] + (pow(-1,L)*V_M_BB[0]) + (pow(-1,S)*V_B_BB[0]) + (pow(-1,T)*V_H_BB[0]));
                      exchange_value[1] = (V_W_BB[1] + (pow(-1,L)*V_M_BB[1]) + (pow(-1,S)*V_B_BB[1]) + (pow(-1,T)*V_H_BB[1]));
                      
		      
		      radial_value[0] = Rosenfeld_quadratures(n_a, l_a, n_b, l_b, L, n_c, l_c, n_d, l_d, L ,GL_order, mu_BB[0] , b_param);
                      radial_value[1] = Rosenfeld_quadratures(n_a, l_a, n_b, l_b, L, n_c, l_c, n_d, l_d, L ,GL_order, mu_BB[1] , b_param);
		      
                      exchange_value[0] = (V_W_BB[0] + (pow(-1,L)*V_M_BB[0]) + (pow(-1,S)*V_B_BB[0]) + (pow(-1,T)*V_H_BB[0]));
                      exchange_value[1] = (V_W_BB[1] + (pow(-1,L)*V_M_BB[1]) + (pow(-1,S)*V_B_BB[1]) + (pow(-1,T)*V_H_BB[1]));
                      */
                      // Add everything
                      Result += (unc_J_total) * ( (radial_value[0] * exchange_value[0]) + (radial_value[1] * exchange_value[1]));

                }


            }// end J1

        }
    }

    //std::cout << "Result = "<< Result << std::endl;
    return Result;

}
//*/

///* 	// COUPLING AFTER UNCOUPLING
double Non_antisimetriced_BBME_by_exchange(const QN_2body_jj_Coupling & Q_Numbers_left,
					   const QN_2body_jj_Coupling & Q_Numbers_right, double b_param){

    /**
     *  This variant of the function  evaluate the BB interaction by the exchange operators, and not by < sigma(1).sigma(2)>/< tau(1).tau(2)>
     * 	in order to simplify the computation and comparation.
     *		1) uncoupling  J and jj scheme
     * 		2) evaluate radial with Moshinsky Method
     * 		3) couple to LS and apply exhange operators
     */
    bool Moshinsky_Method = true;	// Option to use or not this method instead of the Tridimensional decomposition.

        // Passing of the parameters
    int n_a = Q_Numbers_left.n1;
    int l_a = Q_Numbers_left.l1;
    Fraction j_a = Q_Numbers_left.j1;
    int n_b = Q_Numbers_left.n2;
    int l_b = Q_Numbers_left.l2;
    Fraction j_b = Q_Numbers_left.j2;
    Fraction J = Q_Numbers_left.J;
    Fraction M = Q_Numbers_left.M;
    int T = Q_Numbers_left.T;

    int n_c = Q_Numbers_right.n1;
    int l_c = Q_Numbers_right.l1;
    Fraction j_c = Q_Numbers_right.j1;
    int n_d = Q_Numbers_right.n2;
    int l_d = Q_Numbers_right.l2;
    Fraction j_d = Q_Numbers_right.j2;
    Fraction J_q = Q_Numbers_right.J;
    Fraction M_q = Q_Numbers_right.M;
    int T_q = Q_Numbers_right.T;

    // Starting of computation
    if(M != M_q){return 0;}
    if(J != J_q){return 0;}
    if(T != T_q){return 0;}

        // Isospin and total coupled angular momentum must be odd it the states are in the same orbit:
    if(((n_a == n_b) && (l_a == l_b) && (j_a == j_b)) || ((n_c == n_d) && (l_c == l_d) && (j_c == j_d))){
	if((J + T)%2 != 1){return 0;}
    }

    double radial_value[2]   = { 0. , 0. };
    double exchange_value[2] = { 0. , 0. };

    Fraction m_b, m_d;
    int m_l_a, m_l_b, m_l_c, m_l_d;
    double unc_J_total, unc_j, coup_S, coup_L;   	//"unc" from uncoupling
    double Result = 0.;
    QN_1body_radial radial_a, radial_b, radial_c, radial_d;

    for( Fraction m_a = Fraction(-1*j_a.numerator, j_a.denominator); m_a <= j_a ; m_a += 1){
        m_b = M - m_a;
        for( Fraction m_c = Fraction(-1*j_c.numerator, j_c.denominator); m_c <= j_c ; m_c += 1){
            m_d = M - m_c;

            unc_J_total = Clebsh_Gordan(j_a,j_b,J, m_a, m_b, M) * Clebsh_Gordan(j_c,j_d,J, m_c, m_d, M);
	    //std::cout << "Uncoupling J total= " << unc_J_total << std::endl;

            if( 1.e-10 < fabs(unc_J_total)){
                for( Fraction m_s_a = Fraction(1,2); m_s_a >= Fraction(-1,2) ; m_s_a -= 1){
                    m_l_a = m_a - m_s_a;
                    for( Fraction m_s_b = Fraction(1,2); m_s_b >= Fraction(-1,2) ; m_s_b -= 1){
                        m_l_b = m_b - m_s_b;
                        for( Fraction m_s_c = Fraction(1,2); m_s_c >= Fraction(-1,2) ; m_s_c -= 1){
                            m_l_c = m_c - m_s_c;
                            for( Fraction m_s_d = Fraction(1,2); m_s_d >= Fraction(-1,2) ; m_s_d -= 1){
                                m_l_d = m_d - m_s_d;

                                unc_j = Clebsh_Gordan(Fraction(l_a),Fraction(1,2),j_a, Fraction(m_l_a),m_s_a, m_a )
                                        * Clebsh_Gordan(Fraction(l_b),Fraction(1,2),j_b, Fraction(m_l_b),m_s_b, m_b )
                                          * Clebsh_Gordan(Fraction(l_c),Fraction(1,2),j_c, Fraction(m_l_c),m_s_c, m_c )
                                            * Clebsh_Gordan(Fraction(l_d),Fraction(1,2),j_d, Fraction(m_l_d),m_s_d, m_d );
				//std::cout << std::endl;
				//std::cout << "Uncoupling j LS= " << unc_j << "  // Uncoupling J total= " << unc_J_total << std::endl;

                                if( 1.e-10 < fabs(unc_j)){

				    // PARTE RADIAL
				    // MONSHINSKY PROCEDURE ///
				    if(Moshinsky_Method){

					//if(1.e-10 < fabs(spin_value)){

						// If the spin part is non-negative, evaluate the lengthy radial part
						    // Evaluation of the Radial part,
					    radial_a = {n_a, l_a, m_l_a};
					    radial_b = {n_b, l_b, m_l_b};
					    radial_c = {n_c, l_c, m_l_c};
					    radial_d = {n_d, l_d, m_l_d};
					    
					    LAMBDA = (b_param/ mu_BB[1]);
					    radial_value[0] = Radial_BB_Moshinsky(radial_a, radial_b, radial_c, radial_d, mu_BB[0], b_param*sqrt(2.));
					    LAMBDA = (b_param/ mu_BB[1]);
					    radial_value[1] = Radial_BB_Moshinsky(radial_a, radial_b, radial_c, radial_d, mu_BB[1], b_param*sqrt(2.));

					//}
				    }

                                    // PROCEDURE FOR DECOMPOSITION IN nX nY nZ //
                                    else{

					//if(1.e-10 < fabs(spin_value)){
						// If the spin part is non-negative, evaluate the lengthy radial part
						    // Evaluation of the Radial part,
					    radial_a = {n_a, l_a, m_l_a};
					    radial_b = {n_b, l_b, m_l_b};
					    radial_c = {n_c, l_c, m_l_c};
					    radial_d = {n_d, l_d, m_l_d};

					    radial_value[0] = Radial_BB(radial_a, radial_b, radial_c, radial_d, mu_BB[0],b_param);
					    std::cout << std::endl;
					    radial_value[1] = Radial_BB(radial_a, radial_b, radial_c, radial_d, mu_BB[1],b_param);

					//}
				    }
                                    //std::cout << "Radial value [0][1]= " << radial_value[0] << " , " << radial_value[1] << std::endl;

				    exchange_value[0] = 0.;
				    exchange_value[1] = 0.;

				    for(int S = 0; S <= 1; S++){
					for(int M_S = -S; M_S <= S; M_S++){
					    coup_S = Clebsh_Gordan(Fraction(1,2),Fraction(1,2),Fraction(S), m_s_a, m_s_b, Fraction(0) )
						      * Clebsh_Gordan(Fraction(1,2),Fraction(1,2),Fraction(S), m_s_c, m_s_d, Fraction(0) );

					    for(int L = max(abs(l_a - l_b), abs(l_c - l_d)) ; L <=  min(l_a - l_b, l_c - l_d) ; L++){
						for(int M_L = -L; M_L <= L; M_L++){
						    coup_L = Clebsh_Gordan(l_a,l_b,L, m_l_a, m_l_b, M_L )
							      * Clebsh_Gordan(l_c,l_d,L, m_l_c, m_l_d, M_L );

						    exchange_value[0] += coup_L* coup_S * (V_W_BB[0] + (pow(-1,L)*V_M_BB[0]) 
												     + (pow(-1,S)*V_B_BB[0]) + (pow(-1,T)*V_H_BB[0]));
						    exchange_value[1] += coup_L* coup_S * (V_W_BB[1] + (pow(-1,L)*V_M_BB[1]) 
												     + (pow(-1,S)*V_B_BB[1]) + (pow(-1,T)*V_H_BB[1]));

						}
					    }

					}
				    }


                                        // Add the two gaussians
                                    Result += (unc_J_total * unc_j) * ( (radial_value[0] * exchange_value[0]) + (radial_value[1] * exchange_value[1]));

                                } // end if j

                            } // end for m_s
                        }
                    }
                }
            }// end if J

        }
    }

    //std::cout << "Result = "<< Result << std::endl;
    return Result;

}

//*/

double Non_antisimetriced_r_square(const QN_2body_jj_Coupling & Q_Numbers_left,
					   const QN_2body_jj_Coupling & Q_Numbers_right, double b_param){
					     
        // Passing of the parameters
    int n_a = Q_Numbers_left.n1;
    int l_a = Q_Numbers_left.l1;
    Fraction j_a = Q_Numbers_left.j1;
    int n_b = Q_Numbers_left.n2;
    int l_b = Q_Numbers_left.l2;
    Fraction j_b = Q_Numbers_left.j2;
    Fraction J = Q_Numbers_left.J;
    Fraction M = Q_Numbers_left.M;
    int T = Q_Numbers_left.T;

    int n_c = Q_Numbers_right.n1;
    int l_c = Q_Numbers_right.l1;
    Fraction j_c = Q_Numbers_right.j1;
    int n_d = Q_Numbers_right.n2;
    int l_d = Q_Numbers_right.l2;
    Fraction j_d = Q_Numbers_right.j2;
    Fraction J_q = Q_Numbers_right.J;
    Fraction M_q = Q_Numbers_right.M;
    int T_q = Q_Numbers_right.T;

    // Starting of computation
    if(M != M_q){return 0;}
    if(J != J_q){return 0;}
    if(T != T_q){return 0;}

    // Isospin and total coupled angular momentum must be odd it the states are in the same orbit:
    if(((n_a == n_b) && (l_a == l_b) && (j_a == j_b)) || ((n_c == n_d) && (l_c == l_d) && (j_c == j_d))){
	if((J + T)%2 != 1){return 0;}
    }

    double radial_value   = 0.;

    //double nu = M_PROTON * h_bar_omega(A) / pow( H_BAR_C,2);
    Set_Manual_Lambda = true;
    //LAMBDA = 1./(sqrt(nu) * mu_param);


    double unc_J_total, LSjj_1,LSjj_2;   	//"unc" from uncoupling
    double Result = 0.;
    
    for( int L = abs(l_a-l_b); L <= l_a + l_b; L++ ){
        for( int S = 0; S <= 1; S++){

            LSjj_1 = LS_jj_coupling_Coeff( Fraction(l_a), Fraction(1,2),j_a,
                            Fraction(l_b),Fraction(1,2),j_b, Fraction(L),  Fraction(S),  J);
            if(1.e-10 < fabs(LSjj_1) ){

                LSjj_2 = LS_jj_coupling_Coeff( Fraction(l_c), Fraction(1,2),j_c,
				Fraction(l_d),Fraction(1,2),j_d, Fraction(L),  Fraction(S),  J);
                unc_J_total = LSjj_1 * LSjj_2;

                //std::cout << "Uncoupling J total= " << unc_J_total << std::endl;

                if( 1.e-10 < fabs(unc_J_total)){

                      //std::cout << std::endl;
                      //std::cout << "Uncoupling j LS= " << unc_j << "  // Uncoupling J total= " << unc_J_total << std::endl;

                          // If the spin part is non-negative, evaluate the lengthy radial part
                          // Evaluation of the Radial part,
			
                      radial_value = M_E_Central(n_a, l_a, n_b, l_b, L, n_c, l_c, n_d, l_d, L,4, b_param*sqrt(2));
                      
                      // Add everything
                      Result += sqrt(2*L + 1) * unc_J_total * radial_value;

                }


            }// end J1

        }
    }

    //std::cout << "Result = "<< Result << std::endl;
    return Result;
}

double M_E_Brink_Boeker(const QN_2body_jj_Coupling & Q_Numbers_left, const QN_2body_jj_Coupling & Q_Numbers_right, int option, double b_param){
    /*
    std::cout << " V_W[0][1]=" << V_W_BB[0] << "  "<< V_W_BB[1] << std::endl;
    std::cout << " V_M[0][1]=" << V_M_BB[0] << "  "<< V_M_BB[1] << std::endl;
    std::cout << " V_B[0][1]=" << V_B_BB[0] << "  "<< V_B_BB[1] << std::endl;
    std::cout << " V_H[0][1]=" << V_H_BB[0] << "  "<< V_H_BB[1] << std::endl;

    std::cout << " A[0][1]=" << A_BB[0] << "  "<< A_BB[1] << std::endl;
    std::cout << " B[0][1]=" << B_BB[0] << "  "<< B_BB[1] << std::endl;
    std::cout << " C[0][1]=" << C_BB[0] << "  "<< C_BB[1] << std::endl;
    std::cout << " D[0][1]=" << D_BB[0] << "  "<< D_BB[1] << std::endl;
    */


    /* * Normalization and antisimetrization of matrix elements is required,
     * 	 in this function, the permutation is carried out and the  conventional matrix elements
     * 	 (Non_antisimetriced__BB_M_E) is called
     *
     * */
    //int n_a = Q_Numbers_left.n1;
    int l_a = Q_Numbers_left.l1;
    Fraction j_a = Q_Numbers_left.j1;
    //int n_b = Q_Numbers_left.n2;
    int l_b = Q_Numbers_left.l2;
    Fraction j_b = Q_Numbers_left.j2;
    
    Fraction J = Q_Numbers_left.J;
    Fraction M = Q_Numbers_left.M;
    int T = Q_Numbers_left.T;

    int n_c = Q_Numbers_right.n1;
    int l_c = Q_Numbers_right.l1;
    Fraction j_c = Q_Numbers_right.j1;
    int n_d = Q_Numbers_right.n2;
    int l_d = Q_Numbers_right.l2;
    Fraction j_d = Q_Numbers_right.j2;
    
    Fraction J_q = Q_Numbers_right.J;
    Fraction M_q = Q_Numbers_right.M;
    int T_q = Q_Numbers_right.T;

    // Starting of computation
    if(M != M_q){return 0;}
    if(J != J_q){return 0;}
    if(T != T_q){return 0;}

        // Isospin and total coupled angular momentum must be odd it the states are in the same orbit:
    if(((Q_Numbers_left.n1 == Q_Numbers_right.n2) && (Q_Numbers_left.l1 == Q_Numbers_left.l2) && (Q_Numbers_left.j1 == Q_Numbers_left.j2)) 
	|| ((n_c == n_d) && (l_c == l_d) && (j_c == j_d))){
	if((J + T)%2 != 1){return 0;}
    }

    // Create the permuted ket.
    QN_2body_jj_Coupling  Q_Numbers_right_permuted = {n_d,l_d,j_d, n_c,l_c,j_c, J,M, T, 0};

    //std::cout << " Norm[1] =" << Normalization_JT(Q_Numbers_left) << "  Norm[2]"<<  Normalization_JT(Q_Numbers_right)
	//	<< " // Pow= " << pow(-1, j_c + j_d + J + T) << std::endl;
    
    //int option = 1;
    int lamda_max;
    switch (option){
	//  using exchange operators
	    // * Complete Explicit Separation
	case (0):
	    std::cout << "Exch. Sep" << std::endl; 
	    return Normalization_JT(Q_Numbers_left) * Normalization_JT(Q_Numbers_right) * sqrt(2*J+1)*
		    (Non_antisimetriced_BBME_by_exchange( Q_Numbers_left, Q_Numbers_right, b_param) + (pow(-1, j_c + j_d + J + T)*
			Non_antisimetriced_BBME_by_exchange( Q_Numbers_left, Q_Numbers_right_permuted, b_param)));
	break;
	    // * using LS-jj coupling coefficient
	case (1):
	    std::cout << "Exch. LSjj" << std::endl; 
	    return Normalization_JT(Q_Numbers_left) * Normalization_JT(Q_Numbers_right) * sqrt(2*J+1)*
		  (Non_antisimetriced_BBME_by_exchange_LSjj( Q_Numbers_left, Q_Numbers_right, b_param) + (pow(-1, j_c + j_d + J + T)*
		      Non_antisimetriced_BBME_by_exchange_LSjj( Q_Numbers_left, Q_Numbers_right_permuted, b_param)));
	break;
	
	//  using < sigma(1).sigma(2)>/< tau(1).tau(2)>
	case (2):
	    std::cout << "averages" << std::endl; 
	    return Normalization_JT(Q_Numbers_left) * Normalization_JT(Q_Numbers_right) * sqrt(2*J+1)*
		  (Non_antisimetriced_BB_M_E( Q_Numbers_left, Q_Numbers_right, b_param) + (pow(-1, j_c + j_d + J + T)*
		      Non_antisimetriced_BB_M_E( Q_Numbers_left, Q_Numbers_right_permuted, b_param)));
	break;
	//  using Multipolar descomposition
	case (3):
	    lamda_max = max(max(((l_a + l_c),int(j_a + j_c)),int(j_a + j_b)), int(j_c+j_d)); // +4 is a tip
	    std::cout << "Multi   lambda= "<< lamda_max << std::endl; 
 	    return Normalization_JT(Q_Numbers_left) * Normalization_JT(Q_Numbers_right) * sqrt(2*J+1)*
		  (Non_antisimetriced_BB_Multi( Q_Numbers_left, Q_Numbers_right,lamda_max, b_param) + (pow(-1, j_c + j_d + J + T)*
		      Non_antisimetriced_BB_Multi( Q_Numbers_left, Q_Numbers_right_permuted,lamda_max, b_param)));
	break;	  
	//  r square matrix elements
	case (4):
	    std::cout << "-r^2" << std::endl; 
 	    return Normalization_JT(Q_Numbers_left) * Normalization_JT(Q_Numbers_right) * sqrt(2*J+1)*
		  (Non_antisimetriced_r_square( Q_Numbers_left, Q_Numbers_right, b_param) + (pow(-1, j_c + j_d + J + T)*
		      Non_antisimetriced_r_square( Q_Numbers_left, Q_Numbers_right_permuted, b_param)));
	default: 
	    return 0;
    }
    
    return 0;

}

