#include "../include/factorials.h"
#include "../include/Fractions.h"
#include "../include/Angular_Functions.h"
#include "../include/BM_Brackets.h"
#include "../include/Wave_Functions.h"
#include "../include/Integrals.h"

#include "../include/Index_Coefficients.h"
#define _USE_MATH_DEFINES
#include <cmath>
#include <iostream>

using namespace std;

///  To obtain special coefficients needed in the calculus  of BMB
///  and functions to select the interval of  indexes

///////    CONSTANTS    //////////
//////////////////////////////////

double H_BAR = 1.; // 6.58211928e-22; // MeV s

double H_BAR_C = 197.32697; // MeV fm

double M_PROTON = 940. ;// 938.272046; //MeV c^-2

int A  = 4 ;


// Oscillator Energy (MeV)
double h_bar_omega (int A){
    //return 41 * pow(A,-1./3);
    return ((45*pow(A,-1./3)) - (25*pow(A,-2./3)));
}


// Oscillator length (fm)
double b_lenght(int A){
    return 197.33/(sqrt(M_PROTON * h_bar_omega(A)));
    // sqrt(pow(H_BAR_C,2)/ h_bar_omega(A) * M_PROTON );
    //return 1.005* pow(A,1./6);
}

// Estimation of nuclei radius
double R_radius( int A){
    return 1.25 * pow(A,1./3);
}

// Well depth (MeV)
double V_0 = 30.;

// Radius for Rosenfeld interaction (fm)
double R_interaction = 1.48;

// Auxiliary setter for the BB ME, every other method to explicitly change
// the physical values of the Talmi Integral without adding another useless
// argument to all the functions.
bool Set_Manual_Lambda = false;      // 1º Set it as true
double LAMBDA = 1.0;                 // 2º Default value, should be changed

/////////////// BMB METHOD  /////////////////

// Constant coefficient in BMB_00, for extend factorial evaluations
double A_coefficient(int l1, int l, int l2, int L, int x){
    /// used in BMB_00 computation
    int c1,c2,c3,c4,c5;
    // q bigger or equal than (and to x and 0)
    c1 = l1 - l;
    c2 = l2 - L;
    c3 = - x - 1;
    // q lower or equal to
    c4 = l + l1;
    c5 = L + l2;

    int major = min(c4,c5);
    int minor = max(max(max(max(c1,c2),c3),x),0);

    int *Num;
    Num = new int[2];
    int *Den;
    Den = new int[6];
    double sum = 0.;

    for(int q=minor; q<=major;q++){
        if((l + q - l1)%2!=0){}//do nothing
        else{
            Num[0] = l + q - l1;
            Num[1] = L + q - l2;

            Den[0] = (l + q - l1)/2;
            Den[1] = (l + l1 - q)/2;
            Den[2] = q - x;
            Den[3] = q + x + 1;
            Den[4] = (L + q - l2)/2;
            Den[5] = (L + l2 - q)/2;

            if(((l + q - l1)%2!=0) || ((l + l1 - q)%2!=0) || ((L + q - l2)%2!=0) || ((L + l2 - q)%2!=0)){
                std::cout<< " Hay factores que no son pares"<<std::endl;
            }
            sum += pow(-1,(l + q - l1)/2) * exp(factorial_function(Num,Den,2,6));
        }
    }
    // I reuse the pointer array, due to its equal dimension, to save memory during execution
    Den[0] = l1 + l + x + 1;
    Den[1] = l1 + l - x;
    Den[2] = l1 + x - l;
    Den[3] = l2 + L + x + 1;
    Den[4] = l2 + L - x;
    Den[5] = l2 + x - L;

    Num[0] = l + x - l1;
    Num[1] = L + x - l2;

    double aux = factorial_function(Den,Num,6,2);
    delete[] Num;
    delete[] Den;

    return exp(0.5 * aux) * sum;
}

// Recurrence Matrix element for evaluation of general BMB
double ME_ri2_times_BMB(int i ,int lambda,int n,int l,int N,int L,int n_q,int l_q,int N_q,int L_q, int n1, int l1, int n2, int l2){

    if ((n_q<0)||(l_q<0)||(N_q<0)||(L_q<0)){
        //cout<<"un indice negativo"<<endl;
        return 0;}

    // Matrix element -(r_1)^2 used in recurrence relations for monshinsky brackets
    // The matrix element vanish if n'l'N'L' are not this six cases:
    if(n_q == n-1){
        if(l_q == l){
            if(N_q==N){
                if(L_q==L){
                    //std::cout<<"0)"<<std::endl;
                    // does not depend on i
                    return 0.5 * sqrt(n*(n + l + 0.5))* BMB(n_q,l_q,N_q,L_q,lambda, n1, l1, n2, l2);
                }
                else{return 0;}
            }
            else{return 0;}
        }
        else if(l_q == l+1){
            if(N_q==N-1){
                if(L_q==L+1){
                    if(i == 1){
                        //std::cout<<"2)"<<std::endl;
                        return sqrt(n*N*(l+1)*(L+1)) * pow(-1,(lambda + L + l)) * Racah_Coefficient(l,l+1,L,L+1,1,lambda)* BMB(n_q,l_q,N_q,L_q,lambda, n1, l1, n2, l2);

                    }
                    else{// if(i == 2)
                        //std::cout<<"2)"<<std::endl;
                        return sqrt(n*N*(l+1)*(L+1)) * pow(-1,(lambda + L + l + 1)) * Racah_Coefficient(l,l+1,L,L+1,1,lambda)* BMB(n_q,l_q,N_q,L_q,lambda, n1, l1, n2, l2);

                    }

                }
                else {return 0;}
            }
            else if(N_q==N){
                if(L_q==L-1){
                    if(i == 1){
                        //std::cout<<"3)"<<std::endl;
                        return sqrt(n*(N + L + 0.5)*(l+1)*L) * pow(-1,(lambda + L + l)) * Racah_Coefficient(l,l+1,L,L-1,1,lambda)* BMB(n_q,l_q,N_q,L_q,lambda, n1, l1, n2, l2);

                    }
                    else{// i=2
                        //std::cout<<"3)"<<std::endl;
                        return sqrt(n*(N + L + 0.5)*(l+1)*L) * pow(-1,(lambda + L + l + 1)) * Racah_Coefficient(l,l+1,L,L-1,1,lambda)* BMB(n_q,l_q,N_q,L_q,lambda, n1, l1, n2, l2);

                    }
                }
                else{return 0;}
            }
            else{return 0;}
        }
        else{return 0;}
    }

    else if(n_q == n){
        if(l_q == l){
            if(N_q==N-1){
                if(L_q==L){
                    //std::cout<<"1)"<<std::endl;
                    // does not depend on i
                    return 0.5 * sqrt(N*(N + L + 0.5))* BMB(n_q,l_q,N_q,L_q,lambda, n1, l1, n2, l2);
                }
                else{return 0;}
            }
            else{return 0;}
        }
        else if(l_q == l-1){
            if(N_q==N-1){
                if(L_q == L+1){
                    if(i == 1){
                        //std::cout<<"4)"<<std::endl;
                        return sqrt((n + l + 0.5)*N*l*(L + 1)) * pow(-1,(lambda + L + l)) * Racah_Coefficient(l,l-1,L,L+1,1,lambda)* BMB(n_q,l_q,N_q,L_q,lambda, n1, l1, n2, l2);

                    }
                    else{// i=2
                        //std::cout<<"4)"<<std::endl;
                        return sqrt((n + l + 0.5)*N*l*(L + 1)) * pow(-1,(lambda + L + l + 1)) * Racah_Coefficient(l,l-1,L,L+1,1,lambda)* BMB(n_q,l_q,N_q,L_q,lambda, n1, l1, n2, l2);

                    }
                }
                else { return 0;}
            }
            else if(N_q==N){
                if(L_q == L-1){
                    if(i == 1){
                        //std::cout<<"5)"<<std::endl;
                        return sqrt((n + l + 0.5)*(N + L + 0.5)*l*L) * pow(-1,(lambda + L + l)) * Racah_Coefficient(l,l-1,L,L-1,1,lambda)* BMB(n_q,l_q,N_q,L_q,lambda, n1, l1, n2, l2);

                    }
                    else{
                        //std::cout<<"5)"<<std::endl;
                        return sqrt((n + l + 0.5)*(N + L + 0.5)*l*L) * pow(-1,(lambda + L + l + 1)) * Racah_Coefficient(l,l-1,L,L-1,1,lambda)* BMB(n_q,l_q,N_q,L_q,lambda, n1, l1, n2, l2);

                    }
                }
                else {return 0;}
            }
            else {return 0;}
        }
        else {return 0;}
    }
    else{return 0;}

}

// Coefficient for the decomposition of the
// reduced M.Element < nl||V(r)||n'l'> in Talmi integrals (by ME)
/*//
double B_coefficient(int n, int l, int n_q, int l_q, Fraction p_argument, double b_param){

    // Reevaluation (by me) of the coefficient from Monsinsky matrix element implementation.
    /// coefficient used in the sum decomposition  in Talmi integrals of the reduced matrix element of the potential V(r)
    if((n<0)||(n_q<0)||(l<0)||(l_q<0) ){return 0;}
    //cout << "B(n,l,n',l',p)" << n<< " "<< l<<" "<<n_q<<" "<<l_q<<" "<<p<<endl;

    // Parity conservation
    if(((l+l_q)%2)!=0){return 0;}

    int p = int(p_argument);

    double Sum = 0.;

    Fraction * Den;
    Fraction * Num;
    Den = new Fraction[6];
    Num = new Fraction[2];

    for( int k = max(0,p - ((l+l_q)/2) - n_q);
        k <= min(n,(p - ((l+l_q)/2))); k ++){

        Num[0] = k + l + 1;
        Num[1] = p - k + ((l_q - l)/2) + 1;

        Den[0] = k;
        Den[1] = 2*(l + k + 1);
        Den[2] = n - k;
        Den[3] = 2*(p - k + 1)  + l_q - l;
        Den[4] = n_q - p + ((l + l_q)/2) + k;
        Den[5] = p - ((l + l_q)/2) - k;

        Sum += exp(factorial_function(Num,Den,2,6));

    }

    int * Num_C;
    int * Den_C;
    Num_C = new int[4];
    Num_C[0] = n;
    Num_C[1] = n_q;
    Num_C[2] = 2*(n + l + 1);
    Num_C[3] = 2*(n_q + l_q + 1);
    Den_C = new int [2];
    Den_C[0] = n + l + 1;
    Den_C[1] = n_q + l_q + 1;

    int *Num_CF;
    Num_CF = new int[1];
    Num_CF[0] = 2*(p + 1);
    int *Den_CF;
    Den_CF = new int[1];
    Den_CF[0] = p + 1;

    double Const = exp(factorial_function(Num_CF,Den_CF,1,1) - ((n+n_q)*log(2)) +
			  (0.5 * (factorial_function(Num_C,Den_C,4,2))));
    
    delete[] Num_C;
    delete[] Den_C;
    delete[] Num;
    delete[] Den;
    delete[] Num_CF;
    delete[] Den_CF;

    return pow(-1,int(p-((l+l_q)/2))) * Const * Sum / (pow(b_param,3));
}
//*/
/*
double B_coefficient(int n, int l, int n_q, int l_q, Fraction p, double b_param){

	// Extracted from paper: Symbolic algorithms for the computation of Moshinsky brackets
    //and nuclear matrix elements
    /// coefficient used in the sum decomposition  in Talmi integrals of the reduced matrix element of the potential V(r)
    if((n<0)||(n_q<0)||(l<0)||(l_q<0) ){return 0;}
    //cout << "B(n,l,n',l',p)" << n<< " "<< l<<" "<<n_q<<" "<<l_q<<" "<<p<<endl;
    double Sum = 0.;
    
    // Parity conservation
    if(((l+l_q)%2)!=0){return 0;}
    
    Fraction * Den;
    Fraction * Num;
    Den = new Fraction[6];
    Num = new Fraction[2];

    for( Fraction k = max(Fraction(0),(p - Fraction((l+l_q),2) - n_q));
        k<= min(Fraction(n),(p - Fraction((l+l_q),2))); k +=1){
        // k and p are fractions, the addition with them results in a fraction type
        Num[0] = k + l;
        Num[1] = p - Fraction((l-l_q),2) - k;

        Den[0] = k;
        Den[1] = 2*l + 2*k + 1;
        Den[2] = n - k;
        Den[3] = 2*p - l + l_q - 2*k + 1;
        Den[4] = n_q - p + Fraction((l + l_q),2) + k;
        Den[5] = p - Fraction((l + l_q),2) - k;

        Sum += exp(factorial_function(Num,Den,2,6));

    }

    int * Num_C;
    int * Den_C;
    Num_C = new int[4];
    Num_C[0] = n;
    Num_C[1] = n_q;
    Num_C[2] = 2*n + 2*l +1;
    Num_C[3] = 2*n_q + 2*l_q + 1;
    Den_C = new int [2];
    Den_C[0] = n + l;
    Den_C[1] = n_q + l_q;

    Fraction *Num_CF;
    Num_CF = new Fraction[1];
    Num_CF[0] = 2*p + 1;
    Fraction *Den_CF;
    Den_CF = new Fraction[1];
    Den_CF[0] = p;

    double Const = exp(factorial_function(Num_CF,Den_CF,1,1) -
		((n+n_q)*log(2)) + ((factorial_function(Num_C,Den_C,4,2))/2));
    delete[] Num_C;
    delete[] Den_C;
    delete[] Num;
    delete[] Den;
    delete[] Num_CF;
    delete[] Den_CF;

    return pow(-1,int(p-Fraction((l+l_q),2))) * Const * Sum/ (pow(b_param,3));

}
//*/

   // THIS VARIANT IS TAKING THE p INDEX AS INTEGER.
double B_coefficient(int n, int l, int n_q, int l_q, int p, double b_param){

	// Extracted from paper: Symbolic algorithms for the computation of Moshinsky brackets
    //and nuclear matrix elements
    /// coefficient used in the sum decomposition  in Talmi integrals of the reduced matrix element of the potential V(r)
    if((n<0)||(n_q<0)||(l<0)||(l_q<0) ){return 0;}
    //cout << "B(n,l,n',l',p)" << n<< " "<< l<<" "<<n_q<<" "<<l_q<<" "<<p<<endl;
    double Sum = 0.;
    
    // Parity conservation
    if(((l+l_q)%2)!=0){return 0;}
    
    int * Den;
    int * Num;
    Den = new int[6];
    Num = new int[2];

    for( int k = max(0,(p - ((l+l_q)/2) - n_q));
        k<= min(n,(p - ((l+l_q)/2))); k++ ){
        // k and p are fractions, the addition with them results in a fraction type
	
	Num[0] = k + l + 1;
        Num[1] = p - k + ((l_q - l)/2) + 1;

        Den[0] = k;
        Den[1] = 2*(l + k + 1);
        Den[2] = n - k;
        Den[3] = 2*(p - k + 1)  + l_q - l;
        Den[4] = n_q - p + ((l + l_q)/2) + k;
        Den[5] = p - ((l + l_q)/2) - k;

        Sum += exp(factorial_function(Num,Den,2,6));

    }

    int * Num_C;
    int * Den_C;
    Num_C = new int[4];
    Num_C[0] = n;
    Num_C[1] = n_q;
    Num_C[2] = 2*(n + l + 1);
    Num_C[3] = 2*(n_q + l_q + 1);
    Den_C = new int [2];
    Den_C[0] = n + l + 1;
    Den_C[1] = n_q + l_q + 1;

    int *Num_CF;
    Num_CF = new int[1];
    Num_CF[0] = 2*(p + 1);
    int *Den_CF;
    Den_CF = new int[1];
    Den_CF[0] = p + 1;

    double Const = exp(factorial_function(Num_CF,Den_CF,1,1) - ((n+n_q)*log(2)) +
			  (0.5 * (factorial_function(Num_C,Den_C,4,2))));
    
    delete[] Num_C;
    delete[] Den_C;
    delete[] Num;
    delete[] Den;
    delete[] Num_CF;
    delete[] Den_CF;

    return pow(-1,(p-((l+l_q)/2))) * Const * Sum/ (pow(b_param,3));

}
//*/

/////////////////    GAUSS-LEGENGRE NODES AND WEIGHTS FOR QUADRATURES    ////////////////////////////
/////////////////////////////////////////////////////////////////////////////////////////////////////

double **Gauss_Legendre_values = new double*[2];   // [0] -> Roots   	[1] -> Weights 

void run_Gauss_Legendre(int N){
    
    /*
    *  Functions to calculate integration points and weights for Gaussian quadrature
    *		by Mark Newman <mejn@umich.edu>, June 4, 2011 for Python 2/3
    */
    
    Gauss_Legendre_values[0] = new double[N]; // roots
    Gauss_Legendre_values[1] = new double[N]; // weights
    
    // Initial approximation to roots of the Legendre polynomial
    double * a;
    double * x;
    a = new double[N];
    x = new double[N];
    
    double da = 2.*(N-1) / ((N-1) * (2*N+1)); 
    for(int i=0; i<N; i++){
	a[i] = (3./(4*N+2)) + (da*i);
	x[i] = cos( (M_PI*a[i]) + (1./(8*pow(N,2)*tan(a[i]))) );
    }
    
    // Auxiliar  arrays for the roots in the iteration
    double * p_0;
    double * p_1;
    double * d_p;
    double * d_x;
    p_0 = new double[N];
    p_1 = new double[N];
    d_p = new double[N];
    d_x = new double[N];
    
    // Find roots using Newton's method
    double tolerance = 1e-15;
    double delta = 1.0;
    int iter = 0; 
    double aux_0, aux_1;
    while (delta>tolerance){
      
	// Actualization
	for(int i=0; i<N; i++){
	    p_0[i] = 1.; 
	    p_1[i] = x[i];
	}
        
        for( int k=1; k < N; k++){  // in range(1,N):
            //p0,p1 = p1,((2*k+1)*x*p1-k*p0)/(k+1)
	    for(int i=0; i<N; i++){
		aux_0 = p_0[i];
		aux_1 = p_1[i];
		
		p_0[i] = aux_1;
		p_1[i] = (((2*k+1)*x[i]*aux_1) - (k * aux_0)) / (k+1);
	    }
	    
	}
	
	for(int i=0; i<N; i++){
	    //dp = (N+1)*(p0-x*p1)/(1-x*x)
	    d_p[i] = (N+1)*(p_0[i] - (x[i]*p_1[i]))/ (1 - pow(x[i],2));
	    //dx = p1/dp 
	    // x -= dx
	    d_x[i] = p_1[i] / d_p[i];
	    x[i] = x[i] - d_x[i];
	    
	    delta = 0; // always lower than every d_x{i}
	    // take the maximum in absolute value
	    if(i != 0){
		if( fabs(d_x[i]) > delta ){
		    delta = fabs(d_x[i]);
		}
	    }
	    else{
		delta = d_x[i];
	    }
	    
	}
    }
    //for(int i=N-1; i>-1; i--){
    for(int i=0; i<N; i++){
	// Pass Roots and weights to the global variables
	Gauss_Legendre_values[0][i] = x[i]; 
	// Calculate the weights	w = 2*(N+1)*(N+1)/(N*N*(1-x*x)*dp*dp)
	Gauss_Legendre_values[1][i] = 2* pow( (N+1)/(N*d_p[i]) ,2) / (1-pow(x[i],2));
	
	//std::cout <<"Root= "<<Gauss_Legendre_values[0][i] << "  	Weight= " << Gauss_Legendre_values[1][i]<< std::endl;
    }
    
    // delete arrays
    delete [] a;
    delete [] x;
    delete [] p_0;
    delete [] p_1;
    delete [] d_p;
    delete [] d_x;
    //std::cout <<"Gauss-Legendre Values obtained (Order:"<< N <<")"<< std::endl;
}

double Gauss_Legendre_Weights (int ORDER, int k){
    // Returns the weights until the maximum order
    if(ORDER > 25){
	std::cout << "excedeed G-L quadrature maximum divisions ... (1 <= ORDER <= 8)" << std::endl;
	return 0;
    }
    if((k+1) > ORDER ){
	std::cout << "k(="<< k <<") excedeed the maximum interval ORDER(="<< ORDER <<")"<< std::endl;
	return 0;
    }
    switch(ORDER){
	case 1:
	    if     (k == 0){return 1.0000000000;}
	    break;
	case 2:
	    if     (k == 0){return 1.0000000000;}
	    else if(k == 1){return 1.0000000000;}
	    break; 
	case 3:
	    if     (k == 0){return 0.5555555556;}
	    else if(k == 1){return 0.8888888888;}
	    else if(k == 2){return 0.5555555556;}
	    break;
	case 4:
	    if     (k == 0){return 0.3478548451;}
	    else if(k == 1){return 0.6521451549;}
	    else if(k == 2){return 0.6521451549;}
	    else if(k == 3){return 0.3478548451;}
	    break;	    
	case 5:
	    if     (k == 0){return 0.2369268851;}
	    else if(k == 1){return 0.4786286705;}
	    else if(k == 2){return 0.5688888888;}
	    else if(k == 3){return 0.4786286705;}
	    else if(k == 4){return 0.2369268851;}
	    break;	
	case 6:
	    if     (k == 0){return 0.1713244924;}
	    else if(k == 1){return 0.3607615730;}
	    else if(k == 2){return 0.4679139346;}
	    else if(k == 3){return 0.4679139346;}
	    else if(k == 4){return 0.3607615730;}
	    else if(k == 5){return 0.1713244924;}
	    break;
	case 7:
	    if     (k == 0){return 0.1294849662;}
	    else if(k == 1){return 0.2797053915;}
	    else if(k == 2){return 0.3818300505;}
	    else if(k == 3){return 0.4179591837;}
	    else if(k == 4){return 0.3818300505;}
	    else if(k == 5){return 0.2797053915;}
	    else if(k == 6){return 0.1294849662;}
	    break;
	case 8:
	    if     (k == 0){return 0.1012285363;}
	    else if(k == 1){return 0.2223810345;}
	    else if(k == 2){return 0.3137066459;}
	    else if(k == 3){return 0.3626837834;}
	    else if(k == 4){return 0.3626837834;}
	    else if(k == 5){return 0.3137066459;}
	    else if(k == 6){return 0.2223810345;}
	    else if(k == 7){return 0.1012285363;}
	    break;	   
	case 25:
	    if     (k ==  0){return 0.0113937985010263;}
	    else if(k ==  1){return 0.0263549866150321;}
	    else if(k ==  2){return 0.0409391567013063;}
	    else if(k ==  3){return 0.0549046959758352;}
	    else if(k ==  4){return 0.0680383338123569;}
	    else if(k ==  5){return 0.0801407003350010;}
	    else if(k ==  6){return 0.0910282619829637;}
	    else if(k ==  7){return 0.1005359490670506;}
	    else if(k ==  8){return 0.1085196244742637;}
	    else if(k ==  9){return 0.1148582591457116;}
	    else if(k == 10){return 0.1194557635357848;}
	    else if(k == 11){return 0.1222424429903100;}
	    else if(k == 12){return 0.1231760537267154;}
	    else if(k == 13){return 0.1222424429903100;}
	    else if(k == 14){return 0.1194557635357848;}
	    else if(k == 15){return 0.1148582591457116;}
	    else if(k == 16){return 0.1085196244742637;}
	    else if(k == 17){return 0.1005359490670506;}
	    else if(k == 18){return 0.0910282619829637;}
	    else if(k == 19){return 0.0801407003350010;}
	    else if(k == 20){return 0.0680383338123569;}
	    else if(k == 21){return 0.0549046959758352;}
	    else if(k == 22){return 0.0409391567013063;}
	    else if(k == 23){return 0.0263549866150321;}
	    else if(k == 24){return 0.0113937985010263;}
	    break;
	default: 
	    std::cout << " ORDER not in options (Weights)"  << std::endl;
	    return 0;
    }
    
}


double Gauss_Legendre_Nodes (int ORDER, int k){
    // Returns the weights until the maximum order
    if(ORDER > 25){
	std::cout << "excedeed G-L quadrature maximum divisions ... (1 <= ORDER <= 8)" << std::endl;
	return 0;
    }
    if((k+1) > ORDER ){
	std::cout << "k(="<< k <<") excedeed the maximum interval ORDER(="<< ORDER <<")"<< std::endl;
	return 0;
    }
    switch(ORDER){
	case 1:
	    if     (k == 0){return  0.0000000000;}
	    break;
	case 2:
	    if     (k == 0){return -0.5773502692;}
	    else if(k == 1){return  0.5773502692;}
	    break; 
	case 3:
	    if     (k == 0){return -0.7745966692;}
	    else if(k == 1){return  0.0000000000;}
	    else if(k == 2){return  0.7745966692;}
	    break;
	case 4:
	    if     (k == 0){return -0.8611363116;}
	    else if(k == 1){return -0.3399810436;}
	    else if(k == 2){return  0.3399810436;}
	    else if(k == 3){return  0.8611363116;}
	    break;	    
	case 5:
	    if     (k == 0){return -0.9061798459;}
	    else if(k == 1){return -0.5384693101;}
	    else if(k == 2){return  0.0000000000;}
	    else if(k == 3){return  0.5384693101;}
	    else if(k == 4){return  0.9061798459;}
	    break;	
	case 6:
	    if     (k == 0){return -0.9324695142;}
	    else if(k == 1){return -0.6612063865;}
	    else if(k == 2){return -0.2386191861;}
	    else if(k == 3){return  0.2386191861;}
	    else if(k == 4){return  0.6612063865;}
	    else if(k == 5){return  0.9324695142;}
	    break;
	case 7:
	    if     (k == 0){return -0.9491079123;}
	    else if(k == 1){return -0.7415311856;}
	    else if(k == 2){return -0.4058451514;}
	    else if(k == 3){return  0.0000000000;}
	    else if(k == 4){return  0.4058451514;}
	    else if(k == 5){return  0.7415311856;}
	    else if(k == 6){return  0.9491079123;}
	    break;
	case 8:
	    if     (k == 0){return -0.9602898565;}
	    else if(k == 1){return -0.7966664774;}
	    else if(k == 2){return -0.5255324099;}
	    else if(k == 3){return -0.1834346425;}
	    else if(k == 4){return  0.1834346425;}
	    else if(k == 5){return  0.5255324099;}
	    else if(k == 6){return  0.7966664774;}
	    else if(k == 7){return  0.9602898565;}
	    break;
	 case 25:
	    if     (k ==  0){return -0.9955569697904981;}
	    else if(k ==  1){return -0.9766639214595175;}
	    else if(k ==  2){return -0.9429745712289743;}
	    else if(k ==  3){return -0.8949919978782753;}
	    else if(k ==  4){return -0.8334426287608340;}
	    else if(k ==  5){return -0.7592592630373576;}
	    else if(k ==  6){return -0.6735663684734684;}
	    else if(k ==  7){return -0.5776629302412229;}
	    else if(k ==  8){return -0.4730027314457150;}
	    else if(k ==  9){return -0.3611723058093879;}
	    else if(k ==  10){return -0.2438668837209884;}
	    else if(k ==  11){return -0.1228646926107104;}
	    else if(k ==  12){return 0.0000000000000000;}
	    else if(k ==  13){return 0.1228646926107104;}
	    else if(k ==  14){return 0.2438668837209884;}
	    else if(k ==  15){return 0.3611723058093879;}
	    else if(k ==  16){return 0.4730027314457150;}
	    else if(k ==  17){return 0.5776629302412229;}
	    else if(k ==  18){return 0.6735663684734684;}
	    else if(k ==  19){return 0.7592592630373576;}
	    else if(k ==  20){return 0.8334426287608340;}
	    else if(k ==  21){return 0.8949919978782753;}
	    else if(k ==  22){return 0.9429745712289743;}
	    else if(k ==  23){return 0.9766639214595175;}
	    else if(k ==  24){return 0.9955569697904981;}

	default: 
	    std::cout << " ORDER not in options (Nodes)" << std::endl;
	    return 0;
    }   
}

/////////////////////////////////////////////////////////////////////////////////////////////////////

// Minor auxiliary function for the calculation of normalization
int delta(int left, int right){
    if(left == right){return 1;}
    return 0;
}
// Minor auxiliary function to check if all jj Quantum Numbers are the same
int delta_2_1body_j(const QN_1body_jj & WF1,const QN_1body_jj & WF2){
    if(WF1.n == WF2.n){
        if(WF1.j == WF2.j){
            if(WF1.l == WF2.l){
                return 1;
            }
            else{return 0;}
        }
        else{return 0;}
    }
    else{return 0;}
}

// Calculation of the normalization for the Two Body interaction.
double Normalization_JT(const QN_2body_jj_Coupling &WF){
    
    if((WF.n1 == WF.n2) && (WF.l1 == WF.l2) && (WF.j1 == WF.j2)){
        return sqrt(1 - pow(-1,int(WF.J + WF.T)))/2;
    }
    else{
        return 1.;
    }
}


////////// SUHONEN INTERMEDIATE COEFFICIENTS (8.43-8.47) ////////////
///// IN MULTIPOLE DECOMPOSITION FOR M.E          ///////


double Y_Coeff(const QN_1body_jj & a, const QN_1body_jj & c, int lambda){
    /*std::cout << "Y_Coeff = " << Three_j_Symbol(a.j,c.j,Fraction(lambda), Fraction(1,2),Fraction(-1,2),Fraction(0)) *
                 pow(-1,double(c.j - Fraction(1,2) + lambda)) *
                 (1 + pow(-1,(a.l + c.l + lambda))) * sqrt( double((2*a.j + 1)*(2*c.j + 1))/M_PI)/4 << std::endl;*/
    
    return Three_j_Symbol(a.j,c.j,Fraction(lambda), Fraction(1,2),Fraction(-1,2),Fraction(0)) *
                 pow(-1,double(c.j - Fraction(1,2) + lambda)) *
                 (1 + pow(-1,(a.l + c.l + lambda))) * sqrt( double((2*a.j + 1)*(2*c.j + 1))/M_PI)/4 ;
}

double Omega_Coeff(const QN_2body_jj_Coupling &left,const QN_2body_jj_Coupling &right, int lambda, Fraction J){
    QN_1body_jj a,b,c,d;
    a = {left.n1,left.l1,left.j1};
    b = {left.n2,left.l2,left.j2};

    c = {right.n1,right.l1,right.j1};
    d = {right.n2,right.l2,right.j2};

    /*std::cout << "Omega_Coeff = " << pow(-1,double(b.j + c.j + J))*
            Six_j_Coefficient(a.j, b.j, J, d.j, c.j,Fraction(lambda))*
            Y_Coeff(a,c,lambda)*Y_Coeff(b,d,lambda) << std::endl;
	    */
	    
    return pow(-1,double(b.j + c.j + J))*
            Six_j_Coefficient(a.j, b.j, J, d.j, c.j,Fraction(lambda))*
            Y_Coeff(a,c,lambda)*Y_Coeff(b,d,lambda);
}

double Z_Coeff(const QN_1body_jj &a, const QN_1body_jj &c, int lambda, int j_index){
    /*
    std::cout << "Z_Coeff = " << sqrt(double((2*a.j+1)*(2*c.j+1)*(2*j_index+1)) *
                (2*a.l+1)*(2*c.l+1)*(2*lambda+1)/(4*M_PI))*pow(-1,a.l) *
                Three_j_Symbol(a.l,lambda,c.l, 0,0,0) *
                Nine_j_Symbol(Fraction(a.l),Fraction(1,2),a.j,
                              Fraction(c.l),Fraction(1,2),c.j,
                              Fraction(lambda),Fraction(1),Fraction(j_index)) << std::endl;
    */
			      
    return sqrt(double((2*a.j+1)*(2*c.j+1)*(2*j_index+1)) *
                (2*a.l+1)*(2*c.l+1)*(2*lambda+1)/(4*M_PI))*pow(-1,a.l) *
                Three_j_Symbol(a.l,lambda,c.l, 0,0,0) *
                Nine_j_Symbol(Fraction(a.l),Fraction(1,2),a.j,
                              Fraction(c.l),Fraction(1,2),c.j,
                              Fraction(lambda),Fraction(1), Fraction(j_index));
}

double Lambda_Coeff(const QN_2body_jj_Coupling &left, const QN_2body_jj_Coupling &right, int lambda, Fraction J){
    QN_1body_jj a,b,c,d;
    a = {left.n1,left.l1,left.j1};
    b = {left.n2,left.l2,left.j2};

    c = {right.n1,right.l1,right.j1};
    d = {right.n2,right.l2,right.j2};
    double Sum = 0.;
    
    //for(Fraction j_index = max(abs((a.j - c.j)),abs((d.j - b.j)));
    //        j_index <= min((a.j + c.j),(d.j + b.j)); j_index +=1)	// Asummed by me
    
    //int j_index_0 = abs(lambda - 1);
    //if(lambda == 0){ j_index_0 = 0;}
    //else{j_index_0 = lambda - 1;}
    
    // Talmi Limits (Simple Models for complex nuclei)
    for(int j_index = (abs(lambda - 1)); j_index <= (lambda + 1); j_index++){		
      
        Sum += pow(-1, double(j_index)) * Six_j_Coefficient(a.j, b.j, J, d.j, c.j,Fraction(j_index))*
		    Z_Coeff(a,c,lambda,int(j_index)) * Z_Coeff(b,d,lambda,j_index);
    }
    
    //std::cout << "Lambda_Coeff = " << pow(-1,double(b.j + c.j + J) + lambda + 1) * Sum << std::endl;
    return pow(-1,double(b.j + c.j + J) + lambda + 1) * Sum;
}

double U_Coeff(const QN_2body_jj_Coupling &left, const QN_2body_jj_Coupling &right, int lambda, double *Parameters ){
    // Constants for exchange operators and,
    double A = Parameters[0]; // V_0
    double B = Parameters[1]; // V_sigma
    double C = Parameters[2]; // V_tau
    double D = Parameters[3]; // V_tau_sigma
    
    //std::cout <<"A,B,C,D  = "<<A<<" "<<B<<" "<< C<<" "<<D<<std::endl;
    
    int T = left.T;
    Fraction J = left.J;
    
    //std::cout << "U_coeff =" << (A + C*(delta(T,1) - 3*delta(T,0)))* Omega_Coeff(left,right,lambda,J) +
    //        ((B + D*(delta(T,1) - 3*delta(T,0)))* Lambda_Coeff(left,right,lambda,J)) << std::endl;
    
    double U_AC = 0;
    double U_BD = 0;
    if( (fabs(A) > 1e-10) || (fabs(C) > 1e-10) ){
	U_AC = (A + C*(delta(T,1) - 3*delta(T,0)))* Omega_Coeff (left,right,lambda,J);
    }
    if( (fabs(B) > 1e-10) || (fabs(D) > 1e-10) ){
	U_BD = (B + D*(delta(T,1) - 3*delta(T,0)))* Lambda_Coeff(left,right,lambda,J);
    }
    
    return U_AC + U_BD;
}

////////// SUHONEN INTERMEDIATE COEFFICIENTS ////////////
/////        FOR SURFACE DELTA INTERACTION        ///////

double K_Coeff(const QN_2body_jj_Coupling & left, const QN_2body_jj_Coupling &right){
    // This is useful if we use the A_T as a fitting parameter.
    //
    double A_0 = 1.0;  // MeV
    double A_1 = 1.0;  // MeV
    
    int option = 0;
    switch(option){
        case 0:
            {
                // Constant integral value with fitting parameter A_T
                double A_T;
                A_T = (left.T==1)? A_1:A_0;
                return pow(-1, left.n1 + left.n2 + right.n1 + right.n2 + 1)* A_T / 4;
            }

        case 1:
            {
                // Explicit calculation of the integral
                QN_2body_radial ac = {left.n1, left.l1, right.n1, right.l1 };
                QN_2body_radial bd = {left.n2, left.l2, right.n2, right.l2 };

                return - V_0 * Radial_SDI(ac)*Radial_SDI(bd)/(16*M_PI);
            }
    }
}



//////////    TALMAN ARTICLE COEFFICIENTS,     ////////////
/////    EVALUATION OF RADIAL INTEGRALS FOR M.E     ///////

double c_normalization_Tal( int n, int l){
    // Normalization
    return pow(-1,n) * sqrt(exp(factorial(n) + double_factorial(2*(l + n) + 1)) /
                            (pow(M_PI,double(3)/2)*(2*l+1)*pow(2,l+n)));
}

double Big_C_coeff_Tal(const QN_2body_radial & Arguments, int S, int L){
    // This value arise from the coupling of two SHO functions in terms of one SHO
    // l + l_q have the same parity than L
    int s = Arguments.n1;
    int l = Arguments.l1;
    int s_q = Arguments.n2;
    int l_q = Arguments.l2;

    double Aux = pow(-1,s + s_q)* (2*l+1)* (2*l_q+1)* pow(2,double(l+l_q - L)/2) *
                Three_j_Symbol(l,l_q,L, 0,0,0);

    //std::cout << "Big C::=============" << std::endl;
    //std::cout << "s,l,s',l/ S,L' = "<<s<<"," <<l<<"," <<s_q<<","<<l_q<<"/ "<<S<<","<<L<< std::endl;
    //std::cout << "Aux =" << Aux << std::endl;
    double Sum = 0.;

    // The sum goes for all mu, nu index for which the factorials-double factorials are positive
    // the upper limit is when the sum get one of the minimum values of this
    int l_aux = (l + l_q -L )/2;
    //int upper_limit = min(l_aux,l_aux + S); // = mu + nu
    //std::cout << "upper =" << max(max(S - l_aux,0) - s,0) << std::endl;

    int *Num;
    int *Den;
    Num = new int[1];
    Den = new int[5];

    int *D_Num;
    int *D_Den;
    D_Num = new int[1];
    D_Den = new int[2];

    // (the index goes opposite)
    for (int mu = 0; mu <= s; mu++){
        for(int nu = s_q; nu >= max(S - l_aux - mu,0); nu--){

            Num[0] = l_aux + mu + nu;
            Den[0] = l_aux + mu + nu  - S;
            Den[1] = mu;
            Den[2] = nu;
            Den[3] = s - mu;
            Den[4] = s_q - nu;

            D_Num[0] = l + l_q + L + 2*(mu + nu) +1;
            D_Den[0] = 2*(mu + l) +1;
            D_Den[1] = 2*(nu + l_q) +1;

            Sum += pow(-1,mu + nu) *
                    exp(double_factorial_function(D_Num,D_Den,1,2) + factorial_function(Num,Den,1,5));


        }
    }
    //std::cout << "Sum =" << Sum<< std::endl;
    delete []Num;
    delete []Den;
    delete []D_Num;
    delete []D_Den;

    return Aux * Sum;

}

double R_coeff_Tal( int S_1, int S_2, int L, int nu){
    // This coefficient carry the radial integral information
    // Talmi integrals defined in Talman article differ in the normalization constant,
    // that is why the coefficient differs from the paper formula (48)

    // factorial function
    int *Num;
    int *Den;
    Num = new int[1];
    Den = new int[4];

    Num[0] = L + S_1 + S_2;
    Den[0] = S_1;
    Den[1] = S_2;
    Den[2] = nu;
    Den[3] = L + S_1 + S_2 - nu;
    
    // double factorial function
    int *D_Num;
    int *D_Den;
    D_Num = new int[1];
    D_Den = new int[2];

    D_Num[0] = 2*(L + S_1 + S_2) + 1;
    D_Den[0] = 2*(L + S_1) + 1;
    D_Den[1] = 2*(L + S_2) + 1;

    double Aux_factorial_function = factorial_function(Num,Den,1,4) +
                                    double_factorial_function(D_Num,D_Den,1,2);

    delete [] Num;
    delete [] Den;
    delete [] D_Num;
    delete [] D_Den;

    return pow(M_PI,3) * pow(-1, S_1 + S_2 + nu) * (2*L + 1)*
            exp(Aux_factorial_function);
}
