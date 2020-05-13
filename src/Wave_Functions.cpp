#include "../include/Fractions.h"
#include "../include/factorials.h"
#include "../include/Angular_Functions.h"
#include "../include/Index_Coefficients.h"
#define _USE_MATH_DEFINES
#include <cmath>

double SHO_normalization(int n, int l, double b_param){
    // definition Suhonen of the normalization factorial
    //return sqrt(1/(sqrt(M_PI)*pow(b_lenght(100),3))) *
    //    exp(0.5*(factorial(n) + (n+l+2)*log(2) - double_factorial(2*(n + l) - 1)));

    // explicit definition
    //int aux_NULL = NULL;  
    //return sqrt((2 * exp(factorial(n)) * pow(1/b_param,3 + 2*l) ) /
    //    pow(exp(gamma_function(Fraction(2*n + 2*l + 3, 2),aux_NULL)),1) );
    int *Num;
    int *Den;
    Num = new int[2];
    Den = new int[1];
    
    Num[0] = n;
    Num[1] = n + l + 1;
    Den[0] = 2*(n + l + 1);
    
    double aux = factorial_function(Num,Den,2,1);
    delete [] Num;
    delete [] Den;
    
    return sqrt( (exp(aux + ((2*(n+l+1) + 1)*log(2))) ) /
        ( sqrt(M_PI)*pow(b_param,3 + 2*l)) );
}

double SHO_Radial(int n, int l, double r, double b_param){
    
    // This function is design for a Cuadrature calculation for 2body integral
    
    double Aux = 0.;

    // Laguerre Polynomial in terms of the series.
    int *Num;
    int *Den;
    Num = new int[2];
    Den = new int[4];

    for(int k = 0; k<=n; k++){
	Num[0] = 2*(l + n + 1);
	Num[1] = (l + k + 1) ;
	
        Den[0] = 2*(l + k + 1);
        Den[1] = (l + n + 1);
        Den[2] = n - k;
	Den[3] = k;

        Aux += pow(-1,k) * exp( ((2*(k - n))* log(2)) + factorial_function(Num,Den,2,4)) * pow(r / b_param, 2*k);
    }
    delete [] Num;
    delete [] Den;

    // constant part
    Aux *= SHO_normalization(n,l,b_param);

    // other radial part
    return Aux * pow(r,l) * exp( (-0.5) * pow((r/b_param),2));
}

double Legendre_Polynomial (int n, double x){
    
    if (n < 0){
	std::cout <<"(ERROR) Leg Poly (n<0) "<<std::endl; 
	return 0;
    }
    if (fabs(x) > 1.0 ){
	std::cout <<"(ERROR) Leg Poly (|x|>1) "<<std::endl; 
	return 0;      
    }
    
    double Aux = 0.;
    
    int *Num;
    int *Den;
    Num = new int[1];
    Den = new int[3];

    for(int k = 0; k <= floor(double(n)/2); k++){
	Num[0] = 2*(n - k);
	
        Den[0] = (n - (2*k));
        Den[1] = n - k;
	Den[2] = k;

        Aux += pow(-1,k) * exp(factorial_function(Num,Den,1,3)) * pow(x, n - 2*k);
    }
    delete [] Num;
    delete [] Den;
    
    return Aux / pow(2 , n) ;
    
    
}

