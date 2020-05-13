#include <iostream>
#define _USE_MATH_DEFINES
#include <cmath>
#include "../include/factorials.h"
#include "../include/Fractions.h"
#include "../include/Index_Coefficients.h"
#include "../include/Wave_Functions.h"
#include "../include/Angular_Functions.h"

#include "../include/Nuclear_Matrix_Elements_BB.h"

//double oscillator_lenght = 1.;
double Talmi_Integral(int p , int option, double b_param){
    // Talmi tesis, mu = omega*m/h_bar = hbar_omega* mc^2/ (hbar_c)^2
    // (mu_Heyde = mu_Talmi/2)
    bool Set_Manual_Lambda = true;
    
    double lambda;
    if(Set_Manual_Lambda){
	//std::cout << "Manual LAMBDA="<< LAMBDA << std::endl; 
        lambda = LAMBDA;
    }
    else{
	double r_0 = 1.2;  // fm
        lambda = r_0 / b_lenght(A) ;
    }


    // Choose of type potential(each one of them must be multiplied by V constant):
    //int option = 0;
    switch( option ){
        case(0):{
            // V exp(-(r/b)^2)
	    //  Required:  LAMBDA = b_param/ mu_param;
            //return pow(pow(lambda,2)/(1+pow(lambda,2)),double(p + Fraction(3/2))); // this have been derived in units of b_param
	    return 1/pow(1 + pow(lambda,2),double(p + (3./2))); //  pow(b_param,3) 
        }
        // all options are in units of b_param
        case(1):{
            // V exp(-(r/b)^2) / (r/b)
            // normalization used in the Tesis is diferent,
            // so I put directly as a constant
            return ((pow(2, int(p)+1) * factorial(p)) /
                    (sqrt(M_PI)* double_factorial(int(2*p + 1))))*
                    lambda*(pow(lambda,2)/(1+pow(lambda,2)),double(p + 1));
        }
        case(2):{
            // V exp(-(r/b)^2) / (r/b)^2
            return (pow(lambda, 2)* 2 / double((2*p + 1)))*
                    (pow(lambda,2)/(1+pow(lambda,2)),double(p + (1/2)));
        }
	case(3):{
	    // Yukawa potential (not in b units) V = exp(-r/mu)/(r/mu)
	    //  Required:  LAMBDA = b_param/ mu_param;
	    
	    double Sum = 0.;
	    
	    int * Num;
	    Num = new int[1];
	    int * Den;
	    Den = new int[2];
	    int sign = NULL; // not necessary , i >= 0
	    
	    for( int i = 0; i <= (2*p + 1); i++){
		Num[0] = 2*p + 1;
		
		Den[0] = i;
		Den[1] = 2*p + 1 - i;
		
		Sum += exp( gamma_function(Fraction(i+1,2) ,sign) * factorial_function(Num,Den,1,2))* 
			pow( ((-1) * 0.5 * LAMBDA), 2*p + 1 - i );
		
	    }
	    
	    delete [] Num;
	    delete [] Den;
	    
	    return ((pow(b_param,3)* exp( pow((LAMBDA /2),2) + double_factorial(2*p + 1))) / LAMBDA ) 
		    * Sum; 
	    
	}
	
	case(4):
	    // r^2 potential
	    return pow(b_param,5) * (p + 3./2);
    }
}

/*
double Talmi_Integral(Fraction p){
    // Talmi tesis, mu = omega*m/h_bar = hbar_omega* mc^2/ (hbar_c)^2
    // (mu_Heyde = mu_Talmi/2)

    double r_0 = 1.2;  // fm

    double lambda;
    if(Set_Manual_Lambda){
        lambda = LAMBDA;
    }
    else{
        lambda = r_0 / b_lenght(A) ;
    }


    // Choose of type potential(each one of them must be multiplied by V constant):
    int option = 0;
    switch( option ){
        case(0):{
            // V exp(-(r/b)^2)
            //return pow(pow(lambda,2)/(1+pow(lambda,2)),double(p + Fraction(3/2)));
	    return ;	    
        }
        case(1):{
            // V exp(-(r/b)^2) / (r/b)
            // normalization used in the Tesis is diferent,
            // so I put directly as a constant
            return ((pow(2, int(p)+1) * factorial(p)) /
                    (sqrt(M_PI)* double_factorial(int(2*p + 1))))*
                    lambda*(pow(lambda,2)/(1+pow(lambda,2)),double(p + 1));
        }
        case(2):{
            // V exp(-(r/b)^2) / (r/b)^2
            return (pow(lambda, 2)* 2 / double((2*p + 1)))*
                    (pow(lambda,2)/(1+pow(lambda,2)),double(p + Fraction(1/2)));
        }
    }
}
*/

///*
double Multipolar_Rosenfeld(int lambda, int GL_ORDER, double r_1, double r_2, double mu_param){
    //std::cout << "Multipolar_Rosenfeld( lam " << lambda << ", r1 "<< r_1 << ",r2 "<< r_2 <<")"<<std::endl;
    // Evaluation of the kernel v_lambda(r_1,r_2) of the multipolar descomposition
    // of the Rosenfeld central potential (which is a gaussian), 
    //
    // In this case by gaussian quadratures 
    
    double A = pow(r_1,2) + pow(r_2,2);
    double B = 2 * r_1 * r_2;
    //std::cout <<"A,B ="<< A <<" ,"<<B<< std::endl;
    
    double Sum;
    double Sum_Gauss = 0.;
        
    int SUBDIVISIONS = 20; 
    double dx = 2. / SUBDIVISIONS;
    double x, x_min,x_max;
    /*
    // GAUSSIAN QUADRATURE OF THE SUBDIVIDED INTERVALS (more accurate and slower)
    for(int SUB=0; SUB < (SUBDIVISIONS); SUB++){
	x_min = -1. + (SUB*dx);
	x_max = -1. + ((SUB+1)*dx);
	
	Sum = 0.;
	for(int k = 0; k < GL_ORDER; k++){	    
	    x = (0.5*(x_max + x_min )) + (0.5*(x_max - x_min)* Gauss_Legendre_values[0][k]);
	    Sum +=  Gauss_Legendre_values[1][k] * exp(-(A-B*x)/pow(mu_param,2)) * Legendre_Polynomial(lambda,x) ;
	    //std::cout <<"x: "<<x<<"   GLW= "<<Gauss_Legendre_values[1][k]<< " , exp= "<<exp(-(A-B*x)/pow(mu_param,2))<< " ,PLeg "<< Legendre_Polynomial(lambda,x)<<std::endl;
	    //std::cout <<"	Sum= "<<Sum<< std::endl;
	}
	Sum_Gauss += 0.5* dx * Sum;
    }// */
    
    ///*
    // GAUSSIAN QUADRATURE OF THE WHOLE RANGE -1,1
    for(int k = 0; k < GL_ORDER; k++){
	x = Gauss_Legendre_values[0][k];
	Sum_Gauss +=  Gauss_Legendre_values[1][k] * exp(-(A-B*x)/pow(mu_param,2)) * Legendre_Polynomial(lambda,x) ;
	//std::cout <<"	GLW= "<<Gauss_Legendre_values[1][k]<< " , exp= "<<exp(-(A-B*x)/pow(mu_param,2))<< " ,PLeg "<< Legendre_Polynomial(lambda,x)<<std::endl;
	//std::cout <<"	Sum= "<<Sum<< std::endl;
	
    }// * /
    //std::cout << "Rosenfeld Multipolar=" << (2*M_PI)* Sum_Gauss  << std::endl; 
    return (2*M_PI)* Sum_Gauss;
}
//*/

/*
double Multipolar_Rosenfeld(int lambda,int GL_ORDER, double r_1, double r_2, double mu_param){
    //std::cout << "Multipolar_Rosenfeld( lam " << lambda << ", r1 "<< r_1 << ",r2 "<< r_2 <<")"<<std::endl;
  
    //				 ANALITIC INTEGRAL				//
    // Evaluation of the kernel v_lambda(r_1,r_2) of the multipolar 		//
    // descomposition of the Rosenfeld central potential (which is a gaussian), //
    // by direct integration of explicit Legendre polinomials.			//
    
    double A = pow(r_1,2) + pow(r_2,2);
    double B = 2 * r_1 * r_2;
   
    double Sum_1 = 0.;
    double Sum_2 =  0.;
    double Aux;
    double aux_B = B/ pow(mu_param,2);
    double exp_BpA = exp((-1)*(A+B) / pow(mu_param,2));
    double exp_BmA = exp((B-A) / pow(mu_param,2));
    
    //std::cout <<"A, B="<< A << " , " << B << " , mu= "<< mu_param<< "  / BpA , BmA ="<<exp_BpA << " , "<<exp_BmA <<std::endl;
     
    int *Num;
    int *Den;
    Num = new int[1];
    Den = new int[2];
    
    for(int k = 0; k <= floor(lambda/2); k++){
	Num[0] = 2*(lambda - k);
	Den[0] = k;
	Den[1] = lambda - k;
	
	Aux = pow(-1,k) * exp(factorial_function(Num,Den,1,2));
	//std::cout << " Aux=" << Aux << std::endl;
	
	Sum_2 =  0.;
	for(int i = 0; i <= (lambda - 2*k); i++){
	    
	    //if(exp(aux_BmA) < 1.e-15){exp_BmA = 0;}
	    //if(exp((-1)*aux_BpA) < 1.e-15){exp_BpA = 0;}
	    //std::cout << "Sum + Num(" << (exp_BmA - pow(-1, lambda - (2*k) - i)*exp_BpA)<<") / Den("<< 
					  //(pow(aux_B,i) * exp(factorial(lambda - (2*k) - i)))<<")" <<std::endl;
	    //std::cout << "powB = " <<pow(aux_B,i)<< std::endl;
	    //Sum_2 += (pow(-1,i)*(exp(aux_BmA) - pow(-1, lambda - (2*k) - i)*exp((-1)*aux_BpA)))/ (pow(aux_B,i+1) * exp(factorial(lambda - (2*k) - i))); 
	    Sum_2 += (pow(-1,i)*(exp_BmA - pow(-1, lambda - (2*k) - i)*exp_BpA))/ (pow(aux_B,i+1) * exp(factorial(lambda - (2*k) - i))); 
	}
	
	Sum_1 += Aux * Sum_2;
    }
    
    delete [] Num;
    delete [] Den;
    
    //std::cout << "Rosenfeld Multipolar=" << (2*M_PI)* Sum_1 / (pow(2,lambda)) << std::endl;
    bool continua;
    //std::cin >> continua;
    return (2*M_PI)* Sum_1 / (pow(2,lambda));
}
//*/

double Radial_4_order_Subdivided_Quadrature(const QN_2body_radial &Q_Numbers_left,
				       const QN_2body_radial &Q_Numbers_right, int ORDER, 
				       int lambda, double mu_param, double b_param){
    // GAUSS QUADRATURE INTEGRATION of separated subdivisions of the interval,
    // 		the upgrade requires a lower number of nodes, to achieve the accuracy
  
    // GAUSS QUADRATURE INTEGRATION PARAMETERS 
    double R_MIN =  0.; //
    double R_MAX = 25.; // fm
    int Subdivisions = 25; // Number of splitings of the range R_MIN - R_MAX
    
    // passing of 2body Wave_Functions to 1body radial Wave_Functions
    int n_a = Q_Numbers_left.n1;
    int l_a = Q_Numbers_left.l1;
    int n_b = Q_Numbers_left.n2;
    int l_b = Q_Numbers_left.l2;

    int n_c = Q_Numbers_right.n1;
    int l_c = Q_Numbers_right.l1;
    int n_d = Q_Numbers_right.n2;
    int l_d = Q_Numbers_right.l2;
    
    double Aux_1 = 0.;
    double Gauss_Sum = 0.;
    double Subdivision_Sum_1,Subdivision_Sum_2;
    double r_1_k1, r_2_k2;
    double r_1_min,r_1_max,r_2_min,r_2_max;
    
    double d_interval = (R_MAX - R_MIN) / Subdivisions;
    
    // First set of subdivisions
    for(int SUB_1 = 0; SUB_1 < (Subdivisions); SUB_1++){
	
	r_1_min = R_MIN + (SUB_1 * d_interval);
	r_1_max = R_MIN + ((SUB_1 + 1) * d_interval);
	Subdivision_Sum_1 = 0.;
	for(int SUB_2 = 0; SUB_2 < (Subdivisions); SUB_2++){
	      r_2_min = R_MIN + (SUB_2 * d_interval);
	      r_2_max = R_MIN + ((SUB_2 + 1) * d_interval);
	      
	      Subdivision_Sum_2 = 0.;
	      //std::cout << " / r_1_min: " <<r_1_min << " / r_1_max: " <<r_1_max 
		//	<< " / r_2_min: " <<r_2_min << " / r_2_max: " <<r_2_max<<std::endl;
	      for( int k1 = 0; k1 < ORDER; k1++){
		  r_1_k1 = (0.5*(r_1_max + r_1_min )) + (0.5*(r_1_max - r_1_min)* Gauss_Legendre_values[0][k1]);
		  
		  Aux_1 = Gauss_Legendre_values[1][k1] * pow(r_1_k1,2) * 
				SHO_Radial(n_a,l_a, r_1_k1 ,b_param) * SHO_Radial(n_c,l_c, r_1_k1 ,b_param);
				
		  // Second set of subdivisions	  		
		  for( int k2 = 0; k2 < ORDER; k2++){
		      
		      r_2_k2 = (0.5*(r_2_max + r_2_min )) + (0.5*(r_2_max - r_2_min)* Gauss_Legendre_values[0][k2]);
		      
		      Subdivision_Sum_2 += Aux_1 * Gauss_Legendre_values[1][k2] * pow(r_2_k2,2) *
				    Multipolar_Rosenfeld(lambda, ORDER, r_1_k1,  r_2_k2, mu_param) *  
					SHO_Radial(n_b,l_b, r_2_k2 ,b_param) * SHO_Radial(n_d,l_d, r_2_k2 ,b_param);
		  }
	      }
	      
	      Subdivision_Sum_1 += pow(0.5*(r_2_max - r_2_min), 2) * Subdivision_Sum_2;
	}
	Gauss_Sum += Subdivision_Sum_1;
	
    }
    
    return Gauss_Sum;
}
//*/

double Radial_4_order_Gauss_Quadrature(const QN_2body_radial &Q_Numbers_left,
				       const QN_2body_radial &Q_Numbers_right, int GL_ORDER, 
				       int lambda, double mu_param, double b_param){
    // GAUSS QUADRATURE INTEGRATION PARAMETERS 
    double R_MIN =  0.; // fm
    double R_MAX = 25.; // fm
    
    // passing of 2body Wave_Functions to 1body radial Wave_Functions
    int n_a = Q_Numbers_left.n1;
    int l_a = Q_Numbers_left.l1;
    int n_b = Q_Numbers_left.n2;
    int l_b = Q_Numbers_left.l2;

    int n_c = Q_Numbers_right.n1;
    int l_c = Q_Numbers_right.l1;
    int n_d = Q_Numbers_right.n2;
    int l_d = Q_Numbers_right.l2;
    
    double Aux_1 = 0.;
    double Gauss_Sum = 0.;
    double r_1_k1, r_2_k2;
    
    for( int k1 = 0; k1 < GL_ORDER; k1++){
	r_1_k1 = (0.5*(R_MAX + R_MIN )) + (0.5*(R_MAX - R_MIN)* Gauss_Legendre_values[0][k1]);
	
	Aux_1 = Gauss_Legendre_values[1][k1] * pow(r_1_k1,2) * 
		      SHO_Radial(n_a,l_a, r_1_k1 ,b_param) * SHO_Radial(n_c,l_c, r_1_k1 ,b_param);
	
	//std::cout << "*[r1,SHOa,SHOb,r1^2,GLW,Aux] = "<< r_1_k1<<" , "<<SHO_Radial(n_a,l_a, r_1_k1 ,b_param)<<" , "<< SHO_Radial(n_c,l_c, r_1_k1 ,b_param)
	//	  <<" , "<< pow(r_1_k1,2)<<" , "<< Gauss_Legendre_values[1][k1] <<" , "<< Aux_1 <<std::endl;
	for( int k2 = 0; k2 < GL_ORDER; k2++){
	    r_2_k2 = (0.5*(R_MAX + R_MIN )) + (0.5*(R_MAX - R_MIN)* Gauss_Legendre_values[0][k2]);
	    
	//    std::cout << " [r2,SHOb,SHOd,MRos,r2^2,GLW] = "<< r_2_k2 <<" , "<<SHO_Radial(n_b,l_b, r_2_k2 ,b_param)<<" , "<<SHO_Radial(n_d,l_d, r_2_k2 ,b_param)
	//	      <<" , "<< Multipolar_Rosenfeld(lambda, GL_ORDER,  r_1_k1,  r_2_k2, mu_param) <<" , "<< pow(r_2_k2,2)<<" , "<< Gauss_Legendre_values[1][k1] <<std::endl;
	    Gauss_Sum += Aux_1 * Gauss_Legendre_values[1][k2] * pow(r_2_k2,2) *
			   Multipolar_Rosenfeld(lambda, GL_ORDER , r_1_k1,  r_2_k2, mu_param) *  
			      SHO_Radial(n_b,l_b, r_2_k2 ,b_param) * SHO_Radial(n_d,l_d, r_2_k2 ,b_param);
	}
	//std::cout <<"	Sum= "<<Gauss_Sum<< std::endl;
    }
    //std::cout << " Gauss_sum: " <<Gauss_Sum <<std::endl;
    return pow(0.5*(R_MAX - R_MIN), 2) * Gauss_Sum;
    
}

double Rosenfeld_quadratures(int n_a,int l_a,int n_b,int l_b,int L,int n_c,int l_c,int n_d,int l_d,int Lq,int GL_ORDER,double mu_param ,double b_param){
    // GAUSS QUADRATURE INTEGRATION PARAMETERS 
    double R_MIN =  0.; // fm
    double R_MAX = 25.; // fm
    
    double Aux_1 = 0.;
    double Gauss_Sum = 0.;
    double r_1_k1, r_2_k2;
    
    for( int k1 = 0; k1 < GL_ORDER; k1++){
	r_1_k1 = (0.5*(R_MAX + R_MIN )) + (0.5*(R_MAX - R_MIN)* Gauss_Legendre_values[0][k1]);
	
	Aux_1 = Gauss_Legendre_values[1][k1] * pow(r_1_k1,2) * 
		      SHO_Radial(n_a,l_a, r_1_k1 ,b_param) * SHO_Radial(n_c,l_c, r_1_k1 ,b_param);
	
	//std::cout << "*[r1,SHOa,SHOb,r1^2,GLW,Aux] = "<< r_1_k1<<" , "<<SHO_Radial(n_a,l_a, r_1_k1 ,b_param)<<" , "<< SHO_Radial(n_c,l_c, r_1_k1 ,b_param)
	//	  <<" , "<< pow(r_1_k1,2)<<" , "<< Gauss_Legendre_values[1][k1] <<" , "<< Aux_1 <<std::endl;
	for( int k2 = 0; k2 < GL_ORDER; k2++){
	    r_2_k2 = (0.5*(R_MAX + R_MIN )) + (0.5*(R_MAX - R_MIN)* Gauss_Legendre_values[0][k2]);
	    
	//    std::cout << " [r2,SHOb,SHOd,MRos,r2^2,GLW] = "<< r_2_k2 <<" , "<<SHO_Radial(n_b,l_b, r_2_k2 ,b_param)<<" , "<<SHO_Radial(n_d,l_d, r_2_k2 ,b_param)
	//	      <<" , "<< Multipolar_Rosenfeld(lambda, GL_ORDER,  r_1_k1,  r_2_k2, mu_param) <<" , "<< pow(r_2_k2,2)<<" , "<< Gauss_Legendre_values[1][k1] <<std::endl;
	    Gauss_Sum += Aux_1 * Gauss_Legendre_values[1][k2] * pow(r_2_k2,2) *
			   exp(- pow(( r_1_k1 - r_2_k2),2) / pow(mu_param,2)) *  
			      SHO_Radial(n_b,l_b, r_2_k2 ,b_param) * SHO_Radial(n_d,l_d, r_2_k2 ,b_param);
	}
	//std::cout <<"	Sum= "<<Gauss_Sum<< std::endl;
    }
    //std::cout << " Gauss_sum: " <<Gauss_Sum <<std::endl;
    return pow(0.5*(R_MAX - R_MIN), 2) * Gauss_Sum;
}


/*//   EVALUATION METHOD OF THE TWO BODY MATRIX ELEMENTS BY GENERATING FUNCTIONS 
double Radial_4_order(const QN_2body_radial &Q_N_radial_left,const QN_2body_radial &Q_N_radial_right, int lambda){

    // Index passing
    // n >= 0
    int n_a = Q_N_radial_left.n1;
    int l_a = Q_N_radial_left.l1;
    int n_b = Q_N_radial_left.n2;
    int l_b = Q_N_radial_left.l2;

    int n_c = Q_N_radial_right.n1;
    int l_c = Q_N_radial_right.l1;
    int n_d = Q_N_radial_right.n2;
    int l_d = Q_N_radial_right.l2;

    double Sum = 0.;
    double aux_6j,aux_C1,aux_C2;
    QN_2body_radial Q2B_1, Q2B_2;

    bool continua;

    for(int L = max(abs(l_a - l_c),abs(l_b - l_d)); L <= min(l_a + l_c,l_b + l_d); L += 2){
        //std::cout<<"L,lambda="<<L<<","<<lambda<<std::endl;
        aux_6j = Six_j_Coeficient(l_c,l_d,lambda, l_b,l_a,L);

        if(fabs(aux_6j) > 1e-10){
            // S1, and S2 limits come from Big_C (S <= (l+l'-L)/2 + minimum(mu + nu))
            for(int S_1 = max(((l_a + l_c - L)/2),0); S_1 <= n_a + n_c; S_1++){

                Q2B_1 = {n_a,l_a,n_c,l_c};

                aux_C1 = Big_C_coeff_Tal(Q2B_1,S_1,L);


                if(fabs(aux_C1) > 1e-10){
                    for(int S_2 = max(((l_b + l_d - L)/2),0); S_2 <= n_b + n_d; S_2++){

                        Q2B_2 = {n_b,l_b,n_d,l_d};

                        aux_C2 = Big_C_coeff_Tal(Q2B_2,S_2,L);

                        if(fabs(aux_C1) > 1e-10){
                            // nu limit is derived from R(S1,S2,L,nu) integral (48)
                            for(int nu = 0; nu <= (L + S_1 + S_2); nu++){

                                //std::cout << "   C1 =" << aux_C1 << std::endl;
                                //std::cout << "   C2 =" << aux_C2 << std::endl;
                                //std::cout << "  6-j =" << aux_6j << std::endl;
                                //std::cout << "    R =" << R_coeff_Tal(S_1,S_2,L,nu) << std::endl;
                                //std::cout << "Talmi =" << Talmi_Integral(Fraction(nu)) << std::endl;

                                Sum += aux_C1 * aux_C2 * aux_6j *
                                         R_coeff_Tal(S_1,S_2,L,nu) * Talmi_Integral(Fraction(nu), b_lenght(A));
                                //std::cout << "       * Sum    =" << Sum << std::endl;


                            }
                        }
                        else{}

                    }
                }
                else{}

            }
        }
        else{}

    }


    //std::cout << "Radial_4_order=" << Sum << std::endl;
    //std::cin >> continua;
    //std::cout << "-------------------" <<std::endl;
    return Sum * c_normalization_Tal(n_a,l_a)*c_normalization_Tal(n_b,l_b)*
                 c_normalization_Tal(n_c,l_c)*c_normalization_Tal(n_d,l_d)*
                 pow(-1,l_b + l_d - lambda);

}
//*/

/*//
// Developed by me, does not run.
double Radial_4_order( const QN_2body_radial &Q_N_radial_left,const QN_2body_radial &Q_N_radial_right, int lambda){

    // Index passing
    // n >= 0
    int n_a = Q_N_radial_left.n1;
    int l_a = Q_N_radial_left.l1;
    int n_b = Q_N_radial_left.n2;
    int l_b = Q_N_radial_left.l2;

    int n_c = Q_N_radial_right.n1;
    int l_c = Q_N_radial_right.l1;
    int n_d = Q_N_radial_right.n2;
    int l_d = Q_N_radial_right.l2;

    if(((3 + l_b + l_d) < lambda) || ((3 + l_a + l_c) < lambda)){
        std::cout << "Lambda exceed the integral limit, decrease the order" <<std::endl;
        return 0;
    }

    bool continua;
    //std::cout <<" na la nb lb nc lc nd ld "<< n_a << l_a << n_b << l_b << n_c << l_c << n_d << l_d << lambda<<" continua"<<std::endl;
    //std::cin >> continua;

    /// Integral 8.42 Suhonen
    /// int(r2^2dr2 int(r1^2dr1 g_nala(r1) g_nblb(r2) v_lambda(r1,r2) g_nclc(r1) g_ndld(r2) )0-inf)0-inf

     // Useful constant and auxiliary indexes
    double R_aux = 0.;
    double Aux_n_ac, Aux_n_bd, Aux_k, Aux_i, Aux_j;
    double Aux_a, Aux_b;
    int alpha, beta, n;

    double gamma = pow(1/ b_lenght(A), 2) + pow(1/ R_interaction, 2);

    Fraction * Num_k;
    Fraction * Den_k;
    Num_k = new Fraction[1];
    Den_k = new Fraction[3];
    int aux_NULL = NULL;

    // sum i_x for Laguerre decomposition
    for(int i_a = 0; i_a <= n_a ; i_a ++ ){
        Num_k[0] = Fraction(1,2) + n_a + l_a;

        Den_k[0] = i_a;
        Den_k[1] = Fraction(1,2) + i_a + l_a;
        Den_k[2] = n_a - i_a;

        Aux_a = factorial_function(Num_k,Den_k,1,3);
        for(int i_c = 0; i_c <= n_c; i_c ++){

            Num_k[0] = Fraction(1,2) + n_c + l_c;

            Den_k[0] = i_c;
            Den_k[1] = Fraction(1,2) + i_c + l_c;
            Den_k[2] = n_c - i_c;

            Aux_n_ac = pow((-1)/pow(b_lenght(A),2) , i_a + i_c) * exp(Aux_a + factorial_function(Num_k,Den_k,1,3));

            //std:: cout << "i_ac= " << i_a << i_c <<"   ================= "<<std:: endl;
            //std:: cout << "*---------------------------*"<<std::endl;
            //std::cin >> continua;
            for(int i_b = 0; i_b <= n_b; i_b ++){
                Num_k[0] = Fraction(1,2) + n_b + l_b;

                Den_k[0] = i_b;
                Den_k[1] = Fraction(1,2) + i_b + l_b;
                Den_k[2] = n_b - i_b;

                Aux_b = factorial_function(Num_k,Den_k,1,3);
                for(int i_d = 0; i_d <= n_d; i_d ++){
                    // Reusing of the arrays

                    Num_k[0] = Fraction(1,2) + n_d + l_d;


                    Den_k[0] = i_d;
                    Den_k[1] = Fraction(1,2) + i_d + l_d;
                    Den_k[2] = n_d - i_d;

                    Aux_n_bd = pow((-1)/pow(b_lenght(A),2) , i_b + i_d) * exp(Aux_b + factorial_function(Num_k,Den_k,1,3));

                    //std:: cout << "Aux_n_bd= " << Aux_n_bd <<std:: endl;
                    //std:: cout << "*---------------------------*"<<std::endl;

                    // sum k for Legendre decomposition
                    for(int k = 0; k <= floor(lambda/2); k++){
                        // Reusing parts of the arrays
                        Num_k[0] = Fraction(2*(lambda - k));
                        Den_k[0] = Fraction(k);
                        Den_k[1] = Fraction(lambda - k);

                        Aux_k = pow(-1,k) * exp(factorial_function(Num_k,Den_k,1,2));

                        // sum i for x^n exp(-ax) integral
                        for(int i = 0; i <= (lambda - 2*k); i++){

                            Aux_i = pow((-1)* pow(R_interaction, 2)/2,i) / exp(factorial(Fraction(lambda - 2*k - i)));

                            // habia sumados en alpha y beta +2 donde pone 0
                            alpha = 2 + l_a + l_c + 0 + 2*(i_a + i_c) - 1 - i;
                            beta = 2 + l_b + l_d + 0 + 2*(i_b + i_d) - 1 - i;
                            n = lambda - 2*k - i + 1;
                            // and sum in j for binomial formula in the integrals
                            for(int j = 0; j <= beta; j ++){
                                // Reusing part of the arrays
                                Num_k[0] = Fraction(beta);
                                Den_k[0] = Fraction(j);
                                Den_k[1] = Fraction(beta - j);

                                if((n + beta - j + 1)%2 == 1){
                                    //std::cout<<" AUX J = 0"<< std::endl;
                                    Aux_j = 0;
                                }
                                else{
                                    Aux_j = exp(
                                                factorial_function(Num_k,Den_k,1,2) + gamma_function(Fraction(j + 1,2),aux_NULL)
                                                + gamma_function(Fraction(alpha + 2*(beta - j) + 1,2),aux_NULL)
                                                )/
                                            (2 * pow(gamma, beta + (1-j)/2) * pow(R_interaction, 2*(beta - j)) *
                                                pow(gamma - (1/gamma*pow(R_interaction,4)),(alpha + 2*(beta - j) + 1)/2)
                                             );

                                    /*
                                    std::cout<< "------------------"<<std::endl;
                                    std::cout<< " N 1::   "<<exp(
                                                factorial_function(Num_k,Den_k,1,2) +
                                                gamma_function(Fraction(j + 1,2),aux_NULL) + gamma_function(Fraction(alpha + 2*(beta - j) + 1,2),aux_NULL)
                                                )<<std::endl;
                                    std::cout<< " D_1::   "<<(2 * pow(gamma, beta + (1-j)/2))* pow(R_interaction, 2*(beta - j)) <<std::endl;
                                    std::cout<< " D_2::   "<< pow(gamma - (1/gamma*pow(R_interaction,4)),(alpha + 2*(beta - j)+ 1 )/2) <<std::endl;
                                    std::cout<< "------------------"<<std::endl;
                                    std::cout<< "ia_c,b_d= "<<i_a<<" "<<i_c<<" "<<i_b<<" "<<i_d<<"  // i,j,k,lam= "<<i<<" "<<j<<" "<<k<<" "<<lambda<<std::endl;
                                    std::cout<< "n_ac: " <<Aux_n_ac<<" /n_bd: "<<Aux_n_bd<<" /Aux k: "<<Aux_k<<" /Aux i: "<<Aux_i<<" /Aux j: "<<Aux_j<<std::endl;
                                    std::cout<< " Termino total = " << Aux_n_ac * Aux_n_bd * Aux_k * Aux_i * Aux_j << std::endl;
                                    std::cout<< " R_aux cumulativo = "<< R_aux + (Aux_n_ac * Aux_n_bd * Aux_k * Aux_i * Aux_j)<<std::endl;
                                    std::cout<< std::endl;
                                    * /

                                }
                                R_aux += Aux_n_ac * Aux_n_bd * Aux_k * Aux_i * Aux_j;
                                //std::cin>>continua;
                            }

                        }
                    }

                }
            }
        }
    }

    delete [] Num_k;
    delete [] Den_k;


    R_aux *= 2* M_PI * V_0* pow(R_interaction,2) / pow(2, lambda + 1);

    R_aux *= SHO_normalization(n_a, l_a)* SHO_normalization(n_b, l_b)* SHO_normalization(n_c, l_c)* SHO_normalization(n_d, l_d);

    //std::cout<< " Normalizacion todo = " << SHO_normalization(n_a, l_a)* SHO_normalization(n_b, l_b)* SHO_normalization(n_c, l_c)* SHO_normalization(n_d, l_d) <<
    //             "  ////  Constante fisica =  " <<  2* M_PI * V_0* pow(R_interaction,2) / pow(2, lambda + 1) << std::endl;
    //std::cout<< std::endl;
    //std::cin >> continua;
    if(R_aux > 1e2){
        std::cout <<" ---------- RESULTADO ----------- "<<std::endl;
        std::cout <<" na la nb lb nc lc nd ld "<< n_a << l_a << n_b << l_b << n_c << l_c << n_d << l_d << " "<< lambda<<" continua"<<std::endl;
        std::cout << " Error: R_aux = "<< R_aux <<std::endl;
        //std::cin >> continua;
    }
    return R_aux;
}
*/

double Radial_SDI(const QN_2body_radial &WF){

    // This expresion is normalized to sqrt(-V_0),
    // because of this, squaring lead an attractive potential
    return R_radius(A) * SHO_Radial(WF.n1,WF.l1,R_radius(A), b_lenght(A)) *
                         SHO_Radial(WF.n2,WF.l2,R_radius(A), b_lenght(A));
}
