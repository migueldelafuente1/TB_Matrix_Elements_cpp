#include "../include/Fractions.h"
#include "../include/Angular_Functions.h"
#include "../include/BM_Brackets.h"
#include "../include/Integrals.h"
#include "../include/Index_Coefficients.h"
#define _USE_MATH_DEFINES
#include <cmath>

using namespace std;
/// Functions and indexes from Marcos Monshinsky paper
/*
double R_M_E_Potential(int n, int l, int n_q, int l_q, double b_param){
    //cout << "RME(n,l,n',l',p)" << n<< " "<< l<<" "<<n_q<<" "<<l_q<<endl;
    // Expansion in Talmi integrals
    Fraction p;
    double Sum = 0.;
    int rho = 2*n + 2*n_q + l + l_q;

    for (p = Fraction((l+l_q),2); p <= (Fraction((l+l_q),2)+ n + n_q); p+=1){
        // B_coefficient contains the normalization (also in b_param)
	if((int(p) >= 0) && (int(p) <= rho)){
	    
	    Sum += Talmi_Integral(p, b_param) * B_coefficient(n,l,n_q,l_q,p, b_param);	
	    
	    std::cout << "I_" << p << " * " << B_coefficient(n,l,n_q,l_q,p, b_param)<< std::endl;
	}
        else{std::cout <<" out"<<std::endl;}
    }

    //std::cout << "Sum =" << Sum << std::endl;
    return Sum;
}
//*/

double R_M_E_Potential(int n, int l, int n_q, int l_q, int option,  double b_param){
    //cout << "RME(n,l,n',l')=" << n<< " "<< l<<" "<<n_q<<" "<<l_q<<endl;
    // Expansion in Talmi integrals
    int p;
    double Sum = 0.;
    int rho = 2*n + 2*n_q + l + l_q;

    for (p = ((l+l_q)/2); p <= (((l+l_q)/2)+ n + n_q); p++){
        // B_coefficient contains the normalization (also in b_param)
	if((p >= 0) && (p <= rho)){
	    
	    Sum += Talmi_Integral(p, option, b_param) * B_coefficient(n,l,n_q,l_q,p, b_param);	
	    
	    //std::cout << "I_" << p <<"(= "<< Talmi_Integral(p, b_param) << ") * B(= " << B_coefficient(n,l,n_q,l_q,p, b_param) <<")"<< std::endl;
	}
        else{std::cout <<" out"<<std::endl;}
    }

    //std::cout << "Sum =" << Sum << std::endl;
    return Sum;
}
///*
double M_E_Central(int n1, int l1, int n2, int l2, int lambda, int n1_q, int l1_q, int n2_q, int l2_q, int lambda_q, int option, double b_param){
    
    if(lambda != lambda_q){
        return 0;
    }

    int aux = l1_q + l2_q - l1 - l2;
    if(aux%2 != 0){
        return 0; // according to Moshisky's papper, aux must be even
    }

    else{
        int rho = (2*n1 + l1 + 2*n2 + l2);
	int rho_q = (2*n1_q + l1_q + 2*n2_q + l2_q);
	
        int l, n_q;
        double sum = 0.;
        for(int N = 0; N <= floor(rho/2); N++ ){
            for(int n = 0; n <= floor(rho/2) - N ; n++ ){
                for(int L = 0; L <= rho - 2*(n + N); L++ ){
                    l = (rho - L - 2*(N + n));
		    
		    //std::cout << n <<" "<< l <<" "<< N <<" "<< L << std::endl;
                    //cout << "triangular_condition(l,L,lambda)="<< l<<L<<lambda<< "="<<(triangular_condition(l,L,lambda))<< endl;
		    if(rho == (2*(N+l) + L + l)){
			if(l < 0){}
			
			else if(triangular_condition(l,L,lambda)){
			    //n_q = n1_q + n2_q - n1 - n2 +(aux/2) + n;
			    n_q = ((rho - rho_q)/2) + n;
			    
			    double B1 = BMB(n,l,N,L,lambda,n1,l1,n2,l2) ;
			    if(fabs(B1) > 1.e-12){
				//double B2 = BMB((n1_q + n2_q - n1 - n2 +(aux/2) + n), l,N,L,lambda,n1_q,l1_q,n2_q,l2_q);
				
				double B2 = BMB(n_q, l,N,L,lambda,n1_q,l1_q,n2_q,l2_q);
				if(fabs(B2) > 1.e-12){
				    //sum += B1 * B2 * R_M_E_Potential(n,l,(n1_q + n2_q - n1 - n2 +(aux/2) + n),l, b_param);
				    //std::cout << "(n l N L)=" << n << ","<< l << ","<< N << ","<< L << "  ; " <<
				    //	 "(n' l N L)=" << n_q << ","<< l << ","<< N << ","<< L << "  ; " << std::endl; 
				    //std::cout << "BMB: "<< B1<< "    BMB_q: "<< B2 << std::endl;
				    sum += B1 * B2 * R_M_E_Potential(n,l,n_q,l,option, b_param);
				    
				}
			    }

			    // sum += BMB(n,l,N,L,lambda,n1,l1,n2,l2) *
			    //     BMB((n1_q + n2_q - n1 - n2 +(aux/2) + n), l,N,L,lambda,n1_q,l1_q,n2_q,l2_q) *
			    //     R_M_E_Potential(n,l,(n1_q + n2_q - n1 - n2 +(aux/2) + n),l);
			}
		    }

                }
            }
        }
        return sum;
    }
}
//*/

/* ////  THE MOSHINSKY METHOD ALL AT ONCE FOR CENTRAL ALL AT ONCE
double M_E_Central(int n1, int l1, int n2, int l2, int lambda, int n1_q, int l1_q, int n2_q, int l2_q, int lambda_q, double b_param){
    if(lambda != lambda_q){
        return 0;
    }

    int aux = l1_q + l2_q - l1 - l2;
    if(aux%2 != 0){
        return 0; // according to Moshisky's papper, aux must be even
    }

    else{
	// DIAGONAL CASE
        if ((((n1 == n1_q) && (n2 == n2_q) )&& (l1 == l1_q)) && (l2 == l2_q)){
	    
	    int rho = (2*n1 + l1 + 2*n2 + l2);
	    int l;
	    double sum = 0.;
		
	    for(int N = 0; N <= floor(rho/2); N++ ){
		for(int n = 0; n <= floor(rho/2) - N ; n++ ){
		    for(int L = 0; L <= rho - 2*(n + N); L++ ){
			l = (rho - L - 2*(N + n));
			
			//cout << "triangular_condition(l,L,lambda)="<< l<<L<<lambda<< "="<<(triangular_condition(l,L,lambda))<< endl;
			if( l >= 0){
			    
			    if(rho == (2*(N+n) + L + l)){
				
				if(triangular_condition(l,L,lambda)){
				    
				    double B = BMB(n,l,N,L,lambda,n1,l1,n2,l2) ;

					for(Fraction p = Fraction(l); p <= Fraction(2*n + l); p+=1 ){
					    sum += pow(B,2) * B_coefficient(n,l,n,l,p, b_param) * Talmi_Integral(p, b_param);
					    std::cout << "I_(=" <<Talmi_Integral(p, b_param) <<")" << p << " * " <<B_coefficient(n,l,n,l,p, b_param)<< std::endl;
					}   

				}
			    }
			}

		    }
		}
	      
	    }
	    return sum;
	}
	
	// OFF DIAGONAL CASE
	else{
	    int rho = (2*n1 + l1 + 2*n2 + l2);
	    int rho_q = (2*n1_q + l1_q + 2*n2_q + l2_q);

	    int l, n_q;
	    double sum = 0.;
	    
	    for(int N = 0; N <= floor(rho/2); N++ ){
		for(int n = 0; n <= floor(rho/2) - N ; n++ ){
		    for(int L = 0; L <= rho - 2*(n + N); L++ ){
			l = (rho - L - 2*(N + n));

			if(l >= 0){
			
			    if(triangular_condition(l,L,lambda)){
				n_q = ((rho_q - rho)/2) + n;
				 
				if((rho_q - rho)%2 != 0){std::cout << " n' is HI" << std::endl;}
				    
				double B = BMB(n,l,N,L,lambda,n1,l1,n2,l2) ;
				
				if(fabs(B) > 1e-20){
				    //double B2 = BMB((n1_q + n2_q - n1 - n2 +(aux/2) + n), l,N,L,lambda,n1_q,l1_q,n2_q,l2_q);
				    double B_q = BMB(n1_q,l1_q,n2_q,l2_q,lambda, n_q, l,N,L);
				    if(fabs(B_q) > 1e-20){

					for(Fraction p = Fraction(l); p <= Fraction(n + n_q + l); p+=1){
					    sum += B * B_q * B_coefficient(n,l,n_q,l,p, b_param) * Talmi_Integral(p, b_param);
					     std::cout << "I_(=" <<Talmi_Integral(p, b_param) <<")" << p  << " * " <<B_coefficient(n,l,n_q,l,p, b_param)<< std::endl;
					} 
				    }
				}
			      

			    }
			}
			

		    }
		}
	    }
	    return sum;
	}

	
    }
}
//*/

double M_E_Spin_Orbit(int n1, int l1, int n2, int l2, int lambda, int n1_q, int l1_q, int n2_q, int l2_q, int lambda_q, int option, double b_param ){

    int aux = l1_q + l2_q - l1 - l2;
    if((aux)%2 != 0){
        return 0; // according to Monshisky papper
    }
    else{
        int rho = (2*n1 + l1 + 2*n2 + l2);
        int l;
        double sum = 0.;
        for(int N=0; N<=floor(rho/2); N++){
            for(int n=floor(rho/2); n >= 0; n--){
                for(int L=0; L <= rho; L++){
                    l = (rho - L - 2*(N + n));

                    if(l < 0){}
                        //int n_q = n1_q + n2_q - n1 - n2 +(aux/2) + n;
                    else if(triangular_condition(l,L,lambda)){

                        double tol = 1e-20;
                        double Racah = Racah_Coefficient(l,l,lambda,lambda_q,1,L);
                        if ( fabs(Racah) > tol){

                            double B1 = BMB(n,l,N,L,lambda,n1,l1,n2,l2) ;
                            if(fabs(B1) > tol){

                                double B2 = BMB((n1_q + n2_q - n1 - n2 +(aux/2) + n), l,N,L,lambda_q,n1_q,l1_q,n2_q,l2_q);
                                if(fabs(B2) > tol){

                                    double RME = R_M_E_Potential(n,l,(n1_q + n2_q - n1 - n2 +(aux/2) + n),l,option, b_param);
                                    if(fabs(RME) > tol){
                                        sum += B1 * B2 * pow(-1,(L+1-l-lambda)) * RME * Racah *
                                               sqrt(l*(l+1)*(2*l+1));     // [lambda][lambda_q] are in the other function


                                    }
                                }
                            }
                        }

                        //sum += BMB(n,l,N,L,lambda,n1,l1,n2,l2) *
                        //    BMB((n1_q + n2_q - n1 - n2 +(aux/2) + n), l,N,L,lambda_q,n1_q,l1_q,n2_q,l2_q) *
                        //    R_M_E_Potential(n,l,(n1_q + n2_q - n1 - n2 +(aux/2) + n),l) *
                        //    pow(-1,(L+1-l-lambda))  *
                        //    Racah_Coefficient(l,l,lambda,lambda_q,1,L)*
                        //    sqrt(l*(l+1)*(2*l+1)*(2*lambda+1)*(2*lambda_q+1));
                    }
                }
            }
        }
        return sum;
    }
}


double M_E_Spin_Orbit_jj( const QN_2body_jj_Coupling & BRA,const QN_2body_jj_Coupling & KET , int option,double b_param){

    // Passing of the parameters
    int n1 = BRA.n1;
    int l1 = BRA.l1;
    Fraction j1 = BRA.j1;
    int n2 = BRA.n2;
    int l2 = BRA.l2;
    Fraction j2 = BRA.j2;
    Fraction J = BRA.J;
    Fraction M = BRA.M;
    int T = BRA.T;

    int n1_q = KET.n1;
    int l1_q = KET.l1;
    Fraction j1_q = KET.j1;
    int n2_q = KET.n2;
    int l2_q = KET.l2;
    Fraction j2_q = KET.j2;
    Fraction J_q = KET.J;
    Fraction M_q = KET.M;
    int T_q = KET.T;

    // Starting of computation
    if(M != M_q){return 0;}
    if(J != J_q){return 0;}
    if(T != T_q){return 0;}

    // Isospin and total coupled angular momentum must be odd:
    if((J + T)%2 != 1){return 0;}

    // Sum over S,lambda LS-jj Coupling Coefficients (reduced to 9j symbols)
    double Sum = 0.;

    Fraction S,S_q;
    int lambda,lambda_q;
    for (S = Fraction(0); S<=1; S +=1){
        for(S_q = Fraction(0); S_q <=1; S_q+=1){
            if((S_q != S)||(S != 1)){}
            else{
                for(lambda = abs(l1 - l2); lambda <= l1 + l2; lambda ++){
                    //for(lambda_q = abs(l1_q - l2_q); lambda_q <= (l1_q + l2_q); lambda_q++){ // Global case
                    for(lambda_q = lambda - 1; lambda_q <= lambda+1; lambda_q ++){              // Restriction, A.17
                        //cout<< "S, S_q, Lambda, Lambda_q="<<S <<" " <<S_q<<" "<< lambda<< " "<<lambda_q<<endl;
                        if(lambda_q >= 0){

                            double tol = 1e-20;
                            double N1 = LS_jj_coupling_Coeff(Fraction(l1),Fraction(1,2),j1,Fraction(l2),Fraction(1,2),j2,Fraction(lambda),S,J);
                            if(fabs(N1) > tol){

                                double N2 = LS_jj_coupling_Coeff(Fraction(l1_q),Fraction(1,2),j1_q,Fraction(l2_q),Fraction(1,2),j2_q,Fraction(lambda_q),S_q,J_q);
                                if(fabs(N2) > tol){

                                    double MELS = M_E_Spin_Orbit(n1, l1, n2, l2, lambda, n1_q, l1_q, n2_q, l2_q, lambda_q,option,b_param);
                                    if(fabs(MELS) > tol){
                                        Sum += sqrt((((2*j1+ 1)*(2*j2 + 1)*(2*S+1)) *
                                             ((2*j1_q+ 1)*(2*j2_q + 1)*(2*S_q+1)))) *
                                             (2*lambda+1) *(2*lambda_q+1)*
                                             N1 * N2 * MELS;
                                        //cout<<"Suma externa ="<< Sum<<endl;
                                        //bool pause;
                                        //cin >> pause;
                                    }
                                }
                            }
                            //Sum += sqrt((((2*j1+ 1)*(2*j2 + 1)*(2*S+1)*(2*lambda+1)) *
                            //         ((2*j1_q+ 1)*(2*j2_q + 1)*(2*S_q+1)*(2*lambda_q+1)))) *
                            //         Nine_j_Symbol(Fraction(l1),Fraction(1,2),j1,Fraction(l2),Fraction(1,2),j2,Fraction(lambda),S,J) *
                            //         Nine_j_Symbol(Fraction(l1_q),Fraction(1,2),j1_q,Fraction(l2_q),Fraction(1,2),j2_q,Fraction(lambda_q),S_q,J_q) *
                            //         M_E_Spin_Orbit(n1, l1, n2, l2, lambda, n1_q, l1_q, n2_q, l2_q, lambda_q);
                        }

                    }
                }

            }
        }
    }

    return Sum * H_BAR;

}

double M_E_Tensor(int n1, int l1, int n2, int l2, int lambda, int n1_q, int l1_q, int n2_q, int l2_q, int lambda_q, int option,double b_param){

    //int aux = l1_q + l2_q - l1 - l2;
    // This is not assumed explicitly in the paper. But aux is always even
    /*if((aux)%2 != 0){
        return 0; // considering Monshisky papper, aux must be even
    }
    else{*/

    int rho = (2*n1 + l1 + 2*n2 + l2);
    int l,n_q;
    double sum = 0.;
    for(int N=0; N<=floor(rho/2); N++){
        for(int n=floor(rho/2); n >= 0; n--){
            for(int L=0; L <= rho; L++){
                l = (rho - L - 2*(N + n));
                //cout << "triangular_condition(l,L,lambda)="<< l<<L<<lambda<< "="<<(triangular_condition(l,L,lambda))<< endl;
                if(l < 0){}
                else if(triangular_condition(l,L,lambda)){
                    for(int l_q=l-2; l_q <=l+2; l_q += 2){

                        if(l_q >= 0){
                            //n_q = n + (n1_q + n2_q - n1 - n2)+((l1_q + l2_q - l1 - l2)/2) + ((l - l_q)/2);
                            int aux = l1_q + l2_q - l1 - l2;
                            n_q = n + (n1_q + n2_q - n1 - n2)+(aux/2) +((l - l_q)/2) ;
                            /*if((l-l_q % 2 == 1)||(aux%2==1)){
                                cout<<" Necesaria fraccion, l-l_q="<< l-l_q<<" , aux ="<<aux<<endl;
                                bool pause;
                                cin>>pause;
                            }
                            if((l-l_q % 2 == 1)&&(aux%2==1)){
                                cout<<" Entero"<<endl;
                            }*/
                            //cout<<"l,L,lamnda"<<l<<L<<lambda<<" l_q="<< l_q<<endl;

                            double tol = 1e-20;
                            double Racah = Racah_Coefficient(l,l_q,lambda,lambda_q,2,L);
                            //cout<<"Racah="<<Racah<<endl;
                            if(fabs(Racah) > tol){

                                double CG = Clebsh_Gordan(l,2,l_q,0,0,0);
                                //cout<< "CG="<<CG<<endl;
                                if(fabs(CG) > tol){

                                    double B1 = BMB(n,l,N,L,lambda,n1,l1,n2,l2);
                                    //cout<< "B1="<<B1<<endl;
                                    if(fabs(B1) > tol){

                                        double B2 = BMB(n_q, l_q,N,L,lambda_q,n1_q,l1_q,n2_q,l2_q);
                                        //cout<< "B2="<<B2<<endl;
                                        if(fabs(B2) > tol){

                                            double RME = R_M_E_Potential(n,l,n_q,l_q,option,b_param);
                                            //cout<< "RME="<<RME<<endl;
                                            if(fabs(RME) > tol){
                                                sum += B1* B2 * RME * Racah * CG * pow(-1,(L-l_q-lambda)) *
                                                        sqrt((5*M_PI/4)*(2*l+1));   // [lambda][lambda_q] are in the other function

                                                //cout<<"Suma CT="<<sum<<endl;

                                            }
                                        }
                                    }
                                }
                            }

                            //sum += BMB(n,l,N,L,lambda,n1,l1,n2,l2) *
                            //    BMB(n_q, l_q,N,L,lambda_q,n1_q,l1_q,n2_q,l2_q) *
                            //    R_M_E_Potential(n,l,n_q,l_q) *
                            //    pow(-1,(L-l_q-lambda)) * Racah_Coefficient(l,l_q,lambda,lambda_q,2,L)*
                            //    Clebsh_Gordan(l,2,l_q,0,0,0)*sqrt((5*M_PI/4)*(2*l+1)*(2*lambda+1)*(2*lambda_q+1));
                        }

                    }
                }
            }
        }
    }
    return sum;
    // }
}


double M_E_Tensor_jj(const QN_2body_jj_Coupling &BRA, const QN_2body_jj_Coupling &KET,int option, double b_param){

    // Passing of the parameters
    int n1 = BRA.n1;
    int l1 = BRA.l1;
    Fraction j1 = BRA.j1;
    int n2 = BRA.n2;
    int l2 = BRA.l2;
    Fraction j2 = BRA.j2;
    Fraction J = BRA.J;
    Fraction M = BRA.M;
    int T =  BRA.T;

    int n1_q = KET.n1;
    int l1_q = KET.l1;
    Fraction j1_q = KET.j1;
    int n2_q = KET.n2;
    int l2_q = KET.l2;
    Fraction j2_q = KET.j2;
    Fraction J_q = KET.J;
    Fraction M_q = KET.M;
    int T_q = KET.T;

    // Starting of computation
    if(M != M_q){return 0;}
    if(J != J_q){return 0;}
    if(T != T_q){return 0;}

    // Isospin and total coupled angular momentum must be odd:
    if((J + T)%2 != 1){return 0;}

    // Sum over S,lambda LS-jj Coupling Coefficients (reduced to 9j symbols)
    double Sum = 0.;

    Fraction S,S_q;
    int lambda,lambda_q;
    for (S = Fraction(0); S<=1; S +=1){
        for(S_q = Fraction(0); S_q <=1; S_q+=1){
            if((S_q != S)||(S != 1)){}
            else{
                for(lambda = abs(l1 - l2); lambda <= l1 + l2; lambda ++){
                    //for(lambda_q = abs(l1_q - l2_q); lambda_q <= (l1_q + l2_q); lambda_q++){ // Global case
                    for(lambda_q = lambda - 1; lambda_q <= lambda+1; lambda_q ++){              // Restriction, A.17

                        //cout<< "S, S_q, Lambda, Lambda_q="<<S <<" " <<S_q<<" "<< lambda<< " "<<lambda_q<<endl;
                        if(lambda_q >= 0){

                            double tol = 1e-20;

                            double N1 = LS_jj_coupling_Coeff(Fraction(l1),Fraction(1,2),j1,Fraction(l2),Fraction(1,2),j2,Fraction(lambda),S,J);
                            if(fabs(N1) > tol){

                                double N2 = LS_jj_coupling_Coeff(Fraction(l1_q),Fraction(1,2),j1_q,Fraction(l2_q),Fraction(1,2),j2_q,Fraction(lambda_q),S_q,J_q);
                                if(fabs(N2) > tol){
                                    //cout<< " Hola elemento de matriz"<<endl;
                                    double METEN = M_E_Tensor(n1, l1, n2, l2, lambda, n1_q, l1_q, n2_q, l2_q, lambda_q,option, b_param);
                                    //cout<<"METEN = "<<METEN<<endl;
                                    if(fabs(METEN) > tol){
                                        Sum += sqrt((((2*j1+ 1)*(2*j2 + 1)*(2*S+1))*((2*j1_q+ 1)*(2*j2_q + 1)*(2*S_q+1))))
                                                *(2*lambda+1)*(2*lambda_q+1) * N1 * N2 * METEN;
                                        //cout<<"Suma externa ="<<Sum<<endl;
                                        //bool pause;
                                        //cin >> pause;
                                    }
                                }
                            }


                        }

                        //Sum += sqrt((((2*j1+ 1)*(2*j2 + 1)*(2*S+1)*(2*lambda+1)) *
                        //             ((2*j1_q+ 1)*(2*j2_q + 1)*(2*S_q+1)*(2*lambda_q+1)))) *
                        //             Nine_j_Symbol(Fraction(l1),Fraction(1,2),j1,Fraction(l2),Fraction(1,2),j2,Fraction(lambda),S,J) *
                        //             Nine_j_Symbol(Fraction(l1_q),Fraction(1,2),j1_q,Fraction(l2_q),Fraction(1,2),j2_q,Fraction(lambda_q),S_q,J_q) *
                        //             M_E_Tensor(n1, l1, n2, l2, lambda, n1_q, l1_q, n2_q, l2_q, lambda_q);
                    }
                }

            }
        }
    }

    return Sum * pow(H_BAR,2);

}
