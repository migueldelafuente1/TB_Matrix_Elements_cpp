#include "../include/Angular_Functions.h"
#include "../include/Wave_Functions.h"
#include "../include/factorials.h"
#include "../include/Fractions.h"
#include "../include/BM_Brackets.h"
#include "../include/Index_Coefficients.h"
#include "../include/Integrals.h"
#include "../include/Nuclear_Matrix_Elements_Art.h"
#include "../include/Nuclear_Matrix_Elements_M.h"
#include "../include/Nuclear_Matrix_Elements_Suh.h"
#include "../include/Nuclear_Matrix_Elements_BB.h"
#include "../include/output.h"

#include <cmath>
#include <math.h>
#include <fstream>
#include <iostream>
#include <iomanip>
#include <stdlib.h>
#include <time.h>
#include <complex>
using namespace std;


void check_factorial_function(){
    // modify the values of dimension and the factorials introduced
    int * Denominador;
    Denominador = new int[2];
    Denominador[0]=2;
    Denominador[1]=3;

    int * Numerador;
    Numerador = new int[3];
    Numerador[0]=1;
    Numerador[1]=1;
    Numerador[2]=10;

    double b = factorial_function(Numerador,Denominador,3,2);
    std::cout << "factorial_function (int) =" << b << " / exp(b)="<< exp(b) << std::endl;

    Fraction * Den ;
    Den = new Fraction[2];
    Den[0]=Fraction(2);
    Den[1]=Fraction(3);

    Fraction * Num ;
    Num = new Fraction[3];
    Num[0]=Fraction(1);
    Num[1]=Fraction(1);
    Num[2]=Fraction(10);

    b = factorial_function(Num,Den,3,2);
    std::cout << "factorial_function (frac)=" << b << " / exp(b)="<< exp(b) << std::endl;

}

// Realistic checking of non zero direct formula or Recurrence Functions
// Repeat a millon times and divide the time.
//double aux;
//for(int i = 0; i<=1e+06; i++){
//
//    if(i%1000==0){cout << i << endl;}
//    aux = BMB(3,2,3,4,4,2,4,2,6);
// }


//check_Racah1(1,2,0,12);
void check_Racah1(int b, int d,int f_min,int f_max){
    // W(abcd,ef) coefficients have analitical solution for certain values, which are the followings:
    int e=1;
    double tol = 1e-5;
    double analitic, numerical;
    int errors = 0;
    int total = 0;

    std::ofstream f_py;
    f_py.open("data/racah_for_python.txt",std::fstream::out);
    std::ofstream f_salida;
    f_salida.open("data/Racah_cheking.txt",std::fstream::out);
    // I save the numerical value obtained and, in the next line, the code to check it out from Wolfram Alpha.

    std::cout << "a / b / c / d / e / f // W_coeff() // W_an() //   Error  "<< std::endl;
    std::cout << "______________________//___________//________//__________"<< std::endl;
    for(int f=f_min; f<=f_max; f++){
        std::cout << f <<std::endl;
        for(int c=d-1; c<=d+1; c++){
            std::cout << "______________________//___________//________//__________"<< std::endl;
            for(int a=b-1; a<=b+1; a++){
                analitic = Racah_Coeficient_Analitic( a, b, c, d, 1, f);
                numerical = Racah_Coefficient( a, b, c, d, 1, f);

                f_salida <<"======================================================"<<std::endl;
                f_salida <<"W("<<a<<" / "<< b <<" / "<<c<<" / "<<d<<" / "<<e<<" / "<<f<<")="<< numerical << std::endl;
                f_salida <<"SixJSymbol[{"<< a <<","<< b <<","<< e<< "}, {"<< d <<","<< c <<","<< f <<"}]"<< std::endl;

                f_py << a<<" "<< b <<" "<<c<<" "<<d<<" "<<e<<" "<<f<<" "<< numerical<<std::endl;
                // for checking use  SixJSymbol[{j1, j2, j3}, {j4, j5, j6}] in Wolfram Alpha
                if(fabs(analitic-numerical)>tol){
                    std::cout <<a<<" / "<<b<<" / "<<c<<" / "<<d<<" / "<<e<<" / "<<f<<" // "<<numerical<<" // "<<analitic<<" //   ERROR  "<<std::endl;
                    errors++;
                    f_salida<< "Analitical value "<< analitic<<std::endl;
                    f_salida<< "Was Erroneous"<<std::endl;
                }
                else{
                    std::cout <<a<<" / "<<b<<" / "<<c<<" / "<<d<<" / "<<e<<" / "<<f<<" // "<<numerical<<" // "<<analitic<<" //          "<<std::endl;
                }
                total++;
            }
        }
    }
    std::cout << "______________________//___________//________//__________"<< std::endl;
    std::cout << "Number of Errors:" <<errors<<" of "<< total <<std::endl;

    f_salida.close();
    f_py.close();
}



//double R = Racah_Coefficient(0,2,0,2,2,0);
//double R = Racah_Coefficient(0,4,7,3,4,7);
//cout <<"W()="<< R <<endl;

//check_Racah2(4);  // for N=10, there are 4,100,000 coefficients, execution time = 63 s
//check_Racah(20,with N/2)  are  did in  s of time ( minutes)

void check_Racah2(int N, bool half_integers){
    int i = 0;
    if(half_integers){
        std::ofstream f_py;
        f_py.open("data/racah_hf_for_python.txt",std::fstream::out);

        //int N_max = 3*pow(N,6)*pow(N+1,3)*(3*N+1)*pow((2*N+1),2)/pow(2,4);
        //std::cout<< "Length = "<<N_max<<std::endl;
        double numerical;
        for(Fraction a=0; a<=N; a+=Fraction(1,2)){
            for(Fraction b=0; b<=N; b+=Fraction(1,2)){
                for(Fraction d=0; d<=N; d+=Fraction(1,2)){
                    for(Fraction c=0; c<=a+b+d+1; c+=Fraction(1,2)){
                        for(Fraction e=0; e<=a+b+1; e+=Fraction(1,2)){
                            for(Fraction f=0; f<=b+d+1; f+=Fraction(1,2)){
                                numerical = Racah_Coefficient( a, b, c, d, e, f);
                                f_py << a<<" "<< b <<" "<<c<<" "<<d<<" "<<e<<" "<<f<<" "<< numerical<<std::endl;
                                //std::cout <<a<<" / "<<b<<" / "<<c<<" / "<<d<<" / "<<e<<" / "<<f<<" // "<<numerical<<std::endl;
                                i++;
                            }
                        }
                    }
                }
            }
        }
        f_py.close();

    }


    else{
        std::ofstream f_py;
        f_py.open("data/racah_for_python.txt",std::fstream::out);

        //int N_max = 3*pow(N,6)*pow(N+1,3)*(3*N+1)*pow((2*N+1),2)/pow(2,4);
        //std::cout<< "Length = "<<N_max<<std::endl;
        double numerical;
        for(int a=0; a<=N; a++){
            for(int b=0; b<=N; b++){
                for(int d=0; d<=N; d++){
                    for(int c=0; c<=a+b+d+1; c++){
                        for(int e=0; e<=a+b+1; e++){
                            for(int f=0; f<=b+d+1; f++){
                                numerical = Racah_Coefficient( a, b, c, d, e, f);
                                f_py << a<<" "<< b <<" "<<c<<" "<<d<<" "<<e<<" "<<f<<" "<< numerical<<std::endl;
                                //std::cout <<a<<" / "<<b<<" / "<<c<<" / "<<d<<" / "<<e<<" / "<<f<<" // "<<numerical<<std::endl;
                                i++;
                            }
                        }
                    }
                }
            }
        }
        f_py.close();
    }
    std::cout<< i <<" Racah C. have been caculated"<<std::endl;
}

//double cg = Clebsh_Gordan(1,2,1,-2,-2,-4);
//cout <<"cg()="<< cg <<endl;
//check_CCGG(23);    // for 23, there are 4,580,000 coefficients, execution time = 137 s
//check_CCGG(20,with N/2) are 16,800,000 did in 498.4 s of time (8.3 minutes)

// for half integers, I get, N=20 and then
void check_CCGG(int N, bool with_halfintegers){
    int i=0;
    clock_t t;
    t = clock();
    if(with_halfintegers){

        std::fstream f_py;
        f_py.open("data/ccgg_hf_for_python.txt",ios::out);
        if (!f_py.is_open()){exit(0);}
        Fraction m,step;
        step = Fraction(1,2);
        double numerical;
        for(Fraction j1=0; j1<=N; j1+=step){
            for(Fraction j2=j1; j2<=N; j2+=step){
                for(Fraction j=abs(j1-j2); j<=(j1 + j2); j+=step){

                    for(Fraction m1=(-1)*j1; m1<=j1; m1+=1){
                        for(Fraction m2=(-1)*j2; m2<=j2; m2+=1){
                            m = (m1+m2);
                            numerical = Clebsh_Gordan(j1,j2,j,m1,m2,m);
                            i++;
                            f_py <<j1<<" "<<j2<<" "<<j<<" "<<m1<<" "<<m2<<" "<<m<<" "<< numerical<<std::endl;
                        }
                    }
                }
            }
        }

        f_py.close();
    }
    else{
        std::ofstream f_py;
        f_py.open("data/ccgg_for_python.txt",std::fstream::out);

        double numerical;
        for(int j1=0; j1<=N; j1++){
            for(int j2=j1; j2<=N; j2++){
                for(int j=abs(j1-j2); j<=(j1 + j2); j++){

                    for(int m1=-j1; m1<=j1; m1++){
                        for(int m2=-j2; m2<=j2; m2++){
                            numerical = Clebsh_Gordan(j1,j2,j,m1,m2,m1+m2);
                            f_py <<j1<<" "<<j2<<" "<<j<<" "<<m1<<" "<<m2<<" "<<m1+m2<<" "<< numerical<<std::endl;
                            i++;
                        }
                    }
                }
            }
        }

        f_py.close();
    }
    std::cout<< i <<" CG have been caculated"<<std::endl;
    std::cout<< (clock() - t) <<" seconds of processor time consumed by the program"<<std::endl;
}


//double numerical = Nine_j_Symbol(Fraction(1,1),Fraction(1,2),Fraction(1,2),Fraction(1,1),Fraction(1,2),Fraction(1,2),Fraction(1),Fraction(1),Fraction(1));
void check_9j(int N, bool with_halfintegers){
    int i=0;
    clock_t t;
    t = clock();

    if(with_halfintegers){
        std::fstream f_py;
        f_py.open("data/9j_for_python.txt",ios::out);
        if (!f_py.is_open()){exit(0);}
        Fraction step;
        step = Fraction(1,2);
        double numerical;

        for(Fraction s1=Fraction(-1,2); s1<=Fraction(1,2); s1+=1){
            for(Fraction l1=0; l1<=N; l1+=1){
                for(Fraction s2=Fraction(-1,2); s2<=Fraction(1,2); s2+=1){
                    for(Fraction l2=0; l2<=N; l2+=1){
                        for(Fraction j1=abs(l1-s1); j1<=(l1+s1); j1+=step){
                            for(Fraction j2=abs(l2-s2); j2<=(l2+s2); j2+=step){
                                for(Fraction J=abs(j1-j2); J<=(j1 + j2); J+=step){
                                    for(Fraction S=abs(s1-s2); S<=(s1 + s2); S+=step){
                                        for(Fraction L=abs(l1-l2); L<=(l1 + l2); L+=step){
                                            numerical = Nine_j_Symbol(l1,s1,j1,l2,s2,j2,L,S,J);
                                            f_py <<double(l1)<<" "<<double(s1)<<" "<<double(j1)<<" "
                                                <<double(l2)<<" "<<double(s2)<<" "<<double(j2)<<" "
                                                <<double(L)<<" "<<double(S)<<" "<<double(J)<<" "<< numerical<<std::endl;
                                            i++;
                                        }
                                    }
                                }
                            }
                        }
                    }
                }
            }
        }

        f_py.close();
    }

    else{
        std::ofstream f_py;
        f_py.open("data/9j_for_python.txt",std::fstream::out);

        double numerical;

        for(int s1=0; s1<=N; s1++){
            for(int l1=0; l1<=N; l1++){
                for(int s2=0; s2<=N; s2++){
                    for(int l2=0; l2<=N; l2++){
                        for(int j1=abs(l1-s1); j1<=(l1+s1); j1++){
                            for(int j2=abs(l2-s2); j2<=(l2+s2); j2++){
                                for(int J=abs(j1-j2); J<=(j1 + j2); J++){
                                    for(int S=abs(s1-s2); S<=(s1 + s2); S++){
                                        for(int L=abs(l1-l2); L<=(l1 + l2); L++){
                                            numerical = Nine_j_Symbol(l1,s1,j1,l2,s2,j2,L,S,J);
                                            f_py <<l1<<" "<<s1<<" "<<j1<<" "<<l2<<" "<<s2<<" "<<j2<<" "<<L<<" "<<S<<" "<<J<<" "<< numerical<<std::endl;
                                            i++;
                                        }
                                    }
                                }
                            }
                        }
                    }
                }
            }
        }

        f_py.close();
    }
    std::cout<< i <<" 9j symbols have been caculated"<<std::endl;
    std::cout<< (clock() - t) <<" seconds of processor time consumed by the program"<<std::endl;
}
// for N=5 there are 5,200,000 symbols and the computation requires 430 s

void check_LSjj_coupling_coeffs(QN_2body_jj_Coupling BRA, QN_2body_jj_Coupling KET){
    /*std::cout<< "-------------------------------------------------------------"<<std::endl;
    std::cout<< "  Check Devanathan's LSjj coupling coefficients              "<<std::endl;
    std::cout<< "-------------------------------------------------------------"<<std::endl;
    std::cout<< "  Check Orthonormality:  sum_LS = delta(j1,j1')delta(j2,j2') "<<std::endl;
    std::cout<< " for the same l and s QN.  Input j's BRA !=  j's KET  :      "<<std::endl;*/
    
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
    
    Fraction s_a = Fraction(1,2);
    Fraction s_b = s_a;
    /*
    double sum_LS_bra = 0.;
    double sum_LS_ket = 0.;
    
    for(int L = fabs(l_a-l_b); L <= (l_a+l_b); L++){
	for(int S = 0; S <=1; S++){
	    sum_LS_bra += LS_jj_coupling_Coeff(Fraction(l_a),s_a, j_a, Fraction(l_b), (s_b), j_b, Fraction (L), Fraction (S), J)*
			    LS_jj_coupling_Coeff(Fraction(l_a),s_a, j_a, Fraction(l_b), (s_b), j_b, Fraction (L), Fraction (S), J);
	}
    }
    
    for(int L = fabs(l_c-l_d); L <= (l_c+l_d); L++){
	for(int S = 0; S <=1; S++){
	    sum_LS_ket += LS_jj_coupling_Coeff(Fraction(l_c),s_a, j_c, Fraction(l_d), (s_b), j_d, Fraction (L), Fraction (S), J_q)*
			    LS_jj_coupling_Coeff(Fraction(l_c),s_a, j_c, Fraction(l_d), (s_b), j_d, Fraction (L), Fraction (S), J_q);
	}
    }
    std::cout <<" Same WF: sum_LS_bra= "<< sum_LS_bra <<" ,  sum_LS_ket= "<< sum_LS_ket<< std::endl;
    
    sum_LS_bra = 0.;
    for(int L = fabs(l_a-l_b); L <= (l_a+l_b); L++){
	for(int S = 0; S <=1; S++){
	    sum_LS_bra += LS_jj_coupling_Coeff(Fraction(l_a),s_a, j_a, Fraction(l_b), (s_b), j_b, Fraction (L), Fraction (S), J)*
			    LS_jj_coupling_Coeff(Fraction(l_a),s_a, j_c, Fraction(l_b), (s_b), j_d, Fraction (L), Fraction (S), J);
	}
    }
    std::cout <<" Diferent WF (ls and ss from Bra): sum_LS= "<< sum_LS_bra << std::endl;
    
    sum_LS_bra = 0.;  // L=L' S=S'
    sum_LS_ket = 0.;  // L!=L' 
    int L = l_a + l_b;
    int L_q = abs(l_a - l_b);
    for(Fraction j1 = abs(Fraction(l_a)-s_a); j1 <= (Fraction(l_a)+s_a); j1+=1){
	for(Fraction j2 = abs(Fraction(l_b)-s_b); j2 <= (Fraction(l_b)+s_b); j2+=1){
	    sum_LS_bra += LS_jj_coupling_Coeff(Fraction(l_a),s_a, j1, Fraction(l_b), (s_b), j2, Fraction (L), Fraction (1), J)*
			    LS_jj_coupling_Coeff(Fraction(l_a),s_a, j1, Fraction(l_b), (s_b), j2, Fraction (L), Fraction (1), J);
	    std::cout << " Sum_bra ="<< sum_LS_bra<<std::endl;
	}
    }
    std::cout <<" L=L & S=S': sum= "<< sum_LS_bra << std::endl;
    
    sum_LS_bra = 0.;  // L!=L' 
    sum_LS_ket = 0.;  // S!=S'
    for(Fraction j1 = abs(Fraction(l_a)-s_a); j1 <= (Fraction(l_a)+s_a); j1+=1){
	for(Fraction j2 = abs(Fraction(l_b)-s_b); j2 <= (Fraction(l_b)+s_b); j2+=1){
	    sum_LS_bra += LS_jj_coupling_Coeff(Fraction(l_a),s_a, j1, Fraction(l_b), (s_b), j2, Fraction (L), Fraction (1), J)*
			    LS_jj_coupling_Coeff(Fraction(l_a),s_a, j1, Fraction(l_b), (s_b), j2, Fraction (L), Fraction (0), J);
	    sum_LS_ket += LS_jj_coupling_Coeff(Fraction(l_a),s_a, j1, Fraction(l_b), (s_b), j2, Fraction (L), Fraction (0), J)*
			    LS_jj_coupling_Coeff(Fraction(l_a),s_a, j1, Fraction(l_b), (s_b), j2, Fraction (L_q), Fraction (0), J);
	}
    }
    std::cout <<" L != L': sum= "<< sum_LS_bra <<" ,  S != S': sum= "<< sum_LS_ket<< std::endl;
    std::cout<< std::endl;
    */
    std::cout<< "-----------------------------------------------------------------------------"<<std::endl;
    std::cout<< " Comparation of <BRA|KET> by a complete uncoupling and with LSjj (deltas)    "<<std::endl;
    
    // Complete Uncoupling
    Fraction m_b, m_d;
    int m_l_a, m_l_b, m_l_c, m_l_d;
    double unc_J_total, unc_j;   	
    double Result_C_U = 0.;

    for( Fraction m_a = Fraction(-1*j_a.numerator, j_a.denominator); m_a <= j_a ; m_a += 1){
        m_b = M - m_a;
        for ( Fraction m_c = Fraction(-1*j_c.numerator, j_c.denominator); m_c <= j_c ; m_c += 1){
            m_d = M - m_c;

            unc_J_total = Clebsh_Gordan(j_a,j_b,J, m_a, m_b, M) * Clebsh_Gordan(j_c,j_d,J, m_c, m_d, M);
	    //std::cout << "CG [ac]= " <<Clebsh_Gordan(j_a,j_b,J, m_a, m_b, M)<< " ,  CG [bd]= "<< Clebsh_Gordan(j_c,j_d,J, m_c, m_d, M)<< std::endl;
	    //unc_J_total = Clebsh_Gordan(j_a,j_b,J, m_a, m_b, M) * Clebsh_Gordan(j_c,j_d,J, m_c, m_d, M);
	    //std::cout << "Uncoupling J total= " << unc_J_total << std::endl;
	    /*
	    if((j_a == j_c) && (j_b == j_d)){
	      if(m_a == m_c){
		if(m_b == m_d){
		    Result_C_U += (unc_J_total);
		}
	      }
	    }*/
	    
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
				
				// Deltas for the uncoupled w.f.
				if((l_a == l_c) && (m_l_a == m_l_c)){
				  if((l_b == l_d) && (m_l_b == m_l_d)){
				    if(m_s_a == m_s_c){
				      if(m_s_b == m_s_d){
					  Result_C_U += (unc_j * unc_J_total);
				      }
				    }
				  }
				}

				
			    }
			}
		    }
		}
	    }
	}
    }
    //Result_C_U *= sqrt(2*int(J)+1);
    
    // LSjj uncoupling
    double Result_LSjj_U = 0.;
    
    for(int L1 = fabs(l_a-l_b); L1 <= (l_a+l_b); L1++){
	for(int S1 = 0; S1 <=1; S1++){
	    for(int L2 = fabs(l_c-l_d); L2 <= (l_c+l_d); L2++){
		for(int S2 = 0; S2 <=1; S2++){
		    
		  // Deltas for the uncoupled w.f.
		  if((l_a == l_c) && (l_b == l_d)){
		    if( (S1 == S2) && (L1 == L2)){
			Result_LSjj_U += LS_jj_coupling_Coeff(Fraction(l_a),s_a, j_a, Fraction(l_b), (s_b), j_b, Fraction (L1), Fraction (S1), J)*
					  LS_jj_coupling_Coeff(Fraction(l_c),s_a, j_c, Fraction(l_d), (s_b), j_d, Fraction (L2), Fraction (S2), J);
		    }
		  }
		     
		}
	    }
	}
    }
    std::cout <<" Result_C_U= "<< Result_C_U <<" ,  Result_LSjj_U= "<< Result_LSjj_U<< std::endl;
    /*
    std::cout<< "-----------------------------------------------------------------------------"<<std::endl;
    std::cout<< " Comparation of <BRA|BRA> by a complete uncoupling and with LSjj (deltas)    "<<std::endl;
    Result_C_U = 0.;

    for( Fraction m_a = j_a; m_a >= abs(j_a) ; m_a -= 1){
        m_b = M - m_a;
        for ( Fraction m_c = j_c; m_c >= abs(j_a) ; m_c -= 1){
            m_d = M - m_c; 

            unc_J_total = Clebsh_Gordan(j_a,j_b,J, m_a, m_b, M) * Clebsh_Gordan(j_a,j_b,J, m_c, m_d, M);
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
                                          * Clebsh_Gordan(Fraction(l_a),Fraction(1,2),j_a, Fraction(m_l_c),m_s_c, m_c )
                                            * Clebsh_Gordan(Fraction(l_b),Fraction(1,2),j_b, Fraction(m_l_d),m_s_d, m_d );
				//std::cout << std::endl;
				//std::cout << "Uncoupling j LS= " << unc_j << "  // Uncoupling J total= " << unc_J_total << std::endl;
				
				// Deltas for the uncoupled w.f.
				if((l_a == l_a) && (m_l_a == m_l_c)){
				  if((l_b == l_b) && (m_l_b == m_l_d)){
				    if(m_s_a == m_s_c){
				      if(m_s_b == m_s_d){
					  Result_C_U += (unc_j * unc_J_total);
				      }
				    }
				  }
				}

				
			    }
			}
		    }
		}
	    }
	}
    }
    Result_C_U  *= sqrt(2*J + 1);
    
    // LSjj uncoupling
    Result_LSjj_U = 0.;
    
    for(int L1 = fabs(l_a-l_a); L1 <= (l_b+l_b); L1++){
	for(int S1 = 0; S1 <=1; S1++){
	    for(int L2 = fabs(l_b-l_b); L2 <= (l_b+l_b); L2++){
		for(int S2 = 0; S2 <=1; S2++){
		    
		  // Deltas for the uncoupled w.f.
		  if((l_a == l_a) && (l_b == l_b)){
		    if( (S1 == S2) && (L1 == L2)){
			Result_LSjj_U += LS_jj_coupling_Coeff(Fraction(l_a),s_a, j_a, Fraction(l_b), (s_b), j_b, Fraction (L1), Fraction (S1), J)*
				    LS_jj_coupling_Coeff(Fraction(l_a),s_a, j_a, Fraction(l_b), (s_b), j_b, Fraction (L2), Fraction (S2), J);
		    }
		  }
		     
		}
	    }
	}
    }
    std::cout <<" Result_C_U= "<< Result_C_U <<" ,  Result_LSjj_U= "<< Result_LSjj_U<< std::endl;
    */
}



void check_BM_00(){
    // checking the values obtained in the reference (in order of appearance) with my number code
    std::cout<<"BMB_00("<<1306656<<")="<<BM_Bracket_00(1,3,0,6,6,5,6)<<" // correct="<<-0.29514<<endl;
    std::cout<<"BMB_00("<<1315766<<")="<<BM_Bracket_00(1,3,1,5,7,6,6)<<" // correct="<<0.18881<<endl;
    std::cout<<"BMB_00("<<0.000000<<")="<<BM_Bracket_00(0,0,0,0,0,0,0)<<" // correct="<<1.<<endl;
    std::cout<<"BMB_00("<<0.006066<<")="<<BM_Bracket_00(0,0,0,6,0,6,6)<<" // correct="<<0.<<endl;
    std::cout<<"BMB_00("<<0.001415<<")="<<BM_Bracket_00(0,5,0,1,4,1,5)<<" // correct="<<-0.11111<<endl; /**/
    std::cout<<"BMB_00("<<0.661534<<")="<<BM_Bracket_00(0,6,0,1,5,3,4)<<" // correct="<<0.028490<<endl;
    std::cout<<"BMB_00("<<0.404235<<")="<<BM_Bracket_00(0,4,0,4,2,3,5)<<" // correct="<<0.0<<endl;
    std::cout<<"BMB_00("<<0.10001055<<")="<<BM_Bracket_00(0,10,0,0,10,5,5)<<" // correct="<<-0.49607<<endl;
    std::cout<<"BMB_00("<<0.12001266<<")="<<BM_Bracket_00(0,12,0,0,12,6,6)<<" // correct="<<0.47495<<endl;
}

void BMB_symmetry_checking(int n, int l, int N, int L, int lambda, int n1, int l1, int n2, int l2){
    std::cout<<"Symmetry Properties// This values must be the same ========="<<std::endl;
    std::cout<<"<"<<n<<" "<<l<<" "<<N<<" "<<L<<" "<<lambda<<"  |   "<<n1<<"  "<<l1<<"  "<<n2<<"  "<<l2<<"> ::   "<<BMB(n,l,N,L,lambda,n1,l1,n2,l2)<<" = "<<std::endl;
    std::cout<<"<"<<n<<" "<<l<<" "<<N<<" "<<L<<" "<<lambda<<"  |   "<<n2<<"  "<<l2<<"  "<<n1<<"  "<<l1<<"> ::   "<<pow(-1,L-lambda)*BMB(n,l,N,L,lambda,n2,l2,n1,l1)<<" = "<<std::endl;
    std::cout<<"<"<<N<<" "<<L<<" "<<n<<" "<<l<<" "<<lambda<<"  |   "<<n1<<"  "<<l1<<"  "<<n2<<"  "<<l2<<"> ::   "<<pow(-1,l1-lambda)*BMB(N,L,n,l,lambda,n1,l1,n2,l2)<<" = "<<std::endl;
    std::cout<<"<"<<N<<" "<<L<<" "<<n<<" "<<l<<" "<<lambda<<"  |   "<<n2<<"  "<<l2<<"  "<<n1<<"  "<<l1<<"> ::   "<<pow(-1,l1+l)*BMB(N,L,n,l,lambda,n2,l2,n1,l1)<<" = "<<std::endl;
    std::cout<<"<"<<n1<<" "<<l1<<" "<<n2<<" "<<l2<<" "<<lambda<<"  |   "<<n<<"  "<<l<<"  "<<N<<"  "<<L<<"> ::   "<<pow(-1,l2+L)*BMB(n1,l1,n2,l2,lambda,n,l,N,L)<<" = "<<std::endl;
    std::cout<<"============================================================"<<std::endl;
}


void BMB_ortonormality_checking(int n1, int l1, int n2,int l2, int lambda, int n1_q, int l1_q, int n2_q,int l2_q ){
    /// The monshinsky brackets consist in an orthonormal basis of coefficients to expand states in different basis.
    /// Summing up over non-zero brackets n,l,N,L between different n1, n2, l1, l2 lead zero, as a Kronecker delta over them.

    int lenght = 0;
    int rho,N;
    double Total = 0.0;
    double first, second;
    rho = max(2*n1 + 2*n2 + l1 + l2 ,2*n1_q + 2*n2_q + l1_q + l2_q);

    for(int l=0; l<=rho; l++){
        for(int L = abs(lambda-l); L <=lambda+l; L++){
            for(int n = 0; n <= (floor(rho/2)); n++){
                for( N=0; N<=(floor(rho/2)); N++){

                    first = BMB(n,l,N,L,lambda,n1,l1,n2,l2);
                    second = BMB(n,l,N,L,lambda,n1_q,l1_q,n2_q,l2_q);
                    //std::cout<<"Element"<<std::endl;
                    //std::cout<<n1<<" "<<n2<<" "<<l1<<" "<<l2<<" "<<lambda<<" "<<n<<" "<<N<<" "<<l<<" "<<L<<" "<<first<<std::endl;
                    //std::cout<<n1<<" "<<n2_q<<" "<<l1_q<<" "<<l2_q<<" "<<lambda<<" "<<n<<" "<<N<<" "<<l<<" "<<L<<" "<<second<<std::endl;

                    Total += first * second;
                    lenght++;
                }

            }
        }
    }

    std::cout<< "Norm of all n,l,N,L Coefficients for n1, l1,n2,l2 ::"<<Total<<std::endl;
    //std::cout<< "Number of BMBs :"<<lenght<<std::endl;
}


void check_BMB_article_values(){
    std::ifstream f_a;
    f_a.open("data/BMB_article.txt",std::fstream::in);

    int n,l,N,L,lambda,n1,l1,n2,l2;
    double article_value,BMB_numerical;
    double tol = 1e-5;
    int errors = 0;
    std::cout<<"n l N L lambda n1 l1 n2 l2 ::   artc value / my BMB"<<std::endl;
    while(!f_a.eof()){
        f_a >>n>>l>>N>>L>>lambda>>n1>>l1>>n2>>l2>> article_value;
        BMB_numerical = BMB(n,l,N,L,lambda,n1,l1,n2,l2);
        if(abs(BMB_numerical-article_value)<tol){
            std::cout<<"<"<<n<<" "<<l<<" "<<N<<" "<<L<<" "<<lambda<<"      "<<n1<<"  "<<l1<<"  "<<n2<<"  "<<l2<<"> ::   "<<article_value<<" / "<<BMB_numerical<<std::endl;
            BMB_ortonormality_checking(n1,l1,n2,l2, lambda, n1,l1,n2,l2);
            BMB_symmetry_checking(n,l,N,L,lambda,n1,l1,n2,l2);
        }
        else{
            std::cout<<"<"<<n<<" "<<l<<" "<<N<<" "<<L<<" "<<lambda<<"      "<<n1<<"  "<<l1<<"  "<<n2<<"  "<<l2<<"> ::   "<<article_value<<" / "<<BMB_numerical<<" /ERROR" <<std::endl;
            errors++;
        }
    }
    std::cout<<errors<<" Errors"<<std::endl;

}


void check_BMB(int N_input){
    std::cout << "0 to compute BMB from n1=n2=0 to N_input in a non-zero range" << std::endl ;
    std::cout << "1 to prove symmetry properties of one BMB" << std::endl;


    // testing time efficiency in computing values from n1 and n2 = 0 to the N shell
    // assuming only l and lambda from 0(s) to 6(i) and lambdas in range of addition of ls
    int lenght = 0;
    int rho,N;
    double BMB_value;

    std::ofstream f_salida;
    f_salida.open("data/Monshinsky_brackets.txt",std::fstream::out);

    for(int n1=0; n1 <= N_input; n1++){
        for(int n2=0; n2 <= N_input;n2++){
            for (int l1 =0; l1 <=6; l1++){
                for( int l2 =0; l2<=6; l2++){
                    std::cout<< "n1:"<<n1<<" n2:"<<n2<<" l1:"<<l1<<" l2:"<<l2<<std::endl;
                    for(int lambda=abs(l1-l2);lambda <= (l1+l2);lambda++){
                        // The range of coefficients n,l,N,L goes from the minimum to
                        // maximum according to energy-angular conditions
                        rho = 2*n1 + 2*n2 + l1 + l2;

                        for(int l=0; l<=rho; l++){
                            for(int L = abs(lambda-l); L <=lambda+l; L++){
                                for(int n = 0; n <= (floor(rho/2)- L - l); n++){
                                    if(n>0){
                                        if((rho - l - L)<0){
                                           N = ceil((rho - l - L)/2) - n;
                                        }
                                        else{N = floor((rho - l - L)/2) - n;}

                                        if(N>=0){
                                            //std::cout<<n1<<" "<<n2<<" "<<l1<<" "<<l2<<" "<<lambda<<" "<<n<<" "<<N<<" "<<l<<" "<<L<<" i"<<lenght<<std::endl;
                                            BMB_value = BMB(n,l,N,L,lambda,n1,l1,n2,l2);
                                            lenght ++;
                                            //std::cout<<n1<<" "<<n2<<" "<<l1<<" "<<l2<<" "<<lambda<<" "<<n<<" "<<N<<" "<<l<<" "<<L<<" : "<<BMB_value<< "  / i:"<<lenght<<std::endl;
                                            f_salida <<lenght <<" "<< n<<" "<<l<<" "<<N<<" "<<L<<" "<<lambda<<" "<<n1<<" "<<l1<<" "<<n2<<" "<<l2<<" "<<BMB_value<<std::endl;
                                        }
                                    }

                                }
                            }
                        }

                    }
                }
            }
        }
    }
    std::cout<< "Number of BMBs :"<<lenght<<std::endl;
    f_salida.close();
}



//////////////////////////////////////////////////////////////////////////////////
///////////////////////////       MATRIX ELEMENTS      ///////////////////////////
//////////////////////////////////////////////////////////////////////////////////

void check_ME_Art(int N, int n_Max){

    /// This function apply both Matrix elements from Monshinsky Paper and from
    /// MAPLE implementation paper. There are several differences between them,
    /// starting by the arguments to give.
    /// Monshinsky M.E are also apply in the jj Coupling scheme, so I'm  going
    /// to calculate them in the two body jj-coupled SHO Functions.
    /// KET = |n1,l1,j1,(1/2,)n2,l2,j2,(1/2,),J,M,T,M_T>


    QN_2body_jj_Coupling BRA, KET;

    int n1,n2,l1,l2, n1_q,l1_q,n2_q,l2_q;
    Fraction j1,j2,J,M, j1_q,j2_q,J_q,M_q;
    int T,M_T;

    T = 1;
    M_T = 1;
    M = 0;
    M_q = M;

    BRA.T = T;
    BRA.M = M;
    BRA.M_T = M_T;

    KET.T = T;
    KET.M = M_q;
    KET.M_T = M_T;

    double numeric_M_Central, numeric_M_LS, numeric_M_Tensor;
    double numeric_Art_Central, numeric_Art_LS, numeric_Art_Tensor;

    double b_param = b_lenght(30);

    std::ofstream f_salida_1;
    f_salida_1.open("data/ME_Talmi.txt",std::fstream::out);
    std::ofstream f_salida_2;
    f_salida_2.open("data/ME_Talmi_for_Python.txt",std::fstream::out);
    int lenght = 0;
    double tol = 1e-6;

    f_salida_1 <<"<n1 l1  j1 n2  l2  j2 (J M) |V_12| n1'  l1' j1' n2' l2' j2' (J' M')>    ="<<
        "   M_Central      M_LS      M_Tensor  // Art_Central    Art_LS    Art_Tensor "<< std::endl;
    // Radial Quantum number
    for(n1 = 0; n1 <= n_Max; n1++){
        for(n1_q = 0; n1_q <= n1; n1_q++){
            for(n2 = 0; n2 <= n_Max; n2++){
                for(n2_q = 0; n2_q <= n2; n2_q++){
                    std::cout<< "_____________________________________"<<std::endl;
                    std::cout<< "n1,n1',n2,n2' = "<<n1<<n1_q<<n2<< n2_q<< "  / length = "<<lenght<<std::endl;

                    // Major Oscillator Quantum number (N_min = 0) fixed.
                    // In the SHO, l is determined odd or even due to (n,N)
                    //int l1_lim, l1_q_lim, l2_lim,l2_q_lim;
                    //
                    for(l1 = N; l1 >=0; l1 -= 2){
                        for(l2 = N; l2 >=0 ; l2 -= 2){
                            for(l1_q = N; l1_q >=0; l1_q -= 2){
                                for(l2_q = N; l2_q >=0 ; l2_q -= 2){
                                    std::cout<< "l1,l1',l2,l2' = "<<l1<<l1_q<<l2<< l2_q<< "  / length = "<<lenght<<std::endl;

                                    //Possible j with s = 1/2
                                    for( j1 = Fraction(2*l1-1,2); j1 <= Fraction(2*l1+1,2); j1 += 1){
                                        for( j2 = Fraction(2*l2-1,2); j2 <= Fraction(2*l2+1,2); j2 += 1){
                                            for( j1_q = Fraction(2*l1_q-1,2); j1_q <= Fraction(2*l1_q+1,2); j1_q += 1){
                                                for( j2_q = Fraction(2*l2_q-1,2); j2_q <= Fraction(2*l2_q+1,2); j2_q += 1){

                                                    // J possible (as M don't appear in the ME, it remains on 0)
                                                    for( J = abs((j1-j2)); J <= (j1 + j2); J +=1){
                                                        for( J_q = abs((j1_q-j2_q)); J_q <= (j1_q + j2_q); J_q +=1){

                                                            // Fill the Wave Functions
                                                            BRA = {n1,l1,j1, n2,l2,j2, J,M,T,M_T};
                                                            KET = {n1_q,l1_q,j1_q, n2_q,l2_q,j2_q, J_q,M_q,T,M_T};

                                                            //Compute
                                                            numeric_M_LS = M_E_Spin_Orbit_jj(BRA, KET,0, b_param );
                                                            //numeric_M_LS = round(numeric_M_LS * pow(10,order)) / pow(10,order);
                                                            numeric_M_Tensor = M_E_Tensor_jj(BRA, KET,0, b_param);
                                                            //numeric_M_Tensor = round(numeric_M_Tensor * pow(10,order)) / pow(10,order);

                                                            numeric_Art_Central = M_E_Central_jj(BRA, KET,0, b_param);
                                                            //numeric_Art_Central = round(numeric_Art_Central * pow(10,order)) / pow(10,order);
                                                            numeric_Art_LS = M_E_Spin_Orbit(BRA, KET,0, b_param);
                                                            //numeric_Art_LS = round(numeric_Art_LS * pow(10,order)) / pow(10,order);
                                                            numeric_Art_Tensor = M_E_Tensor(BRA, KET,0, b_param);

                                                            //numeric_Art_Tensor = round(numeric_Art_Tensor * pow(10,order)) / pow(10,order);

                                                            numeric_M_Central = numeric_Art_Central;

                                                            // Print

                                                            if((fabs(numeric_M_LS) > tol)||(fabs(numeric_M_Tensor) > tol)||(fabs(numeric_Art_Central) > tol)
                                                               ||(fabs(numeric_Art_LS) > tol)||(fabs(numeric_Art_Tensor) > tol)){
                                                                lenght ++;
                                                                //f_salida_1 << fixed;

                                                                f_salida_1 <<"< "<<n1 <<"  "<< l1 <<" "<< j1 <<"  "<< n2 <<"  "<< l2 <<"  "<< j2 <<" ("<<  J <<" "<< M<<") |V_12|  " <<
                                                                            n1_q <<"   "<< l1_q <<"  "<< j1_q <<"   "<< n2_q <<"   "<< l2_q <<"   "<< j2_q << " ("<<  J_q <<"  "<< M_q<<" )>   =   "<<
                                                                            fixed<< setprecision(7)<<numeric_M_Central << "   " <<fixed<< setprecision(7)<< numeric_M_LS << "   " <<fixed<< setprecision(5)<< numeric_M_Tensor << "  //  " <<fixed<< setprecision(7)<<
                                                                            numeric_Art_Central << "   " <<fixed<< setprecision(7)<< numeric_Art_LS << "   " <<fixed<< setprecision(5)<< numeric_Art_Tensor << std::endl;

                                                                f_salida_2 <<n1 <<" "<< l1 <<" "<< j1 <<" "<< n2 <<" "<< l2 <<" "<< j2 <<" "<<  J <<" "<< M<<" / / "
                                                                           <<n1_q <<" "<< l1_q <<" "<< j1_q <<" "<< n2_q <<" "<< l2_q <<" "<< j2_q << " "<<  J_q <<" "<< M_q
                                                                           <<"   "
                                                                           << lenght<< "   " << numeric_M_Central << "   " << numeric_M_LS << "   " << numeric_M_Tensor
                                                                           << "   " << numeric_Art_Central << "   " << numeric_Art_LS << "   " << numeric_Art_Tensor << std::endl;

                                                            }

                                                        }
                                                    }
                                                    // End J

                                                }
                                            }
                                        }
                                    }
                                    // End j

                                }
                            }
                        }
                    }
                    // End OQN

                }
            }
        }
    }
    // End RQN
    std::cout<< "Number of ME :"<<lenght<<std::endl;
    f_salida_1.close();
    f_salida_2.close();

}



void check_ME_multipolar(QN_1body_jj WF, int order){
    std::cout << " ////////////   MULTIPOLAR CHECKING  ///////////" << std::endl;
    
    std::ofstream f_salida;
    f_salida.open("data/ME_Multipolar.txt",std::fstream::out);

    Fraction J;
    double numeric;
    for(int I = 0 ;I <= order; I++){
        std::cout<< " ===== Maximum Order I = "<< I << " ======" <<std::endl;
        for(J = Fraction(0); J <= Fraction(2*WF.j); J +=1){
            QN_2body_jj_Coupling BRA = {WF.n, WF.l, WF.j,WF.n,WF.l,WF.j, J,0,1,1};
            QN_2body_jj_Coupling KET = {WF.n, WF.l, WF.j,WF.n,WF.l,WF.j, J,0,1,1};
            std::cout << " --- J = " << J << " --- "<<std::endl;
            numeric = M_E_Spin_Isospin_Multipolar(BRA, KET, I , b_lenght(50));
            f_salida << I <<" "<<J<<" "<<numeric<<std::endl;
	    std::cout << "Multipolar M.E = "<< numeric << std::endl;
        }

    }

    f_salida.close();

}


void check_ME_SDI(){

    // Previous Declarations;
    int Table, Number_of_distint_elements;
    QN_1body_jj *SET;

    // Choosing of parameters
    std::cout << " Table to represent: Table 8.";
    std::cin >> Table;
    std::cout<< std::endl;

    // Declare At the star the quantum numbers and notation employed
    std::ofstream f_salida_1;
    f_salida_1.open("data/ME_SDI.txt",std::fstream::out);

    switch(Table){
        case 1 :
            {
                std::cout<< " 1 = 0s_1/2 ; 2 = 0p_3/2 ; 3 = 0p_1/2  "<<std::endl;
                f_salida_1<< " Table 8." << Table <<std::endl;
                f_salida_1 << " 1 = 0s_1/2 ; 2 = 0p_3/2 ; 3 = 0p_1/2  " <<std::endl;
                Number_of_distint_elements = 3;

                SET = new QN_1body_jj [Number_of_distint_elements];
                SET[0] = {0, 0, Fraction(1,2)};
                SET[1] = {0, 1, Fraction(3,2)};
                SET[2] = {0, 1, Fraction(1,2)};
            }
            break;
        case 2 :
            {
                std::cout<< " 1 = 0d_5/2 ; 2 = 1s_1/2 ; 3 = 0d_3/2   "<<std::endl;
                f_salida_1<< " Table 8." << Table <<std::endl;
                f_salida_1 << " 1 = 0d_5/2 ; 2 = 1s_1/2 ; 3 = 0d_3/2   " <<std::endl;
                Number_of_distint_elements = 3;

                SET = new QN_1body_jj [Number_of_distint_elements];
                SET[0] = {0, 2, Fraction(5,2)};
                SET[1] = {1, 0, Fraction(1,2)};
                SET[2] = {0, 2, Fraction(3,2)};
            }
            break;
        case 3 :
            {
                std::cout<< " 1 = 0p_1/2 ; 2 = 0p_3/2 ; 3 = 0d_3/2 ; 4 = 0d_5/2 ; 5 = 1s_1/2 "<<std::endl;
                f_salida_1<< " Table 8." << Table <<std::endl;
                f_salida_1<<" 1 = 0p_1/2 ; 2 = 0p_3/2 ; 3 = 0d_3/2 ; 4 = 0d_5/2 ; 5 = 1s_1/2 "<<std::endl;
                Number_of_distint_elements = 5;

                SET = new QN_1body_jj [Number_of_distint_elements];
                SET[0] = {0, 1, Fraction(1,2)};
                SET[1] = {0, 1, Fraction(3,2)};
                SET[2] = {0, 2, Fraction(3,2)};
                SET[3] = {0, 2, Fraction(5,2)};
                SET[4] = {1, 0, Fraction(1,2)};
            }
            break;
        case 4 :
            {
                std::cout<< " 1 = 0d_3/2 ; 2 = 0d_5/2 ; 3 = 1s_1/2 ; 4 = 0f_7/2"<<std::endl;
                f_salida_1<< " Table 8." << Table <<std::endl;
                f_salida_1<< " 1 = 0d_3/2 ; 2 = 0d_5/2 ; 3 = 1s_1/2 ; 4 = 0f_7/2"<<std::endl;
                Number_of_distint_elements = 4;

                SET = new QN_1body_jj [Number_of_distint_elements];
                SET[0] = {0, 2, Fraction(3,2)};
                SET[1] = {0, 2, Fraction(5,2)};
                SET[2] = {1, 0, Fraction(1,2)};
                SET[3] = {0, 3, Fraction(7,2)};
            }
            break;
        default:
            // Set your own collection of numbers, no limits.
            {
                std::cout << "Set number of different Quantum Numbers : ";
                std::cin >> Number_of_distint_elements;
                std::cout<< std::endl;
                SET = new QN_1body_jj [Number_of_distint_elements];

                int n_input, l_input, j_num;

                f_salida_1 <<" Choosen Shells"<<std::endl;
                f_salida_1 <<" Item  (n,l, j ) "<<std::endl;
                for (int i = 0; i < Number_of_distint_elements; i++){

                    std::cout << " Set n radial number (0,1,...): " ;
                    std::cin >> n_input ;
                    std::cout<< std::endl;
                    std::cout << " Set l Orbital Number (0 = s, 1 = p, ...): " ;
                    std::cin >> l_input;
                    std::cout<< std::endl;
                    std::cout << " Set the numerator (odd) of the total angular mometum j = ?/2: " ;
                    std::cin >> j_num ;
                    std::cout<< std::endl;
                    std::cout << " j = "<< Fraction (j_num,2) <<std::endl;

                    SET[i] = {n_input, l_input, Fraction (j_num,2)};
                    f_salida_1<<"  "<< i+1<<"    ("<<n_input<<","<<l_input<<","<<Fraction(j_num,2)<<")"<<std::endl;
                }
                std::cout << "Passing of Quantum Numbers Completed" << std::endl;
            }
    }

    QN_2body_jj_Coupling BRA, KET;
    QN_1body_jj a,b,c,d;
    double numerical_value_1, numerical_value_2, numerical_value_3;
    double tol = 1e-30;

    f_salida_1 <<std::endl;
    f_salida_1<< "a b c d  (J T)  =  <ab|V(SDI)|cd>   <ab|V(SDI)|cd>(C!=0)  <ab|V(MSDI)|cd>" << std::endl;
    f_salida_1<< "-------------------------------------------------------------------------------" << std::endl;
    // Running of permutations from the set selected.
    for( int i_a = 0; i_a < Number_of_distint_elements; i_a++){
        a = {SET[i_a].n , SET[i_a].l, SET[i_a].j};

        for( int i_b = 0; i_b < Number_of_distint_elements; i_b++){
            b = {SET[i_b].n , SET[i_b].l, SET[i_b].j};

            for( int i_c = 0; i_c < Number_of_distint_elements; i_c++){
                c = {SET[i_c].n , SET[i_c].l, SET[i_c].j};

                for( int i_d = 0; i_d < Number_of_distint_elements; i_d++){
                    d = {SET[i_d].n , SET[i_d].l, SET[i_d].j};
                    f_salida_1 <<std::endl;

                    // Display of all J,J' coupled momentum for the BRA and KET
                    // and evaluate the matrix element only for coincidence.
                    for( Fraction J = abs((a.j - b.j)); J <= (a.j + b.j); J+= 1){
                        for( Fraction J_q = abs((c.j - d.j)); J_q <= (c.j + d.j); J_q+= 1){
                            if(J == J_q){

                                 // Try T = 0,1 and export the values on .txt
                                 for( int T = 0; T <=1; T++){
                                    BRA = {a.n,a.l, a.j,  b.n,b.l, b.j, J,0,T,0};
                                    KET = {c.n,c.l, c.j,  d.n,d.l, d.j, J_q,0,T,0};

                                    numerical_value_1 = M_E_SDI(BRA,  KET,  false);
                                    numerical_value_2 = M_E_SDI(BRA,  KET,  true);
                                    numerical_value_3 = M_E_MSDI(BRA,  KET);
                                    if((fabs(numerical_value_1) > tol) && (fabs(numerical_value_2) > tol) &&(fabs(numerical_value_3) > tol)){
                                        f_salida_1 << i_a+1 <<" "<< i_b+1 <<" "<< i_c+1 <<" "<< i_d+1 <<"  ("
                                                   << J << " " <<T <<")  =  " <<fixed<< setprecision(7)<< numerical_value_1
                                                   << "             "<< fixed<< setprecision(7)<< numerical_value_2
                                                   << "         " << fixed<< setprecision(7) << numerical_value_3 << std::endl;


                                    }

                                 }

                            }

                        }
                    }


                }
            }
        }

    }

    f_salida_1.close();

}


///////////////////////////////////////////////////////////////////////////////////
///////////////		PRINTING QUANTUM NUBERS FROM STRUCTS 	   ////////////////
///////////////////////////////////////////////////////////////////////////////////

void print_QN_2body_jj ( QN_2body_jj_Coupling WF){
    
    QN_1body_jj State_1 = {WF.n1,WF.l1,WF.j1};
    QN_1body_jj State_2 = {WF.n2,WF.l2,WF.j2};
    
    std::cout << "n1,n2 = " << WF.n1 <<" "<< WF.n2 << std::endl;
    std::cout << "l1,l2 = " << WF.l1 <<" "<< WF.l2 << std::endl;
    std::cout << "j1,j2,J,T = " << WF.j1 <<" "<< WF.j2 <<" "<< WF.J << " "<< WF.T << std::endl;
    std::cout << "Antoine(1):  " << index_antoine(State_1) << " //  Antoine(2):  " << index_antoine(State_2) << std::endl;
    
    std::cout << std::endl;

}
void print_QN_1body_jj ( QN_1body_jj WF){

    std::cout << "n,l,j,m = " << WF.n <<" "<< WF.l << " "<< WF.j << " "<< WF.m << std::endl;
    std::cout << "Antoine:  " << index_antoine(WF) << std::endl;
    std::cout << std::endl;

}
void print_QN_2body_rad( QN_2body_radial WF){

    std::cout << "n1,l1,n2,l2 = " << WF.n1 <<" "<< WF.l1 << " "<< WF.n2 << " "<< WF.l2 << std::endl;
    std::cout << "lambda, mu = " << WF.lambda <<" "<< WF.mu << std::endl;
    std::cout << std::endl;

}
void print_QN_1body_rad( QN_1body_radial WF){

    std::cout << "n,l,m_l = " << WF.n <<" "<< WF.l << " "<< WF.m_l << std::endl;
    std::cout << std::endl;

}



/////////////////////////////////////////////////////////////////////////////////////////
// Loop in order to check Chasman's Conversion
/*
    bool Correct;
    int error = 0;
    for(int n=0; n<=10; n++){
      for(int l=0; l<=12; l++){
	for(int m=(-l); m<=l; m++){

	  QN_1body_radial bra = {n,l,m};
	  Correct = check_Chasman_Conversion_Spherical2Cartesian(bra);

	  if(!Correct){
	    error++;
	    std::cout<< "ERROR "<<error <<std::endl;
	  }
	  std::cout<<std::endl;
	}
      }
    }
    std::cout<<error<<"  errors"<<std::endl;
*/

bool check_Chasman_Conversion_Spherical2Cartesian (QN_1body_radial Spherical){

    std::cout << " The Conversion 'Spherical to Cartesian' Checking  begins ... " << std::endl;
    int top_a = 2*Spherical.n + Spherical.l;
    int n_a_z;
    complex< double > aux_coeff_a;

    complex< double > Norm (0.,0.);
    print_QN_1body_rad (Spherical);
    for(int n_a_x = 0; n_a_x <= top_a; n_a_x++ ){
        for(int n_a_y = 0; n_a_y <= (top_a - n_a_x); n_a_y++ ){
            //for(int n_a_z = n_a_y; n_a_z <= (top_a - n_a_x - n_a_y); n_a_z++ ){
		n_a_z = top_a - n_a_y - n_a_x;
		std::cout << "(N= "<< top_a <<") (n_x,n_y,n_z) = "<< n_a_x<<" , "<<n_a_y<<" , "<<n_a_z;//<< std::endl;

		aux_coeff_a = conj(conversion_NLM_XYZ( Spherical , n_a_x, n_a_y, n_a_z));

		std::cout << " //  Coef= " << aux_coeff_a  <<std::endl;
		//std::cout << std::endl;

		Norm += aux_coeff_a * conj(aux_coeff_a);
	}
    }

    bool continua;
    if(fabs(Norm.imag()) > 1.e-10){std::cin>>continua;}

    std::cout <<" **** Norm= " << Norm <<std::endl;
    if(fabs(abs(Norm)-1) > 1.e-3 ){return false;}

    return true;
}

void check_SHO_WF( QN_1body_jj & WF, double R_MAX, double b_param ){
    std::cout << " ////////////   SHO WAVE FUNCTION CHECKING  ///////////" << std::endl;
    
    std::ofstream f_salida;
    f_salida.open("data/SHO.txt",std::fstream::out);

    int N_MAX = 1000;
    double r, psi, dr;
    dr = R_MAX / N_MAX;
    
    for(int i=0; i<N_MAX; i++){
	r = i * dr;
        psi = SHO_Radial(WF.n, WF.l, r, b_param);
        f_salida << r <<" "<<psi<<" "<<pow(psi,2)<<std::endl;

    }

    f_salida.close();
    std::cout << " DONE " << std::endl;
}

void check_Cuadrature_Convergence(QN_2body_radial Q_N_radial_BRA,QN_2body_radial Q_N_radial_KET,double lambda,
				      double mu_param,double b_param,int Max_order, int precision){
    std::cout << "----------------------------------------------------------------------- " << std::endl;
    std::cout << "    Display of the Gauss quadrature of a Rosenfeld multipolar integral  " << std::endl;
    std::cout << " in a (R_Max - R_Min) range for a certain Quadrature order, and  the    " << std::endl;
    std::cout << " same integral of the same function but with the quadratures applied to " << std::endl;
    std::cout << " 20 subdivisions of the same range (same order).		          " << std::endl;
    std::cout << "----------------------------------------------------------------------- " << std::endl;
  
    double gauss;
    double upgrade;
    double difference,difference_convergnece, difference_quad;
    double upgrade_prev = 1.;
    double gauss_prev   = 1.;
    bool cont = true;
    int first = 0;
    for(int i = 2; i<= Max_order; i++){
	
	run_Gauss_Legendre(i);
	gauss = Radial_4_order_Gauss_Quadrature(Q_N_radial_BRA,Q_N_radial_KET,i, lambda, mu_param, b_param);
	
	if(cont){
	    upgrade = Radial_4_order_Subdivided_Quadrature(Q_N_radial_BRA,Q_N_radial_KET,i,lambda, mu_param, b_param) ;
	}
	
	difference_convergnece = fabs(upgrade_prev - upgrade);
	difference = fabs(gauss - upgrade);
	difference_quad = fabs(gauss - gauss_prev);
	
	if(difference_convergnece < pow(0.1, precision)){
	      // There is not necessary more subdivided iterations. 
	     if(first == 0){
		  std::cout <<" ** Subdivided quadratures Precision ACHIVED ** " << std::endl;
		  first +=1;
	     }
	     
	     cont = false;
	}
	if(difference_quad < pow(0.1, precision)){
	     std::cout <<" ** Gauss Quadrature Precision ACHIVED ** " << std::endl;
	     break;
	}
	upgrade_prev = upgrade;
	gauss_prev = gauss;
	
	std::cout <<" Quad. Ord.= "<< i
		  <<"  / * / R4 Gauss= "<< std::fixed << std::setprecision(precision) << gauss
		  <<"  / * / R4 Subdivided= "<< std::fixed << std::setprecision(precision) << upgrade
		  <<"  / * / abs(difference) = "     << difference << std::endl;
    }
    
}

double check_Coupling_BBME_Methods(QN_2body_jj_Coupling Q_Numbers_left, QN_2body_jj_Coupling Q_Numbers_right,
					  double mu_param,double b_param, int Max_order, int precision, bool Moshinsky_Method){
      // Max_order << Order of the Gauss Legendre quadrature
  
    std::cout << "----------------------------------------------------------------------- " << std::endl;
    std::cout << "    Complete execution of the coupling process in order to obtain the   " << std::endl;
    std::cout << " differences between the final state in LS coupled/uncoupled form.      " << std::endl;
    std::cout << "    Second step is to see the diferences of the result with coupling    " << std::endl;
    std::cout << " and uncoupled scheme over the angular LST operators.			  " << std::endl;
    std::cout << "  		Non Antisimetriced Matrix Elements                        " << std::endl;
    std::cout << "----------------------------------------------------------------------- " << std::endl;      
    
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
  
    // Zero Cases
    if(M != M_q){
	std::cout <<"M != M_q"<<std::endl;
	return 0;}
    if(J != J_q){
	std::cout <<"J != J_q"<<std::endl;
	return 0;}
    if(T != T_q){
	std::cout <<"T != T_q"<<std::endl;
	return 0;}
    if(((n_a == n_b) && (l_a == l_b) && (j_a == j_b)) || ((n_c == n_d) && (l_c == l_d) && (j_c == j_d))){
	if((int(J + Fraction(T)))%2 != 1){
	    std::cout <<"(J + T)%2 != 1"<<std::endl;
	    return 0;}
    }
    std::cout << "Not a Zero Case " << std::endl;
    
    double mu_BB[2] = {mu_param , mu_param };
       
    ///////////////////////////////////////////////////////
    /////////    Non_antisimetriced_BB_M_E    /////////////
    ///////////////////////////////////////////////////////
    
    double Isospin_Radial[2] = {-1 , 0} ;
    double Isospin_Spin[2]   = {0 , 0} ;
    double radial_value[2]   = { 0. , 0. };
    
    Fraction m_b, m_d;
    int m_l_a, m_l_b, m_l_c, m_l_d;
    int S;
    double unc_J_total, unc_j, coup_S, spin_value;   	//"unc" from uncoupling
    double Result_coup_averages = 0.;
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
				//std::cout << "Uncoupling j LS= " << unc_j << "  /*/ Uncoupling J total= " << unc_J_total << "  /=/ " << unc_j * unc_J_total << std::endl;

                                if( 1.e-10 < fabs(unc_j)){
                                        // Evaluation of the Spin part (before because it's simpler than the radial one)
                                    spin_value = 0.;

				    for(S = 0; S <= 1; S++){
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
							    * (2*S*(S+1) - 3);; // <S> = (2*S*(S+1) - 3);
					}
				    }
				    
				    //Result_coup_averages += (unc_J_total * unc_j);
                                    //std::cout << std::endl;
                                    //std::cout << "Spin value= " << spin_value << std::endl;
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
					    
					    LAMBDA = (b_param*sqrt(2.)/ mu_BB[0]);
					    radial_value[0] = Radial_BB_Moshinsky(radial_a, radial_b, radial_c, radial_d, mu_BB[0],sqrt(2)*b_param);
					    LAMBDA = (b_param*sqrt(2.)/ mu_BB[1]);
					    radial_value[1] = Radial_BB_Moshinsky(radial_a, radial_b, radial_c, radial_d, mu_BB[1],sqrt(2)*b_param);

					//}
				    }
                                    //* /

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
                                    //std::cout << "Radial value [0][1]= " << radial_value[0] << " , " << radial_value[1] << std::endl;

                                        // Add everything
                                    Result_coup_averages += (unc_J_total * unc_j) *
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

    
    // EXCHANGE OPERATORS /////////////////////////////////
    ///  Non_antisimetriced_BBME_by_exchange_LSjj      ////
    ///////////////////////////////////////////////////////
    double A_BB[2] = { -1 , 0 };     // Central
    double B_BB[2] = { 0 , 0 };     // Spin
    double C_BB[2] = { 0 , 0 };     // Isospin
    double D_BB[2] = { 0 , 0 };     // Spin * Isospin
    //double C_BB[2] = { -3.541 , -3.541};     // Isospin
    //double D_BB[2] = { -8.25053 , -8.25053 };     // Spin * Isospin
    // Inverted expression
    double V_W_BB[2] = { A_BB[0] - B_BB[0] - C_BB[0] + D_BB[0] ,
			  A_BB[1] - B_BB[1] - C_BB[1] + D_BB[1] };    // Wigner, Central
    double V_M_BB[2] = { -4*D_BB[0] , -4*D_BB[1] };                     // Majorana, Spatial Permutation
    double V_B_BB[2] = { (2*B_BB[0]) - (2*D_BB[0]) ,
			  (2*B_BB[1]) - (2*D_BB[1]) };                // Barlett, Spin Permutation
    double V_H_BB[2] = { (2*C_BB[0]) - (2*D_BB[0]) ,
			(2*C_BB[1]) - (2*D_BB[1]) };                // Heisenberg, Isospin Permutation
			
			
    //double radial_value[2]   = { 0. , 0. };
    radial_value[0]   = 0.;
    radial_value[1]   = 0.;
    double exchange_value[2] = { 0. , 0. };
    

    unc_J_total = 0;
    double LSjj_1,LSjj_2, unc_l;   	//"unc" from uncoupling
    double Result_coup_LSjj = 0.;
    //QN_1body_radial radial_a, radial_b, radial_c, radial_d;
    
    
    //for( int L = abs(l_a-l_b); L <= l_a + l_b; L++ ){
    for( int L = max(abs(l_a-l_b), abs(l_c-l_d)); L <= min(l_a + l_b, l_c + l_d); L++ ){
        for( int S = 0; S <= 1; S++){

            LSjj_1 = LS_jj_coupling_Coeff( Fraction(l_a),Fraction(1,2),j_a,  Fraction(l_b), Fraction(1,2),j_b,
						Fraction(L), Fraction(S),J);
					   
	    std::cout << "L,S  LSjj1= "<<L<< " ,"<< S<<" ,"<< LSjj_1<<std::endl; 
            if(1.e-10 < fabs(LSjj_1) ){
		LSjj_2 = LS_jj_coupling_Coeff( Fraction(l_c),Fraction(1,2),j_c,  Fraction(l_d), Fraction(1,2),j_d,
						Fraction(L), Fraction(S),J);

		//std::cout << "L,S  LSjj2= "<<L<< " ,"<< S<<" ,"<< LSjj_2<<std::endl; 
                unc_J_total = LSjj_1 * LSjj_2;

                std::cout << "LSjj1= " << LSjj_1 << " , LSjj2= " << LSjj_2 <<std::endl;

                if( 1.e-10 < fabs(unc_J_total)){

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
		      
  
		      exchange_value[0] = (V_W_BB[0] + (pow(-1,L)*V_M_BB[0]) + (pow(-1,S)*V_B_BB[0]) + (pow(-1,T)*V_H_BB[0]));
		      exchange_value[1] = (V_W_BB[1] + (pow(-1,L)*V_M_BB[1]) + (pow(-1,S)*V_B_BB[1]) + (pow(-1,T)*V_H_BB[1]));
		      
		      std::cout << "Radial[0]= "<<  radial_value[0]<<"  , Radial[1]= "<< radial_value[1] <<std::endl;
		      std::cout << "exch val[0]= "<<  exchange_value[0]<<"  , exch val[1]= "<< exchange_value[1] <<std::endl;
		      
		      // Add everything
		      //Result_coup_LSjj += (unc_J_total);
		      Result_coup_LSjj += (unc_J_total) * ( (radial_value[0] * exchange_value[0]) + (radial_value[1] * exchange_value[1]));
		      

                }


            }// end J1

        }
    }    
    
    
    
    ///////////////////////////////////////////////////////
    //////  Non_antisimetriced_BBME_by_exchange      //////
    ///////////////////////////////////////////////////////
    
    //double radial_value[2]   = { 0. , 0. };
    //double exchange_value[2] = { 0. , 0. };
    radial_value[0]   = 0.;
    radial_value[1]   = 0.;
    exchange_value[0]   = 0.;
    exchange_value[1]   = 0.;

    //Fraction m_b, m_d;
    //int m_l_a, m_l_b, m_l_c, m_l_d;
    //double unc_J_total, unc_j, 
    unc_J_total = 0.;
    unc_j = 0.;
    coup_S = 0;
    double coup_L;   	//"unc" from uncoupling
    double Result_coup_exch = 0.;
    //QN_1body_radial radial_a, radial_b, radial_c, radial_d;

    for( Fraction m_a = Fraction(-1*j_a.numerator, j_a.denominator); m_a <= j_a ; m_a += 1){
        m_b = M - m_a;
        for( Fraction m_c = Fraction(-1*j_c.numerator, j_c.denominator); m_c <= j_c ; m_c += 1){
            m_d = M - m_c;
	    /*std::cout << "========== Uncoup J =========="<<std::endl;
	    std::cout <<"j_a,j_b,J, m_a, m_b, M "<< j_a<<" "<<j_b<<" "<<J<<" "<< m_a<<" "<< m_b<<" "<< M<<" CG1= "<< Clebsh_Gordan(j_a,j_b,J, m_a, m_b, M) << std::endl;
	    std::cout <<"j_c,j_d,J, m_c, m_d, M "<< j_c<<" "<<j_d<<" "<<J<<" "<< m_c<<" "<< m_d<<" "<< M<<" CG2= "<< Clebsh_Gordan(j_c,j_d,J, m_c, m_d, M) << std::endl;
	    */
            unc_J_total = Clebsh_Gordan(j_a,j_b,J, m_a, m_b, M) * Clebsh_Gordan(j_c,j_d,J, m_c, m_d, M);
	    //std::cout << "Uncoupling J total= " << unc_J_total << std::endl;
	    
            if( 1.e-10 < fabs(unc_J_total)){
		//std::cout << " --------Uncoup jj ----------"<<std::endl;
                for( Fraction m_s_a = Fraction(1,2); m_s_a >= Fraction(-1,2) ; m_s_a -= 1){
                    m_l_a = m_a - m_s_a;
                    for( Fraction m_s_b = Fraction(1,2); m_s_b >= Fraction(-1,2) ; m_s_b -= 1){
                        m_l_b = m_b - m_s_b;
                        for( Fraction m_s_c = Fraction(1,2); m_s_c >= Fraction(-1,2) ; m_s_c -= 1){
                            m_l_c = m_c - m_s_c;
                            for( Fraction m_s_d = Fraction(1,2); m_s_d >= Fraction(-1,2) ; m_s_d -= 1){
                                m_l_d = m_d - m_s_d;
				/*
				std::cout <<"(l_a),(1,2),j_a, (m_l_a),m_s_a, m_a: ["<< Fraction(l_a)<<" "<<Fraction(1,2)<<" "<<j_a
				<<" , "<< Fraction(m_l_a)<<" "<< m_s_a<<" "<< m_a<<"] CG_a= "<<  Clebsh_Gordan(Fraction(l_a),Fraction(1,2),j_a, Fraction(m_l_a),m_s_a, m_a ) << std::endl;
				std::cout <<"(l_b),(1,2),j_b, (m_l_b),m_s_b, m_b: ["<< Fraction(l_b)<<" "<<Fraction(1,2)<<" "<<j_b
				<<" , "<< Fraction(m_l_b)<<" "<< m_s_b<<" "<< m_b<<"] CG_b= "<< Clebsh_Gordan(Fraction(l_b),Fraction(1,2),j_b, Fraction(m_l_b),m_s_b, m_b )  << std::endl;
				std::cout <<"(l_c),(1,2),j_c, (m_l_c),m_s_c, m_c: ["<< Fraction(l_c)<<" "<<Fraction(1,2)<<" "<<j_c
				<<" , "<< Fraction(m_l_c)<<" "<< m_s_c<<" "<< m_c <<"] CG_c= "<< Clebsh_Gordan(Fraction(l_c),Fraction(1,2),j_c, Fraction(m_l_c),m_s_c, m_c )  << std::endl;
				std::cout <<"(l_d),(1,2),j_d, (m_l_d),m_s_d, m_d: ["<< Fraction(l_d)<<" "<<Fraction(1,2)<<" "<<j_d
				<<" , "<< Fraction(m_l_d)<<" "<< m_s_d<<" "<< m_d<<"] CG_d= "<< Clebsh_Gordan(Fraction(l_d),Fraction(1,2),j_d, Fraction(m_l_d),m_s_d, m_d )  << std::endl;
				*/
                                unc_j = Clebsh_Gordan(Fraction(l_a),Fraction(1,2),j_a, Fraction(m_l_a),m_s_a, m_a )
                                        * Clebsh_Gordan(Fraction(l_b),Fraction(1,2),j_b, Fraction(m_l_b),m_s_b, m_b )
                                          * Clebsh_Gordan(Fraction(l_c),Fraction(1,2),j_c, Fraction(m_l_c),m_s_c, m_c )
                                            * Clebsh_Gordan(Fraction(l_d),Fraction(1,2),j_d, Fraction(m_l_d),m_s_d, m_d );
				//std::cout << std::endl;
				//std::cout << "Unc j= " << unc_j << "  // Unc J total= " << unc_J_total << std::endl;

                                if( 1.e-10 < fabs(unc_j)){
				    
				    //Result_coup_exch += (unc_J_total * unc_j);
				    //std::cout << "++ Result =" << Result_coup_exch << std::endl;
				    
  
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
					    
					    LAMBDA = (b_param*sqrt(2.)/ mu_BB[0]);
					    radial_value[0] = Radial_BB_Moshinsky(radial_a, radial_b, radial_c, radial_d, mu_BB[0], b_param*sqrt(2.));
					    LAMBDA = (b_param*sqrt(2.)/ mu_BB[1]);
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
					    //std::cout << std::endl;
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
                                    //Result_coup_exch += (unc_J_total * unc_j) * ( (radial_value[0] * exchange_value[0]) + (radial_value[1] * exchange_value[1]));
				    Result_coup_exch += (unc_J_total * unc_j) * ( (radial_value[0] * exchange_value[0]) + (radial_value[1] * exchange_value[1]));
                                    //*/

                                } // end if j

                            } // end for m_s
                        }
                    }
                }
            }// end if J

        }
    }
    
    
    std::cout << " averages=" << Result_coup_averages <<"  , exch= " <<  Result_coup_exch <<"  , LSjj= " << Result_coup_LSjj << std::endl; 
    
    ///////////////////////////////////////////////////////
    //////  Solution by Multipoles and Quadratures   //////
    //////       Non_antisimetriced_BB_Multi         //////
    ///////////////////////////////////////////////////////
    int lambda_max = max(max(max((l_a + l_c),int(j_a + j_c)),int(j_a + j_b)), int(j_c+j_d)) + 4;
    
    // Convergnece test of the radial integral:
    double gauss;
    
    QN_2body_radial Q_N_radial_BRA = {n_a,l_a,n_b,l_b, 0, 0};
    QN_2body_radial Q_N_radial_KET = {n_c,l_c,n_d,l_d, 0, 0};
    /*
    for (int GL_order = 2; GL_order <= Max_order; GL_order++){
	run_Gauss_Legendre(GL_order);
	gauss = Radial_4_order_Gauss_Quadrature(Q_N_radial_BRA,Q_N_radial_KET, GL_order, lambda_max, mu_param, b_param);
	
	std::cout <<" Quad. Ord.= "<< GL_order
		  <<"  / * / R4 Gauss= "<< std::fixed << std::setprecision(precision) << gauss <<std::endl;
    }
     */
    
    // Evaluation of the matrix element
    run_Gauss_Legendre(Max_order);
    double Numeric = 0; //Non_antisimetriced_BB_Multi( Q_Numbers_left,  Q_Numbers_right, lambda_max,  b_param);
    
    double aux_U;
    double rad_aux;
    
    //int GL_order = 50;
    //run_Gauss_Legendre(GL_order);
    
    // Reduced struct argument for the integrals.
    //QN_2body_radial Q_N_radial_BRA;
    //Q_N_radial_BRA = {BRA.n1, BRA.l1, BRA.n2, BRA.l2, 0, 0};// Total angular momentum and 3rd component are irrelevant
    //QN_2body_radial Q_N_radial_KET;
    //Q_N_radial_KET = {KET.n1, KET.l1, KET.n2, KET.l2, 0, 0}; 
    double *Parameters;
    Parameters = new double[4];
    
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
	    aux_U = U_Coeff(Q_Numbers_left, Q_Numbers_right ,lambda, Parameters);
	    //std::cout <<" // U dir = "<<aux_U<<std::endl;
	    
	    if(fabs(aux_U) > 1.0e-10){

		//rad_aux = Radial_4_order(Q_N_radial_BRA, Q_N_radial_KET, lambda,b_param);
		rad_aux = Radial_4_order_Gauss_Quadrature(Q_N_radial_BRA, Q_N_radial_KET, Max_order ,lambda, mu_BB[i],b_param);
		
		Numeric +=  rad_aux * aux_U;
		//std::cout <<"Radial 4= "<< rad_aux << std::endl;
	    }

	    //std::cout <<"Direct: R= "<< Radial_4_order(Q_N_radial_BRA, Q_N_radial_KET, lambda) << " // U = "<<U_Coeff(BRA, KET ,lambda)<<std::endl;
	}
    }

    
    
    Result_coup_exch *= sqrt(2*J+1);
    Result_coup_averages *= sqrt(2*J+1);
    Result_coup_LSjj *= sqrt(2*J+1);
    Numeric *= sqrt(2*J+1);
    
    std::cout << " averages= " << Result_coup_averages <<"  , exch= " <<  Result_coup_exch <<"  , LSjj= " << Result_coup_LSjj << std::endl;
    std::cout << " Multi (correct) = " << Numeric << std::endl;
    
}
