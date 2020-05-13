#include "../include/Fractions.h"
#include "../include/Angular_Functions.h"
#include "../include/BM_Brackets.h"
#include "../include/Nuclear_Matrix_Elements_M.h"
#include "../include/Index_Coefficients.h"
#include <cmath>

using namespace std;
/// Matrix elements from "Symbolic algorithms for the
///computation of Moshinsky brackets and nuclear matrix elements"

//int n1, int l1,Fraction j1, Fraction s1, int n2, int l2, Fraction j2, Fraction s2, Fraction J,Fraction M,
//                 int n1_q, int l1_q, Fraction j1_q, Fraction s1_q, int n2_q, int l2_q, Fraction j2_q, Fraction s2_q, Fraction J_q,Fraction M_q

double M_E_Central_jj(const QN_2body_jj_Coupling & Q_Numbers_left,const QN_2body_jj_Coupling & Q_Numbers_right, int option, double b_param){
    // This generalize the expression to include jj-coupling scheme. Devanathan (6.73)
    // Is better to pass two arrays with the 9 variables.

    // Passing of the parameters
    int n1 = Q_Numbers_left.n1;
    int l1 = Q_Numbers_left.l1;
    Fraction j1 = Q_Numbers_left.j1;
        //Fraction s1 = Q_Numbers_left.s1;
    int n2 = Q_Numbers_left.n2;
    int l2 = Q_Numbers_left.l2;
    Fraction j2 = Q_Numbers_left.j2;
        //Fraction s2 = Q_Numbers_left.s2;
    Fraction J = Q_Numbers_left.J;
    Fraction M = Q_Numbers_left.M;
    int T = Q_Numbers_left.T;

    int n1_q = Q_Numbers_right.n1;
    int l1_q = Q_Numbers_right.l1;
    Fraction j1_q = Q_Numbers_right.j1;
        //Fraction s1_q = Q_Numbers_right.s1;
    int n2_q = Q_Numbers_right.n2;
    int l2_q = Q_Numbers_right.l2;
    Fraction j2_q = Q_Numbers_right.j2;
        //Fraction s2_q = Q_Numbers_right.s2;
    Fraction J_q = Q_Numbers_right.J;
    Fraction M_q = Q_Numbers_right.M;
    int T_q = Q_Numbers_right.T;

    // Starting of computation
    if(M != M_q){return 0;}
    if(J != J_q){return 0;}
    if(T != T_q){return 0;}

    // Isospin and total coupled angular momentum must be odd:
    if((J + T)%2 != 1){return 0;}
    // Sum over S,lambda LS-jj Coupling Coefficients (reduced to 9j symbols)
    double Sum = 0.;

    //int S_min, S_q_min, S_max, S_q_max;
    //S_min = 0;   //abs((s1-s2));
    //S_q_min = 0; //abs((s1_q-s2_q));
    //S_max = 1;   //s1 + s2;
    //S_q_max = 1; //s1_q + s2_q;

    Fraction S,S_q;
    double tol, N1,N2,MEC;

    int lambda,lambda_q;
    for (S = Fraction(0); S<=1; S +=1){
        for(S_q = Fraction(0); S_q <= 1; S_q+=1){
            if(S_q != S){} // delta(S,S')
            else{
                for(lambda = abs(l1 - l2); lambda <= l1 + l2; lambda ++){
                    for(lambda_q = abs(l1_q - l2_q); lambda_q <= (l1_q + l2_q); lambda_q++){
                        if(lambda_q != lambda){} // delta(lambda, lambda')
                        else{
                            /*
                            cout << " Suma lam="<< lambda<< " / S=" << S << endl;
                            cout << " Raices="<<sqrt((((2*j1+ 1)*(2*j2 + 1)*(2*S+1)*(2*lambda+1)) *
                                         ((2*j1_q+ 1)*(2*j2_q + 1)*(2*S_q+1)*(2*lambda_q+1)))) << endl;
                            cout << " 9j=" << Nine_j_Symbol(Fraction(l1),Fraction(1,2),j1,Fraction(l2),Fraction(1,2),j2,Fraction(lambda),S,J) << endl;
                            cout << " 9j=" << Nine_j_Symbol(Fraction(l1_q),Fraction(1,2),j1_q,Fraction(l2_q),Fraction(1,2),j2_q,Fraction(lambda_q),S_q,J_q) << endl;
                            cout << " Central coef =" << M_E_Central(n1, l1, n2, l2, lambda, n1_q, l1_q, n2_q, l2_q, lambda_q) << endl;
                            //cout<<n1<<" "<< l1<<" "<<n2<<" "<<l2<<" "<<lambda<<"/ "<<n1_q<<" "<<l1_q<<" "<< n2_q<<" "<< l2_q <<" "<< lambda_q<<endl;
                            */

                            tol = 1e-9;
                            N1 =  LS_jj_coupling_Coeff(Fraction(l1),Fraction(1,2),j1,Fraction(l2),Fraction(1,2),j2,Fraction(lambda),S,J);
                            if(fabs(N1) > tol){

                                N2 = LS_jj_coupling_Coeff(Fraction(l1_q),Fraction(1,2),j1_q,Fraction(l2_q),Fraction(1,2),j2_q,Fraction(lambda_q),S_q,J_q);
                                if(fabs(N2) > tol){

                                    MEC = M_E_Central(n1, l1, n2, l2, lambda, n1_q, l1_q, n2_q, l2_q, lambda_q,option, b_param);
                                    if (fabs(MEC) > tol){
                                        Sum += sqrt((((2*j1+ 1)*(2*j2 + 1)*(2*S+1)*(2*lambda+1)) *
                                             ((2*j1_q+ 1)*(2*j2_q + 1)*(2*S_q+1)*(2*lambda_q+1)))) *
                                             N1*N2 *MEC;

                                    }
                                }
                            }
                            //Sum += sqrt((((2*j1+ 1)*(2*j2 + 1)*(2*S+1)*(2*lambda+1)) *
                            //             ((2*j1_q+ 1)*(2*j2_q + 1)*(2*S_q+1)*(2*lambda_q+1)))) *
                            //             Nine_j_Symbol(Fraction(l1),Fraction(1,2),j1,Fraction(l2),Fraction(1,2),j2,Fraction(lambda),S,J) *
                            //             Nine_j_Symbol(Fraction(l1_q),Fraction(1,2),j1_q,Fraction(l2_q),Fraction(1,2),j2_q,Fraction(lambda_q),S_q,J_q) *
                            //             M_E_Central(n1, l1, n2, l2, lambda, n1_q, l1_q, n2_q, l2_q, lambda_q);

                        }

                    }
                }

            }
        }
    }


    // Re use of the M_E_Central in the {n,l,lambda} basis.
    return Sum;
}

double Coeff_CLS(int n1,int l1,int  n2, int l2, int lambda, int n1_q, int l1_q, int n2_q,int l2_q,int lambda_q, int option, double b_param){

    int aux = l1_q + l2_q - l1 - l2;
    if((aux)%2 != 0){
        return 0; // considering Monshisky papper
    }
    else{
        int rho = (2*n1 + l1 + 2*n2 + l2);
        int l,n_q;
        double sum = 0.;
        for(int N=0; N<=floor(rho/2); N++){
            for(int n=floor(rho/2); n >= 0; n--){
                for(int L=0; L <= rho; L++){
                    l = (rho - L - 2*(N + n));

                    if(l < 0){}
                    else if(triangular_condition(l,L,lambda)){
                        n_q = n1_q + n2_q - n1 - n2 +(aux/2) + n;

                        double tol = 1e-20;
                        double Racah = Racah_Coefficient(lambda,lambda_q,l,l,1,L);
                        //cout<< "RAcah("<<l<<l<<lambda<<lambda_q<<";"<<2<<L<<") =" << Racah  <<endl;
                        if(fabs(Racah) > tol){
                            double B1 = BMB(n,l,N,L,lambda,n1,l1,n2,l2);
                            //cout<< "BMB="<< B1 <<endl;
                            if(fabs(B1) > tol){
                                double B2 = BMB(n_q ,l,N,L,lambda_q,n1_q,l1_q,n2_q,l2_q);
                                //cout<< "BMB="<<  B2 <<endl;
                                if(fabs(B2) > tol){
                                    double RME = R_M_E_Potential(n,l,n_q,l,option, b_param);
                                    //cout<< "RME=" << RME <<endl;

                                    if(fabs(RME) > tol){
                                        sum +=  B1 * B2 * RME * Racah* sqrt(l*(l+1)*(2*l+1)); // [lambda][lambda_q] are in the other function
                                        //cout<<" sum="<<sum<<endl;

                                        //bool continue_bool;
                                        //cin >> continue_bool;
                                    }
                                }
                            }
                        }


                        //sum +=  BMB(n,l,N,L,lambda,n1,l1,n2,l2)* BMB(n_q ,l,N,L,lambda_q,n1_q,l1_q,n2_q,l2_q)*
                        //        R_M_E_Potential(n,l,n_q,l) * Racah_Coefficient(l,l,lambda,lambda_q,1,L)*
                        //        sqrt(l*(l+1)*(2*l+1)*(2*lambda+1)*(2*lambda_q+1));
                    }
                }
            }
        }
        return sum;
    }

}

double M_E_Spin_Orbit(const QN_2body_jj_Coupling & Q_Numbers_left,const QN_2body_jj_Coupling & Q_Numbers_right,int option, double b_param){

        // Passing of the parameters
    int n1 = Q_Numbers_left.n1;
    int l1 = Q_Numbers_left.l1;
    Fraction j1 = Q_Numbers_left.j1;
    int n2 = Q_Numbers_left.n2;
    int l2 = Q_Numbers_left.l2;
    Fraction j2 = Q_Numbers_left.j2;
    Fraction J = Q_Numbers_left.J;
    Fraction M = Q_Numbers_left.M;
    int T = Q_Numbers_left.T;

    int n1_q = Q_Numbers_right.n1;
    int l1_q = Q_Numbers_right.l1;
    Fraction j1_q = Q_Numbers_right.j1;
    int n2_q = Q_Numbers_right.n2;
    int l2_q = Q_Numbers_right.l2;
    Fraction j2_q = Q_Numbers_right.j2;
    Fraction J_q = Q_Numbers_right.J;
    Fraction M_q = Q_Numbers_right.M;
    int T_q = Q_Numbers_right.T;

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
            if((S_q != S)||(S != 1)){} // delta(S,S')
            else{
                for(lambda = abs(l1 - l2); lambda <= l1 + l2; lambda ++){
                    //for(lambda_q = abs(l1_q - l2_q); lambda_q <= (l1_q + l2_q); lambda_q++){ // Global case
                    for(lambda_q = lambda - 1; lambda_q <= lambda+1; lambda_q ++){              // Restriction, A.17
                        if(lambda_q >= 0){

                            // Encapsulate
                            double tol = 1e-20;
                            double N1 = Nine_j_Symbol(Fraction(l1),Fraction(1,2),j1,Fraction(l2),Fraction(1,2),j2,Fraction(lambda),S,J);
                            //cout<<"09j 1=" << N1 <<endl;
                            if (fabs(N1) > tol){
                                double N2 = Nine_j_Symbol(Fraction(l1_q),Fraction(1,2),j1_q,Fraction(l2_q),Fraction(1,2),j2_q,Fraction(lambda_q),S_q,J_q);
                                //cout<<"09j 2=" << N2 <<endl;
                                if(fabs(N2) > tol){
                                    double Racah = Racah_Coefficient(Fraction(lambda),Fraction(lambda_q),Fraction(1),Fraction(1),Fraction(1),J);
                                    //cout<<"Racah =" << Racah <<endl;
                                    if(fabs(Racah) > tol){
                                       double Coef = Coeff_CLS(n1, l1, n2, l2, lambda, n1_q, l1_q, n2_q, l2_q, lambda_q,option, b_param);
                                       //cout<<"0Coef T " << Coef <<endl;

                                       if(fabs(Coef) > tol){
                                            int rho = 2*(n1 + n2) + (l1 + l2);

                                            // Just if everything is non-zero Sum add up.
                                            Sum += sqrt((((2*j1+ 1)*(2*j2 + 1)*(2*S+1))*
                                                         ((2*j1_q+ 1)*(2*j2_q + 1)*(2*S_q+1)))) *
                                                         (2*lambda+1)*(2*lambda_q+1)*
                                                         N1 *N2 *pow(-1,int((J + rho))) * 2.4494897 *Racah*Coef;
                                            //cout<<"Suma externa ="<<Sum<<endl;

                                            //bool continue_bool;
                                            //cin >> continue_bool;
                                       }
                                    }

                                }
                            }


                            //Sum += sqrt((((2*j1+ 1)*(2*j2 + 1)*(2*S+1)*(2*lambda+1)) *
                            //             ((2*j1_q+ 1)*(2*j2_q + 1)*(2*S_q+1)*(2*lambda_q+1)))) *
                            //             Nine_j_Symbol(Fraction(l1),Fraction(1,2),j1,Fraction(l2),Fraction(1,2),j2,Fraction(lambda),S,J) *
                            //             Nine_j_Symbol(Fraction(l1_q),Fraction(1,2),j1_q,Fraction(l2_q),Fraction(1,2),j2_q,Fraction(lambda_q),S_q,J_q) *
                            //             pow(-1,int((J + rho))) * 2.4494897 *
                            //             Racah_Coefficient(Fraction(lambda),Fraction(lambda_q),Fraction(1),Fraction(1),Fraction(1),J)*
                            //             Coeff_CLS(n1, l1, n2, l2, lambda, n1_q, l1_q, n2_q, l2_q, lambda_q);
                            // the number 2.4494897 = sqrt(6)
                        }
                    }
                }

            }
        }
    }

    return Sum * H_BAR;

}

double Coeff_CT(int n1,int l1,int  n2, int l2, int lambda, int n1_q, int l1_q, int n2_q,int l2_q,int lambda_q ,int option, double b_param){
    int aux = l1_q + l2_q - l1 - l2;
    if(aux%2 != 0){
        return 0;
    }
    else{

        int rho = (2*n1 + l1 + 2*n2 + l2);
        int l,n_q;
        double sum = 0.;

        for(int N=0; N<=floor(rho/2); N++){
            for(int n=floor(rho/2); n >= 0; n--){

                for(int L=0; L <= rho; L++){
                    l = (rho - L - 2*(N + n));
                    //cout<<"N,n,L,l="<<N<<"  "<<n<<"  "<<L<<"  "<<l <<endl;
                    //cin >> continue_bool;
                    //cout << "triangular_condition(l,L,lambda)="<< l<<L<<lambda<< "="<<(triangular_condition(l,L,lambda))<< endl;
                    if(l < 0){}
                    else if(triangular_condition(l,L,lambda)){


                        for(int l_q=l-2; l_q <=l+2; l_q += 2){
                            if(l_q < 0){}
                            else{
                                n_q = n + (n1_q + n2_q - n1- n2)+((aux)/2);

                                // Encapsulated evaluation increase by 5 the evaluations of conditionals for a non zero element,
                                // but reduces the computation excluding the cases in which one of the elements are 0
                                double tol = 1e-20;

                                double Racah = Racah_Coefficient(lambda,lambda_q,l,l_q,2,L);
                                //cout<< "RAcah("<<lambda<<lambda_q<<l<<l_q<<";"<<2<<L<<") =" << Racah  <<endl;
                                if(fabs(Racah) > tol){
                                    double B0 = BMB(n,l,N,L,lambda,n1,l1,n2,l2);
                                    //cout<< "BMB="<<  B0 <<endl;
                                    if(fabs(B0) > tol){
                                        double B1 = BMB(n_q, l_q,N,L,lambda_q,n1_q,l1_q,n2_q,l2_q);
                                        //cout<< "BMB="<< B1 <<endl;
                                        if(fabs(B1) > tol){
                                            double RME = R_M_E_Potential(n_q,l_q,n,l,option, b_param);
                                            //cout<< "RME=" << RME <<endl;
                                            if(fabs(RME) > tol){
                                                double CG = Clebsh_Gordan(l,2,l_q,0,0,0);
                                                //cout<< "CG =" << CG  <<endl;
                                                if(fabs(CG) > tol){

                                                    sum += sqrt(2*l+1) * B0 * B1 * RME * Racah * CG;   // [lambda][lambda_q] are in the other function
                                                    //cout<<" sum="<<sum<<endl;
                                                    //cin >> continue_bool;
                                                }
                                            }
                                        }
                                    }
                                }
                                //cin >> continue_bool;
                                //cout<< "raices=" << sqrt((2*l+1)*(2*lambda+1)*(2*lambda_q+1)) <<endl;

                                //evaluar en orden con el fin de que todo sea no nulo.
                                //sum += BMB(n,l,N,L,lambda,n1,l1,n2,l2) *
                                //    BMB(n_q, l_q,N,L,lambda_q,n1_q,l1_q,n2_q,l2_q) *
                                //    R_M_E_Potential(n_q,l_q,n,l) *
                                //   Racah_Coefficient(lambda,lambda_q,l,l_q,2,L)* Clebsh_Gordan(l,2,l_q,0,0,0) *
                                //    sqrt((2*l+1)*(2*lambda+1)*(2*lambda_q+1)); //reemplazar en la descomposicion para evitar dos raices
                                //cout<<" sum="<<sum<<endl;
                                //cin >> continue_bool;
                            }

                        }
                    }

                }
            }
        }

        return sum;
    }

}

double M_E_Tensor(const QN_2body_jj_Coupling & Q_Numbers_left,const  QN_2body_jj_Coupling & Q_Numbers_right ,int option, double b_param){

        // Passing of the parameters
    int n1 = Q_Numbers_left.n1;
    int l1 = Q_Numbers_left.l1;
    Fraction j1 = Q_Numbers_left.j1;
    int n2 = Q_Numbers_left.n2;
    int l2 = Q_Numbers_left.l2;
    Fraction j2 = Q_Numbers_left.j2;
    Fraction J = Q_Numbers_left.J;
    Fraction M = Q_Numbers_left.M;
    int T = Q_Numbers_left.T;


    int n1_q = Q_Numbers_right.n1;
    int l1_q = Q_Numbers_right.l1;
    Fraction j1_q = Q_Numbers_right.j1;
    int n2_q = Q_Numbers_right.n2;
    int l2_q = Q_Numbers_right.l2;
    Fraction j2_q = Q_Numbers_right.j2;
    Fraction J_q = Q_Numbers_right.J;
    Fraction M_q = Q_Numbers_right.M;
    int T_q = Q_Numbers_right.T;

    // Starting of computation
    if(M != M_q){return 0;}
    if(J != J_q){return 0;}
    if(T != T_q){return 0;}

    // Isospin and total coupled angular momentum must be odd:
    if((J + T)%2 != 1){return 0;}

    // Sum over S,lambda LS-jj Coupling Coefficients (reduced to 9j symbols)
    double Sum = 0.;

    Fraction S,S_q;
    int rho = 2*(n1 + n2) + (l1 + l2);
    int lambda,lambda_q;
    for (S = Fraction(0); S<=1; S +=1){
        for(S_q = Fraction(0); S_q <= 1; S_q+=1){
            if((S_q != S)||(S != 1)){} // delta(S,S')
            else{

                for(lambda = abs(l1 - l2); lambda <= l1 + l2; lambda ++){
                    //for(lambda_q = abs(l1_q - l2_q); lambda_q <= (l1_q + l2_q); lambda_q++){ // Global case
                    for(lambda_q = lambda - 2; lambda_q <= lambda+2; lambda_q ++){              // Restriction, A.17
                        if(lambda_q < 0){}
                        else{
                            //cout<< "S, S_q, Lambda, Lambda_q="<<S <<" " <<S_q<<" "<< lambda<< " "<<lambda_q<<endl;

                            // Encapsulating the parts to be calculated.
                            double tol = 1e-20;

                            double N1 = LS_jj_coupling_Coeff(Fraction(l1),Fraction(1,2),j1,Fraction(l2),Fraction(1,2),j2,Fraction(lambda),S,J);
                            //cout<<"09j 1=" << N1 <<endl;
                            if(fabs(N1) > tol){

                                double N2 = LS_jj_coupling_Coeff(Fraction(l1_q),Fraction(1,2),j1_q,Fraction(l2_q),Fraction(1,2),j2_q,Fraction(lambda_q),S_q,J_q);
                                //cout<<"09j 2=" << N2 <<endl;
                                if(fabs(N2) > tol){

                                    double Racah = Racah_Coefficient(Fraction(lambda),Fraction(lambda_q),Fraction(1),Fraction(1),Fraction(2),J);
                                    //cout<<"0W() =" << Racah <<endl;
                                    if(fabs(Racah) > tol){

                                        // Coeff_CT involves a maximum of 5 evaluations, its cheaper to evaluate it for the last.
                                        double Coef = Coeff_CT(n1, l1, n2, l2, lambda, n1_q, l1_q, n2_q, l2_q, lambda_q,option, b_param);
                                        //cout<<"0Coef T " << Coef <<endl;
                                        if(fabs(Coef) > tol){

                                            // Just if everything is non-zero Sum add up.
                                            Sum += sqrt((((2*j1+ 1)*(2*j2 + 1)*(2*S+1)) *
                                                         ((2*j1_q+ 1)*(2*j2_q + 1)*(2*S_q+1)))) *
                                                         (2*lambda+1)*(2*lambda_q+1) *
                                                         pow(-1,int ((rho - J + 1))) * 10.954451 *
                                                         N1 * N2 * Racah * Coef;
                                            //cout<<"Suma externa ="<<Sum<<endl;
                                            //bool continue_bool;
                                            //cin >> continue_bool;
                                        }
                                    }
                                }
                            }

                            //cout<<"0potencia " <<pow(-1,int ((rho - J + 1))) * 10.954451  <<endl;


                            //Sum += sqrt((((2*j1+ 1)*(2*j2 + 1)*(2*S+1)*(2*lambda+1)) *
                            //             ((2*j1_q+ 1)*(2*j2_q + 1)*(2*S_q+1)*(2*lambda_q+1)))) *
                            //             Nine_j_Symbol(Fraction(l1),Fraction(1,2),j1,Fraction(l2),Fraction(1,2),j2,Fraction(lambda),S,J) *
                            //             Nine_j_Symbol(Fraction(l1_q),Fraction(1,2),j1_q,Fraction(l2_q),Fraction(1,2),j2_q,Fraction(lambda_q),S_q,J_q) *
                            //             pow(-1,int ((rho - J + 1))) * 10.954451 *
                            //             Racah_Coefficient(Fraction(lambda),Fraction(lambda_q),Fraction(1),Fraction(1),Fraction(2),J)*
                            //             Coeff_CT(n1, l1, n2, l2, lambda, n1_q, l1_q, n2_q, l2_q, lambda_q);
                            //                                                              //the number 10.954451 = sqrt(120)
                            //cout<<"Suma externa ="<<Sum<<endl;
                            //cin >> continue_bool;

                        }

                    }
                }

            }
        }
    }

    return Sum * pow(H_BAR,2);

}
