#include "../include/factorials.h"
#include "../include/Fractions.h"
#include "../include/Angular_Functions.h"
#include "../include/Index_Coefficients.h"
#include <cmath>
#include <iostream>

using namespace std;

double BM_Bracket_00(int n, int l, int N, int L, int lambda, int l1, int l2){
    // To initialize the recurrence relations for the monshinsky brackets
    // we first have to compute the value of <n l, N L,lambda | 0 l1 0 l2  lambda>


    if ((n<0)||(l<0)||(N<0)||(L<0)){return 0;}
    // Following conditions should not be necessary to be evaluated at this point,
    // because <BM_Bracket_00> is called from <BMB>, which take care about it.
    // Uncomment if it's going to call <BM_Bracket_00> directly.
    //if ((lambda<0)||(l1<0)||(l2<0)){return 0;}

    double sum = 0.;

    // sum restricted to l1+l and l2+L possible values

    for(int x= max(abs(l-l1),abs(L-l2)); x<=min(l + l1,L +l2); x++){
        sum += (2*x +1) * A_coefficient(l1,l,l2,L,x) * Racah_Coefficient(l,L,l1,l2,lambda,x);
    }

    int *Num;
    Num = new int[4];
    int *Den;
    Den = new int[6];

    Num[0] = l1;
    Num[1] = l2;
    Num[2] = n + l;
    Num[3] = N + L;

    Den[0] = 2*l1;
    Den[1] = 2*l2;
    Den[2] = n;
    Den[3] = 2*n + 2*l + 1;
    Den[4] = N;
    Den[5] = 2*N + 2*L + 1;

    double Coef_00 = exp(0.5 * (log((2*l + 1)*(2*L + 1)) + factorial_function(Num,Den,4,6) - ((l + L)*log(2))));
    delete[] Num;
    delete[] Den;

    return  Coef_00 * pow(-1,n + l + L - lambda) * sum;
}


double BMB(int n, int l, int N, int L, int lambda, int n1, int l1, int n2, int l2){

    // Tis is a recurrence function where we decrease index n1 to 0, then n2 to 0 and compute
    // the six index non-zero value for BM_Bracket_00.

    // to reduce the amount of recursive calls, I stablish the Energy and Angular Conditions at
    //the beginning of the function

    // Condition Already evaluated within the function <ME_ri2_times_BMB>, called every time that
    // the BMB is n1,n2!=0:
    //if ((n<0)||(l<0)||(N<0)||(L<0)){
    //    //cout<<"un indice negativo"<<endl;
    //    return 0;
    // }

    // Because these tree numbers remains unchanged over the computation
    if ((lambda<0)||(l1<0)||(l2<0)){return 0;}

    else if((n1<0)||(n2<0)){return 0;} // non sense arguments

    else if((2*n + l + 2*N + L)==(2*n1 + l1 + 2*n2 + l2)){
        if((abs(l1 - l2)<= lambda) && (lambda <=(l1+l2))){
            if((abs(l - L)<= lambda) && (lambda <=(L+l))){
                // First, reduce the n1 to 0, then, reduce n2 to 0.
                if((n1==0)){
                    //cout<<"n1=0 /"<<n1<<endl;
                    if(n2==0){
                        // n1 = n2 = 0 -> Calculate the BM00
                        //cout<<"n2=0 /"<<n2<<"/ calculo BM 00"<<endl;
                        return BM_Bracket_00(n, l, N, L, lambda, l1, l2);
                    }
                    else{
                        // Decrease the index n2 to 0,
                        //cout<<"Decrease n2 ::n2-1="<<n2-1<<endl;
                        double suma = 0.;
                        suma += ME_ri2_times_BMB(2 ,lambda, n, l, N, L,    n, l, N-1, L,    n1,  l1,  n2-1, l2);
                        suma += ME_ri2_times_BMB(2 ,lambda, n, l, N, L,    n, l-1, N-1, L+1,    n1,  l1,  n2-1, l2);
                        suma += ME_ri2_times_BMB(2 ,lambda, n, l, N, L,    n, l-1, N, L-1,    n1,  l1,  n2-1, l2);
                        suma += ME_ri2_times_BMB(2 ,lambda, n, l, N, L,    n-1, l, N, L,    n1,  l1,  n2-1, l2);
                        suma += ME_ri2_times_BMB(2 ,lambda, n, l, N, L,    n-1, l+1, N-1, L+1,    n1,  l1,  n2-1, l2);
                        suma += ME_ri2_times_BMB(2 ,lambda, n, l, N, L,    n-1, l+1, N, L-1,    n1,  l1,  n2-1, l2);

                        return sqrt(1/(n2*(n2 + l2 + 0.5))) * suma;
                    }

                }
                else{
                    // decrease the index n1 to 0, keeping the n2 value unchanged
                    //cout<<"Decrease n1 ::n1-1="<<n1-1<<endl;
                    double suma = 0.;
                    suma += ME_ri2_times_BMB(1 ,lambda, n, l, N, L,    n, l, N-1, L,    n1-1,  l1,  n2, l2);
                    suma += ME_ri2_times_BMB(1 ,lambda, n, l, N, L,    n, l-1, N-1, L+1,    n1-1,  l1,  n2, l2);
                    suma += ME_ri2_times_BMB(1 ,lambda, n, l, N, L,    n, l-1, N, L-1,    n1-1,  l1,  n2, l2);
                    suma += ME_ri2_times_BMB(1 ,lambda, n, l, N, L,    n-1, l, N, L,    n1-1,  l1,  n2, l2);
                    suma += ME_ri2_times_BMB(1 ,lambda, n, l, N, L,    n-1, l+1, N-1, L+1,    n1-1,  l1,  n2, l2);
                    suma += ME_ri2_times_BMB(1 ,lambda, n, l, N, L,    n-1, l+1, N, L-1,    n1-1,  l1,  n2, l2);

                    return sqrt(1/(n1*(n1 + l1 + 0.5))) * suma;
                }
            }
            else{return 0;}
        }
        else{return 0;}
    }
    else{return 0;}

}
