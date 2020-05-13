#include "../include/Fractions.h"

#include<math.h>
#define _USE_MATH_DEFINES
#include<iostream>

int Factorial_Maximum_Dimension = 1000;
double *Factorial_int_list = new double[Factorial_Maximum_Dimension + 1];
double *Factorial_hfi_list = new double[2*Factorial_Maximum_Dimension + 1];
double *Double_Factorial_list = new double[Factorial_Maximum_Dimension + 1];

void run_Factotial_Base(){
    /// Execute and fill a factorial pre-calculated base of factorials.
    /// First run in main()

    int N_Max = Factorial_Maximum_Dimension;

    /// 0 and 1 elements
    Factorial_int_list[0] = 0;
    Factorial_int_list[1] = 0;

    Factorial_hfi_list[0] = 0;
    Factorial_hfi_list[1] = log(sqrt(M_PI)/2);

    Double_Factorial_list[0] = 0;
    Double_Factorial_list[1] = 0;

    for(int i = 1; i < N_Max;i++){

        // the integer factorials
        Factorial_int_list[i] = Factorial_int_list[i-1] + log(i);

        // the half-integer factorials
        Factorial_hfi_list[2*i + 1] = Factorial_hfi_list[2*i - 1] + log(2*i+1) - log(2);

        // the half-integer factorials which are also integers
        Factorial_hfi_list[2*(i+1)] = Factorial_hfi_list[2*i] + log(i+1);

        // double integer factorials
        if( i%2 == 0){
            Double_Factorial_list[i] = Double_Factorial_list[i-2] + log(i);
        }
        else{
            Double_Factorial_list[i] = Double_Factorial_list[i-2] + log(i);
        }
        // double half-integer factorials (non-negative-integer values)
    }
    Factorial_int_list[N_Max] = Factorial_int_list[N_Max-1] + log(N_Max);
    std::cout<<"Factorial Basis finished (Dim:"<< Factorial_Maximum_Dimension <<")"<<std::endl;
}

double factorial(Fraction n){
    return Factorial_hfi_list[int(2*n)];     // logarithm value returned
}

double factorial(int n){
    return Factorial_int_list[ n ];
}

double factorial_function(int * numerator, int * denominator,int N_num,int N_den){
    double suma = 0.;
    // Add the numerator
    for(int i=0; i<N_num; i++){
        if(numerator[i]<0){
            std::cout<<" Negative Factorial"<<std::endl;
        }
        else{
            if(numerator[i] > Factorial_Maximum_Dimension){
                std::cout<<"Int. num. Factorial(" <<denominator[i] << ") exceeds the base ("
			  << Factorial_Maximum_Dimension<<"), increment its dimension "<<std::endl;
                return -1;
            }
            else{
                suma += Factorial_int_list[numerator[i]];
            }

        }
    }
    // Subtract the denominator
    for(int i=0; i<N_den; i++){
        if(denominator[i]<0){
            std::cout<<" Negative Factorial"<<std::endl;
        }
        else{
            if(denominator[i] > Factorial_Maximum_Dimension){
                std::cout<<"Int. den. Factorial(" <<denominator[i] << ") exceeds the base ("
			  << Factorial_Maximum_Dimension<<"), increment its dimension "<<std::endl;
                return -1;
            }
            else{
                suma -= Factorial_int_list[denominator[i]];
            }

        }
    }
    return suma; // logarithm value returned
}

// Overloading of factorial_function

double factorial_function(Fraction * num, Fraction * den,int N_num,int N_den){
    double suma = 0.;
    // Add the numerator
    for(int i= 0; i<N_num; i++){
        if(num[i]<0){
            std::cout<<" Negative Factorial"<<std::endl;
            //bool continua;
            //std::cin >> continua;
            return -1;
        }

        else if( num[i].denominator == 1){
            if(num[i].numerator > Factorial_Maximum_Dimension){
                std::cout<<" Factorial exceeds the base, increment its dimension "<<std::endl;
                return -1;
            }
            else{
                suma += Factorial_hfi_list[2*num[i].numerator];
            }
        }
        else{  //if( num[i].denominator == 2)
            if(num[i].numerator > Factorial_Maximum_Dimension){
                std::cout<<" Factorial exceeds the base, increment its dimension "<<std::endl;
                return -1;
            }
            else{
                suma += Factorial_hfi_list[num[i].numerator];
            }
        }
    }

    // Subtract the denominator
    for(int i= 0; i<N_den; i++){
        if(den[i]<0){
            std::cout<<" Negative Factorial"<<std::endl;
            //bool continua;
            //std::cin >> continua;
            return -1;
        }
        else if( den[i].denominator == 1){
            if(den[i].numerator > Factorial_Maximum_Dimension){
                std::cout<<" Factorial exceeds the base, increment its dimension "<<std::endl;
                return -1;
            }
            else{
                suma -= Factorial_hfi_list[2*den[i].numerator];
            }
        }
        else{  //( den[i].denominator == 2)
            if(den[i].numerator > 2*Factorial_Maximum_Dimension){
                std::cout<<" Factorial exceeds the base, increment its dimension "<<std::endl;
                return -1;
            }
            else{
                suma -= Factorial_hfi_list[den[i].numerator];
            }
        }
    }
    return suma; // logarithm value returned
}



double double_factorial(int n){
    double Aux = 0.;
    int i_0;

    if(n % 2 == 0){i_0 = 2;}
    else{i_0 = 1;}

    for( int i = i_0; i <= n+1; i += 2){
        Aux += log(i);
    }
    return Aux;// logarithm value returned
}


double double_factorial_function(int * numerator, int * denominator,int N_num,int N_den){

    double suma = 0.;
    // Add the numerator
    for(int i=0; i<N_num; i++){
        if(numerator[i]<0){
            std::cout<<" Negative Factorial"<<std::endl;
            bool continua;
            std::cin >> continua;
        }
        else{
            if(numerator[i] > Factorial_Maximum_Dimension){
                std::cout<<" Factorial exceeds the base, increment its dimension "<<std::endl;
                return -1;
            }
            else{
                suma += Double_Factorial_list[numerator[i]];
            }

        }
    }
    // Subtract the denominator
    for(int i=0; i<N_den; i++){
        if(denominator[i]<0){
            std::cout<<" Negative Factorial"<<std::endl;
            //bool continua;
            //std::cin >> continua;
        }
        else{
            if(denominator[i] > Factorial_Maximum_Dimension){
                std::cout<<" Factorial exceeds the base, increment its dimension "<<std::endl;
                return -1;
            }
            else{
                suma -= Double_Factorial_list[denominator[i]];
            }

        }
    }
    return suma; // logarithm value returned

}


double gamma_function(Fraction x, int &sign){

    // factorial related results for gamma function.
    if(x.denominator == 1){
        if(x.numerator <= 0){
          std::cout<< "Gamma (0 o x<0) = inf" << std::endl;
          sign = 1;
          return -1e30;
        }

        return factorial(x-1);  // logarithmic output
    }

    // Negative arguments of Gamma function have a phase, not interchangeable by
    // a logarithm value, so there is a dummy variable that changes if the function
    // give a negative sign. So it is only necessary to declare it if the Gamma
    // could be applied over negative arguments.
    //
    // sign variable have to be placed after the call to the function.
    // >> exp(gamma_function(Fraction(-1,2),sign))* sign
    //
    // This value is employed to prevent future problems about argument form.
    else{
        // negative half integer factorials follow the same formula than positive ones.
        // but there is not evaluated for
        if( x.numerator < 0){

            //std::cout << " Gamma(x < 0), "<<std::endl;
            // Change of sign
            sign = pow(-1, (1 - x.numerator)/2);
            //std::cout<<"Sign"<<sign<<std::endl;
            // Logarithmic output
            return log( M_PI / exp(factorial(Fraction( - x.numerator, x.denominator))));
        }
        else{
            if(x.numerator == 1){
                sign = 1;
                return log(sqrt(M_PI));  // logarithmic output
            }
            else{
                sign = 1;
                return factorial(x - 1); // logarithmic output
            }
        }

    }

}



