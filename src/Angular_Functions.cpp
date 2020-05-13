#include "../include/Fractions.h"
#include "../include/factorials.h"
//#include <math.h>
#include<iostream>
#include <cmath>
//#pragma once
using namespace std;

inline int absolute(int a, int b){
    if (a>b)
        return a-b;
    else if (a<=b)
        return b-a;
}

bool triangular_condition(int a, int b, int c){
    if (((abs(a-b) <= c) && (c <= (a+b)))||((abs(c-b) <= a) && (a <= (c+b)))||((abs(a-c) <= b) && (b <= (a+c)))){
        return true;
    }
    else{return false;}
}
bool triangular_condition(Fraction a, Fraction b, Fraction c){
    if (((abs(a-b) <= c) && (c <= (a+b)))||((abs(c-b) <= a) && (a <= (c+b)))||((abs(a-c) <= b) && (b <= (a+c)))){
        return true;
    }
    else{return false;}
}

double Delta(Fraction a, Fraction b, Fraction c){
    // auxiliar function of Racah_Coefficients, no necessary integer type

    // triangular condition is not applicable in the way of the following conditionals,
    // so, there is needed to check that " half integrability" is kept.
    if ((a.denominator + b.denominator + c.denominator)%2==0){return 2;}

    if (((abs(a-b) <= c) && (c <= (a+b)))||((abs(c-b) <= a) && (a <= (c+b)))||((abs(a-c) <= b) && (b <= (a+c)))){ //||((abs(c-b) <= a) && (a <= (c+b)))||((abs(a-c) <= b) && (b <= (a+c)))
        // this condition ensure that the factorial arguments are non negative
        Fraction * numerator;
        numerator = new Fraction [3];
        numerator[0] = (a+b-c);
        numerator[1] = (a-b+c);
        numerator[2] = (b+c-a);

        if((numerator[0]<0) || (numerator[1]<0) || (numerator[2]<0)){
            std::cout<<"numerador negativo:"<<numerator[0]<<","<<numerator[1]<<","<<numerator[2]<<std::endl;
        }
        // evaluate the half-integer factorial
        Fraction * denominator;
        denominator = new Fraction[1];
        denominator[0] = (a+b+c+1);

        double aux = factorial_function(numerator, denominator,3,1);

        delete[] numerator;
        delete[] denominator;

        return aux;
    }
    else{
        // triangular condition must be satisfied between coeficients. Else Delta is 0

        // It's demonstrable that, if coefficients satisfy triangular condition, then,
        // the above function of factorials is less than 1, and then, aux is lower than 0.
        // return value 1 indicate to the following function to return 0 value to the Racah coeficient
        return 2; // there is aux = 0 (for example, for W(0,2,0,2,2,0) or W(1,1,0,0,0,1))
    }
}

double Delta(int a, int b, int c){
    // auxiliar function of Racah_Coefficients

    if (((abs(a-b) <= c) && (c <= (a+b)))||((abs(c-b) <= a) && (a <= (c+b)))||((abs(a-c) <= b) && (b <= (a+c)))){
      //||((abs(c-b) <= a) && (a <= (c+b)))||((abs(a-c) <= b) && (b <= (a+c)))
        // this condition ensure that the factorial arguments are non negative
        int * numerator;
        numerator = new int [3];
        numerator[0] = (a+b-c);
        numerator[1] = (a-b+c);
        numerator[2] = (b+c-a);

        if((numerator[0]<0) || (numerator[1]<0) || (numerator[2]<0)){
            std::cout<<"numerador negativo:"<<numerator[0]<<","<<numerator[1]<<","<<numerator[2]<<std::endl;
        }

        int * denominator;
        denominator = new int[1];
        denominator[0] = (a+b+c+1);

        double aux = factorial_function(numerator, denominator,3,1);

        delete[] numerator;
        delete[] denominator;

        return aux;
    }
    else{
        // triangular condition must be satisfied between coeficients. Else Delta is 0

        // It's demonstrable that, if coefficients satisfy triangular condition, then,
        // the above function of factorials is less than 1, and then, aux is lower than 0.
        // return value 1 indicate to the following function to return 0 value to the Racah coeficient
        return 2; // there is aux = 0 (for example, for W(0,2,0,2,2,0) or W(1,1,0,0,0,1))
    }
}

double Racah_Coefficient(Fraction a,Fraction b, Fraction c, Fraction d, Fraction e, Fraction f){
    // For the moment, this function runs over every int-half integer values,
    // as LS or jj coupling angular momentum does. Only is needed to change
    // the call to the function Delta.

    // Condition for the strictly non negative values
    if((a<0) || (b<0) || (c<0) || (d<0)|| (e<0)|| (f<0)){return 0;}
  
    //evaluation of the deltas before the sum and the product speed the code up.
    double D1 = Delta(a,b,e);

    //Deltas are lower or equal to 0 when triangular condition is fulfilled
    if (D1<1){
        double D2 = Delta(d,c,e);//Delta(c,d,e); from devanathan
        if(D2<1){
            double D3 = Delta(a,c,f);
            if(D3<1){
                double D4 = Delta(d,b,f);//Delta(b,d,f); from devanathan
                if(D4<1){

                    int c1, c2, c3, c4, c5, c6, c7;
                    // c quantities and then x are integer values, because at this conditional
                    // point, a,b,c ... have to fulfill triangular conditions, which means that
                    // e and f are the result of the sum of the others (in a certain order)
                    // ensuring that the integer and half-integer character of a,b,c,d are fixed
                    // to them and easily leads that c_i are always integer type

                    // x higher or equal to (and -1)
                    c1 = int(a + b + e);
                    c2 = int(c + d + e);
                    c3 = int(a + c + f);
                    c4 = int(b + d + f);
                    // x lower or equal to
                    c5 = int(a + b + c + d);
                    c6 = int(a + d + e + f);
                    c7 = int(b + c + e + f);

                    // sumatory goes from the lower of the "lower than" to the highest of "higher than"
                    int major = min(min(c5,c6),c7);
                    // if x is lower than the lowest value of the "lower than" coefficients, then it is for the rest
                    int minor = max(max(max(max(c1,c2),c3),c4),-1) ;
                    // if x is bigger than the biggest "higher than" coefficient, then it is for the rest

                    //std::cout << "x max=" << mayor << "  x min" << minor << std::endl;

                    int *numerator;
                    numerator = new int[1];
                    int *denominator;
                    denominator = new int[7];
                    double aux;
                    double sum = 0.;

                    for(int x= min(minor, major); x<=max(minor,major);x++){
                        numerator[0] =  x + 1;
                        denominator[0] = x - c1;
                        denominator[1] = x - c2;
                        denominator[2] = x - c3;
                        denominator[3] = x - c4;
                        denominator[4] = c5 - x;
                        denominator[5] = c6 - x;
                        denominator[6] = c7 - x;

                        aux = factorial_function(numerator, denominator,1,7);
                        sum += pow(-1,x) * exp(aux);
                    }
                    delete[] numerator;
                    delete[] denominator;

                    return  pow(-1,int(a+b+c+d)) * exp((D1+D2+D3+D4)/2) * sum;
                }
                else{
                    return 0;
                }
            }
            else{
                return 0;
            }
        }
        else{
            return 0;
        }
    }
    else{
        return 0;
    }

    return -1;
    //return  (Delta(a,b,e)*Delta(c,d,e)*Delta(a,c,f)*Delta(b,d,f))* sum;

}

double Racah_Coefficient(int a,int b, int c, int d, int e, int f){
    // For the moment, this function runs over integer values,
    // as L orbital angular momentum does.
    // To use the jj-coupled scheme or semi-integer values in the parameters,
    // the function has to be develop over
    //fractional objects of predefined class.
    
    // Condition for the strictly non negative values
    if((a<0) || (b<0) || (c<0) || (d<0)|| (e<0)|| (f<0)){return 0;}
  
    //evaluation of the deltas before the sum and the product speed the code up.
    double D1 = Delta(a,b,e);

    //Deltas are lower or equal to 0 when triangular condition is fulfilled
    if (D1<1){
        double D2 = Delta(d,c,e);//Delta(c,d,e); from devanathan
        if(D2<1){
            double D3 = Delta(a,c,f);
            if(D3<1){
                double D4 = Delta(d,b,f);//Delta(b,d,f); from devanathan
                if(D4<1){
                    double suma = 0.;
                    int c1, c2, c3, c4, c5, c6, c7;
                    // x higher or equal to (and -1)
                    c1 = a + b + e;
                    c2 = c + d + e;
                    c3 = a + c + f;
                    c4 = b + d + f;
                    // x lower or equal to
                    c5 = a + b + c + d;
                    c6 = a + d + e + f;
                    c7 = b + c + e + f;

                    // sumatory goes from the lower of the "lower than" to the highest of "higher than"
                    int major = min(min(c5,c6),c7);
                    // if x is lower than the lowest value of the "lower than" coefficients, then it is for the rest
                    int minor = max(max(max(max(c1,c2),c3),c4),-1) ;
                    // if x is bigger than the biggest "higher than" coefficient, then it is for the rest

                    //std::cout << "x max=" << mayor << "  x min" << minor << std::endl;

                    int *numerator;
                    numerator = new int[1];
                    int *denominator;
                    denominator = new int[7];
                    double aux;

                    for(int x=minor; x<=major;x++){
                        numerator[0] =  x + 1;
                        denominator[0] = x - c1;
                        denominator[1] = x - c2;
                        denominator[2] = x - c3;
                        denominator[3] = x - c4;
                        denominator[4] = c5 - x;
                        denominator[5] = c6 - x;
                        denominator[6] = c7 - x;

                        aux = factorial_function(numerator, denominator,1,7);
                        suma += pow(-1,x) * exp(aux);

                    }
                    delete[] numerator;
                    delete[] denominator;

                    return  pow(-1,a+b+c+d) * exp((D1+D2+D3+D4)/2) * suma;
                }
                else{
                    return 0;
                }
            }
            else{
                return 0;
            }
        }
        else{
            return 0;
        }
    }
    else{
        return 0;
    }

    return -1;
    //return  (Delta(a,b,e)*Delta(c,d,e)*Delta(a,c,f)*Delta(b,d,f))* sum;

}

double Six_j_Coefficient(int j1,int j2, int j3, int J1, int J2, int J3){
// Suhonen definition
///Six_j_Coeficients={   j1,     j2,  j3=j12}
///                  { J1=j3,  J2=j,  J3=j23)
    
    // Condition for the strictly non negative values
    if((J1<0) || (J2<0) || (J3<0) || (j1<0)|| (j2<0)|| (j3<0)){return 0;}
    
    //cout << "{"<< a <<" "<<b<<" "<< e <<"}"<<endl;
    //cout << "{"<< d <<" "<<c<<" "<< f <<"}"<<" = "<<pow(-1,a+b+c+d)*result <<endl;
    return pow(-1,j1+j2+J2+J1) * Racah_Coefficient( j1, j2, J2, J1, j3, J3);
    // for checking use  SixJSymbol[{j1, j2, j3}, {j4, j5, j6}] in Wolfram AlphA

}

double Six_j_Coefficient(Fraction j1,Fraction j2, Fraction j3, Fraction J1, Fraction J2, Fraction J3){
// Suhonen definition
//Six_j_Coeficients={int j1,int j2, int j12}
//                  {int j3, int j, int j23)
    // Condition for the strictly non negative values
    if((J1<0) || (J2<0) || (J3<0) || (j1<0)|| (j2<0)|| (j3<0)){return 0;}
    
    //cout << "{"<< a <<" "<<b<<" "<< e <<"}"<<endl;
    //cout << "{"<< d <<" "<<c<<" "<< f <<"}"<<" = "<<pow(-1,a+b+c+d)*result <<endl;
    return pow(-1,int(j1+j2+J2+J1)) * Racah_Coefficient( j1, j2, J2, J1, j3, J3);
    // for checking use  SixJSymbol[{j1, j2, j3}, {j4, j5, j6}] in Wolfram AlphA

}

double U_Coefficient(int a,int b, int c, int d, int e, int f){
    return sqrt( (2*e+1) * (2*f+1) ) * Racah_Coefficient(a,b,c,d,e,f);
}

double U_Coefficient(Fraction a,Fraction b, Fraction c, Fraction d, Fraction e, Fraction f){
    return sqrt( double(2*e+1) * double(2*f+1) ) * Racah_Coefficient(a,b,c,d,e,f);
}

/*  NO VA BIEN*/

double Racah_Coeficient_Analitic(int a_int,int b_int, int c_int, int d_int, int e_int, int f_int){
    //checking function for the W(abcd,ef), from Devanathan's book appendix and already checked.
    std::cout<<"Racah Analitical"<<std::endl;
    if (!(((abs(a_int-b_int) <= c_int) && (c_int <= (a_int+b_int)))||
          ((abs(c_int-b_int) <= a_int) && (a_int <= (c_int+b_int)))|| ((abs(a_int-c_int) <= b_int) && (b_int <= (a_int+c_int))))){
        return 0;
    }
    /*if(Delta(a_int,b_int,e_int)<1){return 0*pow(-1,a_int+b_int+c_int+d_int);}
    if(Delta(c_int,d_int,e_int)<1){return 0*pow(-1,a_int+b_int+c_int+d_int);}
    if(Delta(a_int,c_int,f_int)<1){return 0*pow(-1,a_int+b_int+c_int+d_int);}
    if(Delta(b_int,d_int,f_int)<1){return 0*pow(-1,a_int+b_int+c_int+d_int);}*/

    if(e_int == 1){
        double b,d,f;           //Necessary data type conversion for the division and the return
        b = (double)b_int;
        d = (double)d_int;
        f = (double)f_int;

        if(c_int == d_int+1){
            if(a_int == b_int+1){
                return pow(-1,b+d-f) * sqrt(((b+d+f+3)*(b+d+f+2)*(b+d-f+2)*(b+d-f+1))/(4*(2*b+1)*(2*b+3)*(b+1)*(2*d+1)*(2*d+3)*(d+1)));
            }
            if(a_int == b_int){
                return pow(-1,b+d-f) * sqrt(((b+d+f+2)*(b+d-f+1)*(d-b+f+1)*(b-d+f))/(4*b*(2*b+1)*(b+1)*(2*d+1)*(d+1)*(2*d+3)));
            }
            if(a_int == b_int-1){
                return pow(-1,b+d-f) * sqrt(((b-d+f)*(b-d+f-1)*(d+f-b+2)*(d+f-b+1))/(4*(2*b+1)*(2*b-1)*b*(d+1)*(2*d+1)*(2*d+3) ));
            }
            else{
                cout<<" not analitic"<<endl;
                return Racah_Coefficient( a_int, b_int,  c_int,  d_int,  e_int,  f_int);
            }
        }
        if(c_int == d_int){
            if(a_int == b_int+1){
                return pow(-1,b+d-f) * sqrt(((b+d+f+2)*(b-d+f+1)*(b+d-f+1)*(d+f-b))/(4*(2*b+1)*(2*b+3)*(b+1)*(2*d+1)*d*(d+1)));
            }
            if(a_int == b_int){
                return pow(-1,b+d-f) * sqrt((b*(b+1)+(d*d+1)+(f*(f+1)))/(4*b*(2*b+1)*(b+1)*(2*d+1)*d*(d+1)));
            }
            if(a_int == b_int-1){
                return pow(-1,b+d-f-1) * sqrt(((b+d+f+1)*(b+d-f)*(b-d+f)*(d+f-b+1))/(4*(2*b+1)*(2*b-1)*b*(d+1)*d*(2*d+1) ));
            }
            else{
                cout<<" not analitic"<<endl;
                return Racah_Coefficient( a_int, b_int,  c_int,  d_int,  e_int,  f_int);
            }
        }
        if(c_int == d_int-1){
            if(a_int == b_int+1){
                return pow(-1,b+d-f) * sqrt(((-b+d+f)*(d+f-b-1)*(b-d+f+2)*(b-d+f+1)) /(4*d*(2*b+1)*(2*b+3)*(b+1)*(2*d-1)*(2*d+1)));
            }
            if(a_int == b_int){
                return pow(-1,b+d-f-1) * sqrt(((b+d+f+1)*(b+f-d+1)*(d+f-b)*(b+d-f)) /(4*b*d*(2*b+1)*(b+1)*(2*d+1)*(2*d-1)));
            }
            if(a_int == b_int-1){
                return pow(-1,b+d-f) * sqrt(((b+d+f+1)*(b+d+f)*(b+d-f)*(b+d-f-1)) /(4*b*d*(2*b+1)*(2*b-1)*(2*d+1)*(2*d-1)));
            }
            else{
                cout<<" not analitic"<<endl;
                return Racah_Coefficient( a_int, b_int,  c_int,  d_int,  e_int,  f_int);
            }
        }
        else{
            cout<<" not analitic"<<endl;
            return Racah_Coefficient( a_int, b_int,  c_int,  d_int,  e_int,  f_int);
        }
    }
    else{
        cout<<" not analitic"<<endl;
        return Racah_Coefficient( a_int, b_int,  c_int,  d_int,  e_int,  f_int);
    }

}


double Clebsh_Gordan(int j1, int j2, int j, int m1, int m2, int m){
    // Devanathan definition, only for integer j1, j2j j3 values
    
    if((j1 < 0) || (j2 < 0) || (j < 0)) {return 0;}
    //conservation of the third component
    if ((m-m1-m2)!=0){

        return 0;
    }
    // Triangular conditions
    else if((abs(m1)>j1) ||(abs(m2)>j2)|| (abs(m)>j)){
        return 0;
    }
    //else if((abs(j1-j2)>j) || (j > j1+j2)){
    else if(!triangular_condition(j1,j2,j)){
        return 0;
    }
    // Fulfilling the basic conditions ensure properties, like have positive arguments for the factorials AB
    else{

        int *Num;
        Num = new int[9];
        int *Den;
        Den = new int[1];

        Den[0] = j1 + j2 + j +1 ;
        // A: j1j2j dependent constant
        Num[0] = j1 + j2 - j;
        Num[1] = j + j1 - j2;
        Num[2] = j2 + j - j1;
        // B: j1j2j and m1m2m dependent constant
        Num[3] = j1 + m1;
        Num[4] = j1 - m1;
        Num[5] = j2 + m2;
        Num[6] = j2 - m2;
        Num[7] = j + m;
        Num[8] = j - m;

        double AB = factorial_function(Num,Den,9,1);

        int c4,c5;
        //int c1,c2,c3,c4,c5;

        // lower than
        //c1 = Num[0]; //j1 + j2 - j;
        //c2 = Num[4]; //j1 - m1;
        //c3 = Num[5]; //j2 + m2;

        //higher than (and than 0)
        c4 = j2 - j - m1; // in the sum these goes with -
        c5 = j1 + m2 - j; // in the sum these goes with -

        int major = max(max(c4,c5),0);
        // as x must be bigger than these, only the maximum coefficient fulfill the condition for the rest
        int minor = min(min(Num[0],Num[4]),Num[5]);
        // the same but searching the minimum value for x to be lower than the coefficients

        double sum = 0.;
        double aux;

        int *Num_C;
        Num_C = new int[1];
        Num_C[0] = 0;

        int *Den_C;
        Den_C = new int[6];

        for(int x=minor; x<= major; x++){
            Den_C[0] = Num[0] - x;
            Den_C[1] = Num[4] - x;
            Den_C[2] = Num[5] - x;
            Den_C[3] = x - c4;
            Den_C[4] = x - c5;
            Den_C[5] = x;

            aux = factorial_function(Num_C,Den_C,1,6);
            sum += exp(aux) * pow(-1,x);
        }
        delete[] Num;
        delete[] Num_C;
        delete[] Den;
        delete[] Den_C;

        return sqrt((2*j+1) * exp (AB)) * sum;
    }
    return 5;
}

double Clebsh_Gordan(Fraction j1, Fraction j2, Fraction j, Fraction m1, Fraction m2, Fraction m){
    // Devanathan definition, extended to Fractions
    
    if((j1 < 0) || (j2 < 0) || (j < 0)) {return 0;}
    
    //conservation of the third component
    if ((m1+m2)!= m){
        return 0;
    }
    // Triangular conditions
    else if((abs(m1)>j1) ||(abs(m2)>j2)|| (abs(m)>j)){
        return 0;
    }
    //else if((abs(j1-j2)>j) || (j > j1+j2)){
    else if(!triangular_condition(j1,j2,j)){
        return 0;
    }
    // Fulfilling the basic conditions ensure properties, like have positive arguments for the factorials AB
    else{

        int *Num;
        Num = new int[9];
        int *Den;
        Den = new int[1];

        Den[0] = int(j1 + j2 + j +1) ;
        // A: j1j2j dependent constant
        Num[0] = int(j1 + j2 - j);
        Num[1] = int(j + j1 - j2);
        Num[2] = int(j2 + j - j1);
        // B: j1j2j and m1m2m dependent constant
        Num[3] = int(j1 + m1);
        Num[4] = int(j1 - m1);
        Num[5] = int(j2 + m2);
        Num[6] = int(j2 - m2);
        Num[7] = int(j + m);
        Num[8] = int(j - m);

        //print_array(Num,9);
        //print_array(Den,1);
        double AB = factorial_function(Num,Den,9,1);

        int c4,c5;
        //int c1,c2,c3,c4,c5;

        // lower than
        //c1 = Num[0]; //j1 + j2 - j;
        //c2 = Num[4]; //j1 - m1;
        //c3 = Num[5]; //j2 + m2;

        //higher than (and than 0)
        c4 = int(j2 - j - m1); // in the sum these goes with -
        c5 = int(j1 + m2 - j); // in the sum these goes with -

        int major = max(max(c4,c5),0);
        // as x must be bigger than these, only the maximum coefficient fulfill the condition for the rest
        int minor = min(min(Num[0],Num[4]),Num[5]);
        // the same but searching the minimum value for x to be lower than the coefficients
        //cout<<"minor"<<minor<<"major"<<major<<endl;
        double sum = 0.;
        double aux;

        int *Num_C;
        Num_C = new int[1];
        Num_C[0] = 0;

        int *Den_C;
        Den_C = new int[6];

        for(int x = min(minor,major); x <= max(minor,major); x++){
            Den_C[0] = Num[0] - x;
            Den_C[1] = Num[4] - x;
            Den_C[2] = Num[5] - x;
            Den_C[3] = x - c4;
            Den_C[4] = x - c5;
            Den_C[5] = x;
            //print_array(Den_C,6);

            aux = factorial_function(Num_C,Den_C,1,6);
            //cout<<"aux ff="<<aux<<endl;
            sum += exp(aux) * pow(-1,x);
        }
        delete[] Num;
        delete[] Num_C;
        delete[] Den;
        delete[] Den_C;
        //cout<<"suma="<<sum<<endl;
        return sqrt((int(2*j)+1) * exp (AB)) * sum;
    }
    return 5;
}

double Three_j_Symbol(int j1, int j2, int j, int m1, int m2, int m){
    // only for integers
    // The true condition of the 3j is not the same for CCCG
    if((m1 + m2) != -m){return 0;}
    
    if((j1 < 0) || (j2 < 0) || (j < 0)) {
      std::cout << "Negative " << std::endl;
      return 0;}
    
    return pow(-1, j1 - j2 - m) * Clebsh_Gordan(j1,j2,j,m1,m2,-m) / sqrt(2*j+1);   // Mathematica definition
}

double Three_j_Symbol(Fraction j1, Fraction j2, Fraction j, Fraction m1, Fraction m2, Fraction m){
    // only for integers
    if((m1 + m2) != Fraction(-m.numerator,m.denominator)){return 0;}
    
    // Condition for the strictly non negative values
    if((j1 < 0) || (j2 < 0) || (j < 0)) {
      std::cout << "Negative " << std::endl;
      return 0;}
    
    return pow(-1, int(((j1 - j2) - m))) * Clebsh_Gordan(j1,j2,j,m1,m2,Fraction(-m.numerator, m.denominator) ) / sqrt(2*j+1);   // Mathematica definition
}


double Nine_j_Symbol(Fraction l1, Fraction s1, Fraction j1, Fraction l2, Fraction s2, Fraction j2, Fraction L, Fraction S, Fraction J){
    // definition from Devanathan
    // {l1, s1, j1}
    // {l2, s2, j2}
    // {L,  S,  J }
    
    // Condition for the strictly non negative values
    if((l1<0) || (s1<0) || (j1<0) || (l2<0)|| (s2<0)|| (j2<0)|| (L<0)|| (S<0)|| (J<0)){return 0;}
    
    if( triangular_condition(l1,s1,j1) && triangular_condition(l2,s2,j2) && triangular_condition( L, S, J) &&
	triangular_condition(l1,l2, L) && triangular_condition(s1,s2, S) && triangular_condition(j1,j2, J)){
	double sum = 0.;
    
	int minimum = 0; 
	int maximum = min(min( (l1+J) ,(l2+S)), (s1+j2));
	//Fraction minimum = 0;
	//Fraction maximum = min(min( (l1+J) ,(l2+S)), (s1+j2));
	for( int t=minimum; t<= maximum; t++){
	    sum += (2*t + 1)* Racah_Coefficient(l1,l2,J,S,L,Fraction(t))*
		Racah_Coefficient(l1,s1,J,j2,j1,Fraction(t))* Racah_Coefficient(s1,S,j2,l2,s2,Fraction(t));
	}
	return sum;
    }
    else{return 0;}

}
/*
double Nine_j_Symbol(Fraction l1, Fraction s1, Fraction j1, Fraction l2, Fraction s2, Fraction j2, Fraction L, Fraction S, Fraction J){
    // definition from Devanathan
    // {l1, s1, j1}
    // {l2, s2, j2}
    // {L,  S,  L }
    
    // Condition for the strictly non negative values
    if((l1<0) || (s1<0) || (j1<0) || (l2<0)|| (s2<0)|| (j2<0)|| (L<0)|| (S<0)|| (J<0)){return 0;}

    double sum = 0.;
    Fraction minimum = 0;
    Fraction maximum = min(min( (l1+J) ,(l2+S)), (s1+j2));
    for( int t=minimum; t<= maximum; t++){
        sum += (2*t + 1)* Racah_Coefficient(l1,l2,J,S,L,Fraction(t))*
            Racah_Coefficient(l1,s1,J,j2,j1,Fraction(t))* Racah_Coefficient(s1,S,j2,l2,s2,Fraction(t));
    }
    return sum;

}// */

double Nine_j_Symbol(int l1, int s1, int j1, int l2, int s2, int j2, int L, int S, int J){
    // definition from Devanathan, for this program ss must be integer non negative values
    // {l1, s1, j1}
    // {l2, s2, j2}
    // {L,  S,  J }
    
    // Condition for the strictly non negative values
    if((l1<0) || (s1<0) || (j1<0) || (l2<0)|| (s2<0)|| (j2<0)|| (L<0)|| (S<0)|| (J<0)){return 0;}
    
    if( triangular_condition(l1,s1,j1) && triangular_condition(l2,s2,j2) && triangular_condition( L, S, J) &&
	triangular_condition(l1,l2, L) && triangular_condition(s1,s2, S) && triangular_condition(j1,j2, J)){
      
	double sum = 0.;
	int minimum = 0; 
	int maximum = min(min( (l1+J) ,(l2+S)), (s1+j2));
	//int minimum = max(max(max( abs(l1-J) ,abs(l2-S)), abs(s1-j2)), abs(S-J));;
	//int maximum = min(min(min( (l1+J) ,(l2+S)), (s1+j2)), (S+J));
	for( int t=minimum; t<= maximum; t++){
	    sum += (2*t + 1)* Racah_Coefficient(l1,l2,J,S,L,t)*
		Racah_Coefficient(l1,s1,J,j2,j1,t)* Racah_Coefficient(s1,S,j2,l2,s2,t);
	}
	return sum;
    }
    else{return 0;}
}



double LS_jj_coupling_Coeff(Fraction l1, Fraction s1, Fraction j1, Fraction l2, Fraction s2, Fraction j2, Fraction L, Fraction S, Fraction J){
    
    // definition from Devanathan, differs from 9j on a multiply factor
    // {l1, s1, j1}
    // {l2, s2, j2}
    // {L,  S,  L }

    // Condition for the strictly non negative values
    if((l1<0) || (s1<0) || (j1<0) || (l2<0)|| (s2<0)|| (j2<0)|| (L<0)|| (S<0)|| (J<0)){return 0;}
    
    /*
    double sum = 0.;
    Fraction minimum = 0;
    Fraction maximum = min(min( (l1+J) ,(l2+S)), (s1+j2));
    for( int t=minimum; t<= maximum; t++){
        sum += (2*t + 1)* Racah_Coefficient(l1,l2,J,S,L,Fraction(t))*
            Racah_Coefficient(l1,s1,J,j2,j1,Fraction(t))* Racah_Coefficient(s1,S,j2,l2,s2,Fraction(t));
    }*/
    
    
    return sqrt((2*j1 + 1) * (2*j2 + 1) * (2*L + 1) * (2*S + 1)) * Nine_j_Symbol(l1,s1,j1, l2,s2,j2, L,S,J);

}

double LS_jj_coupling_Coeff(int l1, int s1, int j1, int l2, int s2, int j2, int L, int S, int J){
    
    // definition from Devanathan, for this program ss must be integer non negative values
    // {l1, s1, j1}
    // {l2, s2, j2}
    // {L,  S,  L }
    
    // Condition for the strictly non negative values
    if((l1<0) || (s1<0) || (j1<0) || (l2<0)|| (s2<0)|| (j2<0)|| (L<0)|| (S<0)|| (J<0)){return 0;}

    /*
    double sum = 0.;
    int minimum = 0;
    int maximum = min(min( (l1+J) ,(l2+S)), (s1+j2));
    for( int t=minimum; t<= maximum; t++){
        sum += (2*t + 1)* Racah_Coefficient(l1,l2,J,S,L,t)*
            Racah_Coefficient(l1,s1,J,j2,j1,t)* Racah_Coefficient(s1,S,j2,l2,s2,t);
    }*/
    return sqrt((2*j1 + 1) * (2*j2 + 1) * (2*L + 1) * (2*S + 1)) * Nine_j_Symbol(l1,s1,j1, l2,s2,j2, L,S,J);

}

