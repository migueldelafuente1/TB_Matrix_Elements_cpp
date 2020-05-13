#include "../include/Fractions.h"
#include <iostream>
#include <cmath>
#include <fstream>
#include <string>
#include <sstream>
#include <stdlib.h>

using namespace std;

// Calculates the greates common divisor with
// Euclid's algorithm
// both arguments have to be positive
int Fraction::gcd(int a, int b)
{
    while (a != b) {
        if (a > b) {
            a -= b;
        } else {
            b -= a;
        }
    }
    return a;
}

Fraction::~Fraction(){}

Fraction::Fraction() {
    numerator = 0;
    denominator = 1;
}
Fraction::Fraction(int n) {
    numerator = n;
    denominator = 1;
}

Fraction::Fraction(int n, int d) {
    if (d==0) {
        cerr << "Denominator may not be 0." << endl;
        exit(0);
    } else if (n == 0) {
        numerator = 0;
        denominator = 1;
    } else {
        int sign = 1;
        if (n < 0) {
            sign *= -1;
            n *= -1;
        }
        if (d < 0) {
            sign *= -1;
            d *= -1;
        }

        int tmp = gcd(n, d);
        numerator = n/tmp*sign;
        denominator = d/tmp;
    }
}
void Fraction::irreducible(){
    Fraction tmp(numerator,denominator);
    numerator = tmp.numerator;
    denominator = tmp.denominator;
}

Fraction::operator int() {return (numerator)/denominator;}
Fraction::operator float() {return ((float)numerator)/denominator;}
Fraction::operator double() {return ((double)numerator)/denominator;}

Fraction operator+(const Fraction& lhs, const Fraction& rhs) {
    Fraction tmp(lhs.numerator*rhs.denominator
                +rhs.numerator*lhs.denominator,
                lhs.denominator*rhs.denominator);
    return tmp;
}
Fraction operator+(const Fraction& lhs, const int& rhs) {
    Fraction tmp(lhs.numerator+rhs*lhs.denominator,
                lhs.denominator);
    return tmp;
}
Fraction operator+(const int& lhs, const Fraction& rhs) {
    Fraction tmp(rhs.numerator+lhs*rhs.denominator,
                rhs.denominator);
    return tmp;
}

void operator +=(Fraction& lhs, const Fraction& rhs) {
    Fraction tmp(lhs.numerator*rhs.denominator
                +rhs.numerator*lhs.denominator,
                lhs.denominator*rhs.denominator);
    lhs.numerator = tmp.numerator;
    lhs.denominator = tmp.denominator;
}
void operator +=(Fraction& lhs, const int& rhs) {
    Fraction tmp(lhs.numerator + rhs*lhs.denominator,
                lhs.denominator);
    lhs.numerator = tmp.numerator;
    lhs.denominator = tmp.denominator;
}

/*Fraction operator++(Fraction& lhs) {
    Fraction tmp(lhs.numerator+lhs.denominator, lhs.denominator);
    return tmp;
    //lhs.numerator = tmp.numerator;
    //lhs.denominator = tmp.denominator;
}
 */

Fraction operator-(const Fraction& lhs, const Fraction& rhs) {
    Fraction tmp(lhs.numerator*rhs.denominator
                -rhs.numerator*lhs.denominator,
                lhs.denominator*rhs.denominator);
    return tmp;
}
Fraction operator-(const Fraction& lhs, const int& rhs) {
    Fraction tmp(lhs.numerator-rhs*lhs.denominator,
                lhs.denominator);
    return tmp;
}
Fraction operator-(const int& lhs, const Fraction& rhs) {
    Fraction tmp(lhs*rhs.denominator-rhs.numerator,
                rhs.denominator);
    return tmp;
}

void operator-=(Fraction& lhs, const Fraction& rhs) {
    Fraction tmp(lhs.numerator*rhs.denominator
                -rhs.numerator*lhs.denominator,
                lhs.denominator*rhs.denominator);
    lhs.numerator = tmp.numerator;
    lhs.denominator = tmp.denominator;
}
void operator-=(Fraction& lhs, const int& rhs) {
    Fraction tmp(lhs.numerator -rhs*lhs.denominator ,
                lhs.denominator);
    lhs.numerator = tmp.numerator;
    lhs.denominator = tmp.denominator;
}

/*void operator--(Fraction& lhs) {
    Fraction tmp(lhs.numerator-lhs.denominator,
                 lhs.denominator);
    lhs.numerator = tmp.numerator;
    lhs.denominator = tmp.denominator;
}
 */

Fraction operator*(const Fraction& lhs, const Fraction& rhs) {
    Fraction tmp(lhs.numerator*rhs.numerator,
               lhs.denominator*rhs.denominator);
    return tmp;
}
Fraction operator*(int lhs, const Fraction& rhs) {
    Fraction tmp(lhs*rhs.numerator,rhs.denominator);
    return tmp;
}
Fraction operator*(const Fraction& rhs, int lhs) {
    Fraction tmp(lhs*rhs.numerator,rhs.denominator);
    return tmp;
}

void operator*=(Fraction& lhs, const Fraction& rhs) {
    Fraction tmp(lhs.numerator*rhs.numerator,
               lhs.denominator*rhs.denominator);
    lhs.numerator = tmp.numerator;
    lhs.denominator = tmp.denominator;
    //return lhs;
}
void operator*=(Fraction& lhs, const int& rhs) {
    Fraction tmp(lhs.numerator*rhs,lhs.denominator);
    lhs.numerator = tmp.numerator;
    lhs.denominator = tmp.denominator;
    //return lhs;
}


Fraction operator/(const Fraction& lhs, const Fraction& rhs) {
    Fraction tmp(lhs.numerator*rhs.denominator,
                 lhs.denominator*rhs.numerator);
    return tmp;
}
Fraction operator/(const int& lhs, const Fraction& rhs) {
    Fraction tmp(lhs*rhs.denominator,rhs.numerator);
    return tmp;
}
Fraction operator/(const Fraction& lhs, const int& rhs) {
    Fraction tmp(lhs.numerator, lhs.denominator*rhs);
    return tmp;
}

void operator/=(Fraction& lhs, const Fraction& rhs) {
    Fraction tmp(lhs.numerator*rhs.denominator,
               lhs.denominator*rhs.numerator);
    lhs.numerator = tmp.numerator;
    lhs.denominator = tmp.denominator;
    //return lhs;
}
void operator/=(Fraction& lhs, const int& rhs) {
    Fraction tmp(lhs.numerator,lhs.denominator*rhs);
    lhs.numerator = tmp.numerator;
    lhs.denominator = tmp.denominator;
    //return lhs;
}


bool operator!=(const Fraction& lhs, const Fraction& rhs) {
    if((lhs.numerator!=rhs.numerator) || (lhs.denominator!=rhs.denominator)){
        return true;
    }
    else{
        return false;
    }
}
bool operator!=(const int& lhs, const Fraction& rhs) {
    if((lhs!=rhs.numerator) || (rhs.denominator!=1)){
        return true;
    }
    else{
        return false;
    }
}
bool operator!=(const Fraction& lhs, const int& rhs) {
    if((lhs.numerator!=rhs) || (lhs.denominator!=1)){
        return true;
    }
    else{
        return false;
    }
}

bool operator==(const Fraction& lhs, const Fraction& rhs) {
    if((lhs.numerator==rhs.numerator) && (lhs.denominator==rhs.denominator)){
        return true;
    }
    else{
        return false;
    }
}
bool operator==(const int& lhs, const Fraction& rhs) {
    if((lhs==rhs.numerator) && (rhs.denominator==1)){
        return true;
    }
    else{
        return false;
    }
}
bool operator==(const Fraction& lhs, const int& rhs) {
    if((lhs.numerator==rhs) && (lhs.denominator==1)){
        return true;
    }
    else{
        return false;
    }
}

bool operator <=(const Fraction& lhs, const Fraction& rhs){
    if((lhs.numerator*rhs.denominator)<=(rhs.numerator*lhs.denominator)){
        return true;
    }
    else{return false;}
}
bool operator <=(const Fraction& lhs, const int& rhs){
    if((lhs.numerator)<=(rhs*lhs.denominator)){
        return true;
    }
    else{return false;}
}
bool operator <=(const int& lhs, const Fraction& rhs){
    if((lhs*rhs.denominator)<=(rhs.numerator)){
        return true;
    }
    else{return false;}
}

bool operator <(const Fraction& lhs, const Fraction& rhs){
    if((lhs.numerator*rhs.denominator)<(rhs.numerator*lhs.denominator)){
        return true;
    }
    else{return false;}
}
bool operator <(const Fraction& lhs, const int& rhs){
    if((lhs.numerator)<(rhs*lhs.denominator)){
        return true;
    }
    else{return false;}
}
bool operator <(const int& lhs, const Fraction& rhs){
    if((lhs*rhs.denominator)<(rhs.numerator)){
        return true;
    }
    else{return false;}
}

bool operator >=(const Fraction& lhs, const Fraction& rhs){
    if((lhs.numerator*rhs.denominator)>=(rhs.numerator*lhs.denominator)){
        return true;
    }
    else{return false;}
}
bool operator >=(const Fraction& lhs, const int& rhs){
    if((lhs.numerator)>=(rhs*lhs.denominator)){
        return true;
    }
    else{return false;}
}
bool operator >=(const int& lhs, const Fraction& rhs){
    if((lhs*rhs.denominator)>=(rhs.numerator)){
        return true;
    }
    else{return false;}
}

bool operator >(const Fraction& lhs, const Fraction& rhs){
    if((lhs.numerator*rhs.denominator)>(rhs.numerator*lhs.denominator)){
        return true;
    }
    else{return false;}
}
bool operator >(const Fraction& lhs, const int& rhs){
    if((lhs.numerator)>(rhs*lhs.denominator)){
        return true;
    }
    else{return false;}
}
bool operator >(const int& lhs, const Fraction& rhs){
    if((lhs*rhs.denominator)>(rhs.numerator)){
        return true;
    }
    else{return false;}
}

Fraction abs(const Fraction frac){
    int sign_n = 1;
    int sign_d = 1;
    if (frac.numerator < 0) {
        sign_n *= -1;
    }
    if (frac.denominator < 0) {
        sign_d *= -1;
    }
    return Fraction(frac.numerator*sign_n , frac.denominator*sign_d);
}
/*
double pow(Fraction base, Fraction exponent){
    if(base.numerator<0){
        std::cout<<" complex "<<std::endl;
        exit(0);
    }
    else{
        return  pow(base.numerator,double(exponent)),pow(base.denominator,double(exponent));
    }
}
*/

Fraction min(const Fraction lhs, const Fraction rhs){
    //std::cout << "de Fracciones" << std::endl;
    if((lhs.numerator*rhs.denominator)<(rhs.numerator*lhs.denominator)){
        return lhs;
    }
    else{
        return rhs;
    }
}

Fraction max(const Fraction lhs, const Fraction rhs){
    //std::cout << "de Fracciones" << std::endl;
    if((lhs.numerator*rhs.denominator)>(rhs.numerator*lhs.denominator)){
        return lhs;
    }
    else{
        return rhs;
    }
}

std::ostream& operator<<(std::ostream &strm, const Fraction &a) {
    // to run properly, lhs Fraction value doesn't be an operation.
    if (a.denominator == 1) {
        strm << a.numerator;
    } else {
        strm << a.numerator << "/" << a.denominator;
    }
    return strm;
}

void print_array(int *Arr, int lenght){
    for(int i=0; i<lenght;i++){
        std::cout<< Arr[i]<<std::endl;
    }
}

void print_array(Fraction *Arr, int lenght){
    for(int i=0; i<lenght;i++){
        std::cout<< Arr[i]<<std::endl;
    }
}

double sqrt(Fraction x){
    return sqrt(double(x));
}
