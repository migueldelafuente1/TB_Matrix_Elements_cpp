#include <iostream>

#ifndef FRACTIONS_H_INCLUDED
#define FRACTIONS_H_INCLUDED
class Fraction {
    private:
        int gcd(int a, int b);

    public:
        int numerator, denominator;

        Fraction();
        ~Fraction();
        Fraction(int n) ;
        Fraction(int n, int d) ;
    	void irreducible();

        operator int() ;
        operator float() ;
        operator double() ;

};
Fraction operator+(const Fraction& lhs, const Fraction& rhs);
//addition of an integer
Fraction operator+(const Fraction& lhs, const int& rhs) ;
Fraction operator+(const int& lhs, const Fraction& rhs);
void operator+=(Fraction& lhs, const Fraction& rhs) ;
void operator+=(Fraction& lhs, const int& rhs) ;
//Fraction operator++(Fraction& lhs);// since c++20


Fraction operator-(const Fraction& lhs, const Fraction& rhs);
Fraction operator-(const Fraction& lhs, const int& rhs);
Fraction operator-(const int& lhs, const Fraction& rhs);
void operator-=(Fraction& lhs, const Fraction& rhs) ;
void operator-=(Fraction& lhs, const int& rhs) ;
//void operator--(Fraction& lhs); // since c++20

Fraction operator*(const Fraction& lhs, const Fraction& rhs) ;
Fraction operator*(int lhs, const Fraction& rhs) ;
Fraction operator*(const Fraction& rhs, int lhs) ;
void operator*=(Fraction& lhs, const Fraction& rhs);
void operator*=(Fraction& lhs, const int& rhs);

Fraction operator/(const Fraction& lhs, const Fraction& rhs);
Fraction operator/(const int& lhs, const Fraction& rhs);
Fraction operator/(const Fraction& lhs, const int& rhs);

void operator/=(Fraction& lhs, const Fraction& rhs);
void operator/=(Fraction& lhs, const int& rhs);

bool operator!=(const Fraction& lhs, const Fraction& rhs) ;
bool operator!=(const int& lhs, const Fraction& rhs) ;
bool operator!=(const Fraction& lhs, const int& rhs) ;
bool operator==(const Fraction& lhs, const Fraction& rhs) ;
bool operator==(const int& lhs, const Fraction& rhs) ;
bool operator==(const Fraction& lhs, const int& rhs) ;

bool operator <=(const Fraction& lhs, const Fraction& rhs);
bool operator <=(const Fraction& lhs, const int& rhs);
bool operator <=(const int& lhs, const Fraction& rhs);

bool operator<(const Fraction& lhs, const Fraction& rhs);
bool operator<(const Fraction& lhs, const int& rhs);
bool operator<(const int& lhs, const Fraction& rhs);

bool operator>=(const Fraction& lhs, const Fraction& rhs);
bool operator>=(const Fraction& lhs, const int& rhs);
bool operator>=(const int& lhs, const Fraction& rhs);

bool operator>(const Fraction& lhs, const Fraction& rhs);
bool operator>(const Fraction& lhs, const int& rhs);
bool operator>(const int& lhs, const Fraction& rhs);


Fraction abs(const Fraction frac);
//double pow(Fraction frac, Fraction exponent);
Fraction min(const Fraction lhs, const Fraction rhs);
Fraction max(const Fraction lhs, const Fraction rhs);
std::ostream& operator<<(std::ostream &strm, const Fraction &a);

void print_array(int *Arr, int lenght);
void print_array(Fraction *Arr, int lenght);

double sqrt(Fraction x);
#endif // FUNCTIONS_H_INCLUDED
