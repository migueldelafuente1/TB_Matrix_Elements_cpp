#ifndef INDEX_COEFFICIENTS_H_INCLUDED
#define INDEX_COEFFICIENTS_H_INCLUDED
//////////            CONSTANTS           //////////////
extern double H_BAR;
extern double H_BAR_C;
extern double M_PROTON ;
extern int A;

double b_lenght(int A);
double h_bar_omega(int A);
double R_radius(int A);
extern double V_0;
extern double R_interaction;

extern bool Set_Manual_Lambda;

//////////  MONSHINSKY NECESARY COEFFICIENTS  //////////
// Constant coefficient in BMB_00, for extend factorial evaluations
double A_coefficient(int l1, int l, int l2, int L, int x);

//double M_E_minus_square_r_i(int i ,int lambda,int n,int l,int N,int L,int n_q,int l_q,int N_q,int L_q);

// Recurrence Matrix element for evaluation of general BMB
double ME_ri2_times_BMB(int i ,int lambda,int n,int l,int N,int L,int n_q,int l_q,int N_q,int L_q, int n1, int l1, int n2, int l2);
// Coefficient for the decomposition of the
// reduced M.Element < nl||V(r)||n'l'> in Talmi integrals
double B_coefficient(int n, int l, int n_q, int l_q, int p, double b_param);


/////////// STRUCT TO PASS QUANTUM NUMBERS /////////////
// * Can be initialized as  Struct Object = {<item 1>, <item 2>, ...} in standar C++11
// * left without indications on declarations lead 0, but it reads following the order,
// * It can pass directly the parameters between them just with: Object1 = Object2;

struct QN_2body_jj_Coupling{
    // = {n1,l1,j1, n2,l2,j2, J,M,T,MT};
    int n1;
    int l1;
    Fraction j1;
                    //Fraction s1;
    int n2;
    int l2;
    Fraction j2;
                    //Fraction s2;
    Fraction J;
    Fraction M;
    int T;
    int M_T;
};

struct QN_1body_jj{
	// = {n, l, j m};
    int n;
    int l;
    Fraction j;
    Fraction m;
                    //Fraction s1;
};


struct QN_2body_radial{
	// = {n1,l1,n2,l2,lamdba, mu};
    int n1;
    int l1;
    int n2;
    int l2;
    int lambda;
    int mu;
};

struct QN_1body_radial{
	// = {n, l, m_l}
    int n;
    int l;
    int m_l;
};

int delta(int left, int right);

int delta_2_1body_j(const QN_1body_jj &WF1,const QN_1body_jj &WF2);

double Normalization_JT(const QN_2body_jj_Coupling &WF);

////////// SUHONEN INTERMEDIATE COEFFICIENTS ////////////
///// IN MULTIPOLE DECOMPOSITION FOR M.E          ///////

double Y_Coeff(const QN_1body_jj &left, const QN_1body_jj &right, int lambda, Fraction j);
double Omega_Coeff(const QN_2body_jj_Coupling &left, const QN_2body_jj_Coupling &right, int lambda, Fraction J);

double Z_Coeff(const QN_1body_jj &left, const QN_1body_jj &right, int lambda, Fraction j);
double Lambda_Coeff(const QN_2body_jj_Coupling &left, const QN_2body_jj_Coupling &right, int lambda, Fraction J);

double U_Coeff(const QN_2body_jj_Coupling &left, const QN_2body_jj_Coupling &right, int lambda, double *Parameters);

///// GAUSS-LEGENGRE NODES AND WEIGHTS FOR QUADRATURES /////

extern double ** Gauss_Legendre_values;
void run_Gauss_Legendre(int ORDER);
double Gauss_Legendre_Weights (int ORDER, int k); 
double Gauss_Legendre_Nodes (int ORDER, int k);

////////// SUHONEN INTERMEDIATE COEFFICIENTS ////////////
/////        FOR SURFACE DELTA INTERACTION        ///////

double K_Coeff(const QN_2body_jj_Coupling &left,const  QN_2body_jj_Coupling &right);


//////////    TALMAN ARTICLE COEFFICIENTS,     ////////////
/////    EVALUATION OF RADIAL INTEGRALS FOR M.E     ///////
double c_normalization_Tal(int n, int l);

double Big_C_coeff_Tal(const QN_2body_radial & Arguments, int S, int L);

double R_coeff_Tal( int S_1, int S_2, int L, int nu);

#endif // INDEX_COEFFICIENTS_H_INCLUDED
