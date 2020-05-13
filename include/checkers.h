#ifndef CHECKERS_H_INCLUDED
#define CHECKERS_H_INCLUDED

void check_factorial_function();

void check_Racah1(int b, int d,int f_min,int f_max);
void check_Racah2(int N, bool half_integers);

void check_CCGG(int N, bool with_halfintegers);
void check_9j(int N, bool with_halfintegers);
void check_LSjj_coupling_coeffs(QN_2body_jj_Coupling BRA, QN_2body_jj_Coupling KET);

void check_BM_00();

void BMB_symmetry_checking(int n, int l, int N, int L, int lambda, int n1, int l1, int n2, int l2);

void check_BMB_article_values();

void check_BMB(int N_input);

void BMB_ortonormality_checking( int n1, int l1, int n2,int l2, int lambda, int n1_q, int l1_q, int n2_q,int l2_q );


void check_ME_Art(int N, int n_Max);

void check_ME_multipolar(QN_1body_jj WF, int order);

void check_ME_SDI();

void print_QN_2body_jj ( QN_2body_jj_Coupling WF);
void print_QN_1body_jj ( QN_1body_jj WF);
void print_QN_2body_rad( QN_2body_radial WF);
void print_QN_1body_rad( QN_1body_radial WF);

bool check_Chasman_Conversion_Spherical2Cartesian (QN_1body_radial Spherical);
void check_SHO_WF( QN_1body_jj & WF, double R_MAX, double b_param );
void check_Cuadrature_Convergence(QN_2body_radial Q_N_radial_BRA,QN_2body_radial Q_N_radial_KET,double lambda,
				      double mu_param,double b_param,int Max_order, int precision);

double check_Coupling_BBME_Methods(QN_2body_jj_Coupling Q_Numbers_left, QN_2body_jj_Coupling Q_Numbers_right,
					  double mu_param,double b_param,int Max_order, int precision, bool Moshinsky_Method);
#endif // CHECKERS_H_INCLUDED
