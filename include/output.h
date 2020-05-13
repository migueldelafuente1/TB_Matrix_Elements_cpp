#ifndef OUTPUT_H_INCLUDED
#define OUTPUT_H_INCLUDED

std::string IntToString (int a);
int index_antoine(QN_1body_jj SET);
std::string name_index_antoine(QN_1body_jj SET);
std::string name_index_antoine(QN_2body_jj_Coupling WF);

void antoine_output(int n);

int shell_identification (int *SHELLS, QN_1body_jj *VALENCE, int Number_of_distint_elements);
int core_shell(QN_1body_jj *VALENCE, int Number_of_distint_elements);
int core_number( int core_shell);
#endif // OUTPUT_H_INCLUDED
