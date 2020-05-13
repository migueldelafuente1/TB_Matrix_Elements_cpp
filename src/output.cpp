#include "../include/factorials.h"
#include "../include/Index_Coefficients.h"
#include "../include/Angular_Functions.h"
#include "../include/checkers.h"
#include "../include/BM_Brackets.h"
#include "../include/Fractions.h"
#include "../include/Integrals.h"

#include "../include/Nuclear_Matrix_Elements_M.h"
#include "../include/Nuclear_Matrix_Elements_Art.h"
#include "../include/Nuclear_Matrix_Elements_Suh.h"
#include "../include/Nuclear_Matrix_Elements_BB.h"

#include <iostream>
#include <fstream>
#include <math.h>
#include <cmath>
#include <iomanip>
#include <stdlib.h>
#include <time.h>
#include <string>
#include <sstream>

#define WRITESTRING(x) #x

std::string IntToString (int a){
    // Converting integers to string, if the compiler accept the flag -std=c++11, it can be used << to_string()
    std::ostringstream temp;
    temp<<a;
    return temp.str();
}

int index_antoine(QN_1body_jj SET){
    // Antoine's way of reading the quantum numbers
    return 1000 * SET.n + 100 * SET.l + int(2 * SET.j);
}

std::string name_index_antoine(QN_1body_jj SET){
    // This function is necessary to return the antoine index for the state 0s_{1/2}
    // as 001 and not just by 1.
    if(index_antoine(SET) == 1){
        return WRITESTRING(001);
    }
    else{
        return IntToString(index_antoine(SET));
    }
}

std::string name_index_antoine(QN_2body_jj_Coupling WF){
    // This function is necessary to return the antoine index for the state 0s_{1/2}
    // as 001 and not just by 1.
    
    QN_1body_jj SET_1 = {WF.n1, WF.l1, WF.j1, 0};
    QN_1body_jj SET_2 = {WF.n2, WF.l2, WF.j2, 0};
    
    if( (index_antoine(SET_1) == 1) || (index_antoine(SET_2) == 1)){
	if( (index_antoine(SET_1) == 1) && (index_antoine(SET_2) == 1)){
	    return WRITESTRING(001   001);
	}
	else if(index_antoine(SET_1) == 1){
	    return WRITESTRING(001   ) + IntToString(index_antoine(SET_2));
	}
	else if(index_antoine(SET_2) == 1){
	    return IntToString(index_antoine(SET_2)) + "   " + WRITESTRING(   001) ;
	}
    }

    else{
        return IntToString(index_antoine(SET_1)) + WRITESTRING(.  .) + IntToString(index_antoine(SET_2));
    }
}

int shell_identification (int *SHELLS, QN_1body_jj *VALENCE, int Number_of_distint_elements){

    ///     This process  the different shells that appear in the valence space
    /// according to the normal levels of the shell model.
    ///     The return is the lower level of the shells, in order to indicate the core.
    int N, L;
    Fraction J;

    // Lower shell
    int Lower_Shell = 1e5; // arbitrary high

    for( int i = 0; i<Number_of_distint_elements; i++){
        N = VALENCE[i].n;
        L = VALENCE[i].l;
        J = VALENCE[i].j;

        if(N == 0){
            if (L==0){ // S shell
                SHELLS[i] = 1;
            }
            else if ( L==1){ // P shell
                SHELLS[i] = 2;
            }
            else if (L==2){ // D
                SHELLS[i] = 3;
            }
            else if (L==3){ // F
                if(J == Fraction(7,2)){
                    SHELLS[i] = 4;
                }
                else if(J == Fraction(5,2)){
                    SHELLS[i] = 5;
                }
            }
            else if (L==4){ // G
                if(J == Fraction(9,2)){
                    SHELLS[i] = 5;
                }
                else if(J == Fraction(7,2)){
                    SHELLS[i] = 6;
                }
            }
            else if (L==5){ // H
                if(J == Fraction(11,2)){
                    SHELLS[i] = 6;
                }
                else if(J == Fraction(9,2)){
                    SHELLS[i] = 7;
                }
            }
            else if (L==6){ // I
                if(J == Fraction(13,2)){
                    SHELLS[i] = 7;
                }
                else if(J == Fraction(11,2)){
                    SHELLS[i] = 8;
                }
            }
            else if (L==7){ // J
                if(J == Fraction(15,2)){
                    SHELLS[i] = 8;
                }
                else if(J == Fraction(13,2)){
                    SHELLS[i] = 9;
                }
            }
        }
        else if(N == 1){
            if (L==0){      //S
                SHELLS[i] = 3;
            }
            else if ( L==1){ // P
                SHELLS[i] = 5;
            }
            else if (L==2){ // D
                SHELLS[i] = 6;
            }
            else if (L==3){ // F
                SHELLS[i] = 7;
            }
            else if (L==4){ // G
                SHELLS[i] = 8;
            }
        }
        else if(N == 2){
            if (L==0){ //S
                SHELLS[i] = 6;
            }
            else if ( L==1){ // P
                SHELLS[i] = 7;
            }
            else if (L==2){ // D
                SHELLS[i] = 8;
            }
        }
        else if(N == 3){
            SHELLS[i] = 8; // S
        }
        else{
            QN_1body_jj Error_Shell = {N,L,J};
            std::cout << " (ERROR) Shell "<< index_antoine(Error_Shell) << " cannot be identified." << std::endl;
            SHELLS[i] = 1e5;
        }

        //   Check the lower level,
        if(SHELLS[i] < Lower_Shell){
            Lower_Shell =  SHELLS[i];
        }

    }

    return Lower_Shell;
}

int core_shell(QN_1body_jj *VALENCE, int Number_of_distint_elements){
    /// Algorithm to set the full filled shells bellow the valence space (VALENCE)
    if(Number_of_distint_elements == 0){
        std::cout << "(ERROR) There is no Valence space, restart the program"<< std::endl;
        return 0;
    }
    /// find the maximum and minimum "n" (n >= 0)
    int Max_n_Valence = 0;
    int Min_n_Valence = 1e5; // arbitrary large number
    for(int i= 0; i<Number_of_distint_elements; i++){
        if( VALENCE[i].n < Min_n_Valence){
            Min_n_Valence = VALENCE[i].n;
        }
        if( VALENCE[i].n > Max_n_Valence){
            Max_n_Valence = VALENCE[i].n;
        }
    }
    /// find the maximum and minimum "l" (l >= 0)
    int Max_l_Valence = 0;
    int Min_l_Valence = 1e5; // arbitrary large number
    int i_min = 0;           // index of the low orbital state
    for(int i= 0; i<Number_of_distint_elements; i++){
        if( VALENCE[i].l < Min_l_Valence){
            Min_l_Valence = VALENCE[i].l;
            i_min = i;
        }
        if( VALENCE[i].l > Max_l_Valence){
            Max_l_Valence = VALENCE[i].l;
        }
    }
    std::cout << "Min n and l valence:" << VALENCE[i_min].n << "  "<< Min_l_Valence << std::endl;
    std::cout << "Max n and l valence:" << Max_n_Valence << "  "<< Max_l_Valence << std::endl;
    /// Selection of the bounds of the valence space (and then of the core shells)
    /// Shells are numerated from 1 (S) to 7(PFHI)
    int Max_Shell_Valence = 0;
    int Min_Shell_Valence = 0;
    int Core_Shell = 0;

    /// Set the maximum valence space shell //////////////////////////////
    switch (Max_n_Valence)
    {
        // If n (max) = 0 there is 3 possibilities of 1 l orbital momentum
        case 0:
        {
            switch(Max_l_Valence){
                case 0:{Max_Shell_Valence = 1;}  // S shell
                break;
                case 1:{Max_Shell_Valence = 2;}  // P shell
                break;
                case 3:{Max_Shell_Valence = 4;}  // F shell
                break;
                default:{std::cout << " (ERROR 0) Not a Shell Model Maximum Shell"<< std::endl;}
                }
                break;
        }
        break;
        // For n max = 1 or 2 case, the lower l orbital momentums of the shells switch parities,
        // this make easier to identify by the lower l = S(+) or P(-)
        case 1:
        {
            switch(Min_l_Valence){
                case 0:{Max_Shell_Valence = 3;}  // SD shell
                break;
                case 1:{Max_Shell_Valence = 5;}  // PF(G) shell
                break;
                default:{std::cout << " (ERROR 1) Not a Shell Model Maximum Shell"<< std::endl;}
                }
                break;
        }
        break;
        case 2:
        {
            switch(Min_l_Valence){
                case 0:{Max_Shell_Valence = 6;}  // SDG(H) shell
                break;
                case 1:{Max_Shell_Valence = 7;}  // PFH(I) shell
                break;
                default:{std::cout << " (ERROR 2) Not a Shell Model Maximum Shell"<< std::endl;}
                }
                break;
        }
        break;
        case 3:
        {
            Max_Shell_Valence = 8;  // SDGI(J) shell
        }
        break;
        default:{std::cout << " (ERROR 3) Not a Shell Model Maximum Shell"<< std::endl;}
    }


    /// Now the minimum valence space shell //////////////////////////////
    switch (Min_l_Valence)
    {
        case 0:
        {
            if(VALENCE[i_min].n == 3){
                Min_Shell_Valence = 8; //   SDGI(J) Shell
            }
            if(VALENCE[i_min].n == 2){
                Min_Shell_Valence = 6; //   SDG(H) Shell
            }
            if(VALENCE[i_min].n == 1){
                Min_Shell_Valence = 3; //   SD Shell
            }
            if(VALENCE[i_min].n == 0){
                Min_Shell_Valence = 1; //   S Shell
            }

        }
        break;
        case 1:
        {
            if(VALENCE[i_min].n == 2){
                Min_Shell_Valence = 7; //   PFH(I) Shell
            }
            if(VALENCE[i_min].n == 1){
                Min_Shell_Valence = 5; //   PF(G) Shell
            }
            if(VALENCE[i_min].n == 0){
                Min_Shell_Valence = 2; //   P Shell
            }
        }
        break;
        case 3:
        {
            Min_Shell_Valence = 4; // F Shell
        }
        break;
        default:
        {
            std::cout << " (ERROR 0) Not a Shell Model Minimum Shell"<< std::endl;
            Min_Shell_Valence = Max_Shell_Valence;
        }
    }

    /*
    switch (Min_n_Valence)
    {
        // If n (max) = 0 there is 3 possibilities of 1 l orbital momentum
        case 0:
        {
            switch(Max_l_Valence){
                case 0:{Min_Shell_Valence = 1;}  // S shell
                break;
                case 1:{Min_Shell_Valence = 2;}  // P shell
                break;
                case 3:{Min_Shell_Valence = 4;}  // F shell
                break;
                default:
                    {
                        std::cout << " (ERROR 0) Not a Shell Model Minimum Shell"<< std::endl;
                        Min_Shell_Valence = Max_Shell_Valence - 1;
                    }
                }
                break;
        }
        break;
        // For n max = 1 or 2 case, the lower l orbital momentums of the shells switch parities,
        // this make easier to identify by the lower l = S(+) or P(-)
        case 1:
        {
            switch(Min_l_Valence){
                case 0:{Min_Shell_Valence = 3;}  // SD shell
                break;
                case 1:{Min_Shell_Valence = 5;}  // PF(G) shell
                break;
                default:
                    {
                        std::cout << " (ERROR 1) Not a Shell Model Minimum Shell"<< std::endl;
                        Min_Shell_Valence = Max_Shell_Valence - 1;
                    }
                }
                break;
        }
        break;
        case 2:
        {
            switch(Min_l_Valence){
                case 0:{Min_Shell_Valence = 6;}  // SDG(H) shell
                break;
                case 1:{Min_Shell_Valence = 7;}  // PFH(I) shell
                break;
                default:
                    {
                        std::cout << " (ERROR 2) Not a Shell Model Minimum Shell"<< std::endl;
                        Min_Shell_Valence = Max_Shell_Valence - 1;
                    }
                }
                break;
        }
        break;
        case 3:
        {
            Min_Shell_Valence = 8;  // SDGI(J) shell
        }
        break;
        default:
            {
                std::cout << " (ERROR 3) Not a Shell Model Minimum Shell"<< std::endl;
                Min_Shell_Valence = Max_Shell_Valence - 1;
            }
    }
    */

    if (Min_Shell_Valence > 0){
        Core_Shell = Min_Shell_Valence - 1;
    }
    else
    {
        Min_Shell_Valence = 0;
    }

    return Core_Shell;
}

int core_number(int core_shell){
    /// Return the total number of particles which fills the core of the nucleus
    /// (magic numbers)
    switch (core_shell){
        case 1:
        {
           return 2;
        }
        break;
        case 2:
        {
            return 8;
        }
        break;
        case 3:
        {
            return 20;
        }
        break;
        case 4:
        {
            return 28;
        }
        break;
        case 5:
        {
            return 50;
        }
        break;
        case 6:
        {
            return 82;
        }
        break;
        case 7:
        {
            return 126;
        }
        break;
        default:
        {
            return 0;
        }
    }
}


void antoine_output(int nothing){
    /// Antoine's feeder of matrix elements, function to choose different shells to evaluate 
    /// all matrix elements possible and give in the appropriate way the input to Antoine.

    bool repeat;
    do{ //repetition loop

        /// Data file making (different file type)
        std::ofstream f_salida;

        std::string M_E_seletion;
        std::cout << " Select a type of interaction:" <<std::endl;
        std::cout << " -> Surface Delta Interaction                  (type SDI)" <<std::endl;
        std::cout << " -> Modified Surface Delta Interaction         (type MSDI)" <<std::endl;
        std::cout << " -> Multipolar Decomposition of a General Central" <<std::endl;
        std::cout << "               Exchange Interaction            (type MULTI)" <<std::endl;
        std::cout << " -> Central Interaction                        (type CENTRAL)" <<std::endl;
        std::cout << " -> Spin Orbit Interaction                     (type SO)" <<std::endl;
        std::cout << " -> Tensor Interaction                         (type TENSOR)" <<std::endl;
        std::cout << " -> Brink Boeker Interaction                   (type BB)" <<std::endl;
        std::cin >> M_E_seletion;
        std::cout << std::endl;

        bool spin_exchange;
	int option_BB;
        int order;
        std::string L_shell_name, filename, ss;

        //case "SDI" :
        if(M_E_seletion == "SDI")
            {
                filename = "data/sdi_";

                std::cout<< "Isospin Exchange  (type 1 for true or 0 for false) ? ";
                std::cin >> spin_exchange;

                if(spin_exchange)   {filename = filename + "ie1_";}
                else                {filename = filename + "ie0_";}

                std::cout << " Name the different Orbital l shells to use (in capital letters and in one word. f.e: SPD...): " ;
                std::cin >> L_shell_name;
                filename = filename + L_shell_name;

                f_salida.open(filename.c_str(),std::fstream::out);
                f_salida<<"  SDI INT. (";
                f_salida << L_shell_name << " SHELL)"<<std::endl;

            }
        //case "MSDI" :
        else if(M_E_seletion == "MSDI")
            {
                filename = "data/msdi_";

                std::cout << " Name the different Orbital l shells to use (in capital letters and in one word f.e: SPD...): " ;
                std::cin >> L_shell_name;
                filename = filename + L_shell_name;

                f_salida.open(filename.c_str(),std::fstream::out);
                f_salida<<"  MSDI INT. (";
                f_salida << L_shell_name << " SHELL)"<<std::endl;
            }
        //case "MULTI" :
        else if(M_E_seletion == "MULTI")
            {
                filename = "data/multi_";

                std::cout << " Name  the different Orbital l shells to use (in capital letters and in one word f.e: SPD...): " ;
                std::cin >> L_shell_name;

                std::cout<< "Set maximum multipole order for the decomposition: ";
                std::cin >> order;

                std::ostringstream convert;
                convert << order;
                ss = convert.str();

                filename = filename + L_shell_name + ss;
                f_salida.open(filename.c_str(),std::fstream::out);
		f_salida<<"  MULTIPOLAR DECOMPOSITION. (";
                f_salida << L_shell_name << " SHELL) (Order "<< order <<")"<<std::endl;
            }
        //case "CENTRAL" :
        else if(M_E_seletion == "CENTRAL")
            {
                filename = "data/central_";

                std::cout << " Name the different Orbital l shells to use (in capital letters and in one word f.e: SPD...): " ;
                std::cin >> L_shell_name;
                filename = filename + L_shell_name;

                f_salida.open(filename.c_str(),std::fstream::out);
                f_salida<<"  CENTRAL INT. (";
                f_salida << L_shell_name << " SHELL)"<<std::endl;
            }
        //case "Spin Orbit" :
        else if(M_E_seletion == "SO")
            {
                filename = "data/spin_orbit_";

                std::cout << " Name the different Orbital l shells to use (in capital letters and in one word f.e: SPD...): " ;
                std::cin >> L_shell_name;
                filename = filename + L_shell_name;

                f_salida.open(filename.c_str(),std::fstream::out);
                f_salida<<"  SO INT. (";
                f_salida << L_shell_name << " SHELL)"<<std::endl;
            }
        //case "Tensor"
        else if(M_E_seletion == "TENSOR")
            {
                filename = "data/tensor_";

                std::cout << " Name the different Orbital l shells to use (in capital letters and in one word f.e: SPD...): " ;
                std::cin >> L_shell_name;
                filename = filename + L_shell_name;

                f_salida.open(filename.c_str(),std::fstream::out);
                f_salida<<"  TENSOR INT. (";
                f_salida << L_shell_name << " SHELL)"<<std::endl;
            }
        //case "BB"
        else if(M_E_seletion == "BB")
            {
                filename = "data/bb_";

                std::cout << " Name the different Orbital l shells to use (in capital letters and in one word f.e: SPD...): " ;
                std::cin >> L_shell_name;
                filename = filename + L_shell_name;

                f_salida.open(filename.c_str(),std::fstream::out);
                f_salida<<"  BB INT. (";
                f_salida << L_shell_name << " SHELL)"<<std::endl;
		
		std::cout << " Option to calculate BB :"<<std::endl;
		std::cout << " 	 0  ->  Exchange operators; full recoupling"<<std::endl;
		std::cout << " 	 1  ->  Exchange operators; recoupling with LSjj coefficients"<<std::endl;
		std::cout << " 	 2  ->  Spin & Isospin dependence (full recoupling)"<<std::endl;
		std::cout << " 	 3  ->  Multipolar descomposition with GL Quadratures"<<std::endl;
		std::cout << " 	 ...  ";
		std::cin >> option_BB;
		while((option_BB!=0) && (option_BB!=1) && (option_BB!=2) && (option_BB!=3)){
		    std::cout << " No!, select 0,1,2 or 3 :  ";
		    std::cin >> option_BB;
		}
            }


        /// Shell specification (numbers) and proton-neutron distinction
        int Number_of_distint_elements;
        bool charge_symetry = true;

        std::cout << "Set number of different Quantum Numbers : ";
        std::cin >> Number_of_distint_elements;
        std::cout<< std::endl;
        /*if(L_shell_name.length() > Number_of_distint_elements)
        {
            std::cout << " (Error), More distinct elements needed. The program have finished, try again"<<std::endl;
            return;
        }*/

        f_salida<<" "<< charge_symetry <<" "<<Number_of_distint_elements;

            // Specification of manual or default input and the election of the Harmonic oscillator energy parameter.
        bool sp_energy_manual_input = false;   // default
        bool h_bar_omega_maunal_input = false; // default
        double h_bar_omega_input = 0.0;
	
        std::cout << " Input of Single particle energies, (type 1 to do manually or 0 by default) : ";
        std::cin >> sp_energy_manual_input;
	while((sp_energy_manual_input!=0) && (sp_energy_manual_input!=1) ){
	    std::cout<< std::endl;
	    std::cout << " No!, select 0 or 1 :  ";
	    std::cin >> sp_energy_manual_input;
	    std::cout<< std::endl;
	}
        std::cout<< std::endl;
	if(!sp_energy_manual_input){
	    std::cout << " Input of Harmonic oscillator Energy ''h_bar omega'', (type 1 to do manually or 0 by default) : ";
	    std::cin >> h_bar_omega_maunal_input;
	    while((h_bar_omega_maunal_input!=0) && (h_bar_omega_maunal_input!=1) ){
		std::cout<< std::endl;
		std::cout << " No!, select 0 or 1 :  ";
		std::cin >> h_bar_omega_maunal_input;
		std::cout<< std::endl;
	    }
	    
	    if(h_bar_omega_maunal_input){
		std::cout << " 	>> Then, Introduce the ''h_bar omega'' value (MeV): ";
		std::cin >> h_bar_omega_input;
	    	std::cout<< std::endl;
	    }
	    
	    
	}

        QN_1body_jj *SET;
        SET = new QN_1body_jj [Number_of_distint_elements];
        double *Single_particle_energies;
        Single_particle_energies = new double[Number_of_distint_elements];

        int n_input, l_input, j_num;
        double sp_energy;
	
	std::cout << "-----------------------------------------------------------------------------------" << std::endl;
	std::cout << " Start the Setting of Quantum Orbitals for the Matrix Elements ..." << std::endl;
	std::cout << std::endl;
        for (int i = 0; i < Number_of_distint_elements; i++){
	    
	    std::cout << " QN (" << i+1 <<"/"<< Number_of_distint_elements << "):" << std::endl;
            std::cout << " Set n radial number (0,1, ...): " ;
            std::cin >> n_input ;
            //std::cout<< std::endl;
            std::cout << " Set l Orbital Number (0 = s, 1 = p, ...): " ;
            std::cin >> l_input;
            //std::cout<< std::endl;
            std::cout << " Set the numerator (odd value) for the total angular momentum j = ?/2: " ;
            std::cin >> j_num ;
            while(j_num % 2 == 0){
                std::cout << " Even value inserted, put an odd number j = ?/2: " ;
                std::cin >> j_num ;
            }
            std::cout << " j = "<< Fraction (j_num,2);
            
            SET[i] = {n_input, l_input, Fraction (j_num,2)};
	    std::cout<< "       Antoine index: "<< name_index_antoine(SET[i]) << std::endl;
	    std::cout<< std::endl;
	    
            // Single particle energies manual or default
            if(sp_energy_manual_input){
                std::cout << " Set single particle energy (MeV): " ;
                std::cin >> sp_energy;
                std::cout<< std::endl;
                Single_particle_energies[i] = sp_energy;
            }
            else{
		if(h_bar_omega_maunal_input){
		    Single_particle_energies[i] = h_bar_omega_input* (2*SET[i].n + SET[i].l + (3./2));
		}
                else{
		    Single_particle_energies[i] = h_bar_omega(A)* (2*SET[i].n + SET[i].l + (3./2));
		}
            }

            f_salida << " " << name_index_antoine(SET[i]);

        }
        f_salida<<std::endl;
        std::cout << "Passing of Quantum Numbers Completed" << std::endl;

        /// print of the single particle energy
        f_salida<<"   ";
        for(int i = 0; i < Number_of_distint_elements; i++){
            f_salida << Single_particle_energies[i] << "  ";
        }
        f_salida<<std::endl;

        /// shell fill, for density dependent interactions
        /// fix the last number to 0 to not apply:
        bool default_core = true;
        bool density_effect = false;
        int density_just_for_2B_ME = 0;
        int core_number_of_particles = 0;

        int * SHELLS;
        SHELLS = new int[Number_of_distint_elements];
	
	std::cout << std::endl;
        std::cout << "Core and density effects setting:" << std::endl;
        std::cout << "-----------------------------------------------------------------------------------" << std::endl;
        /* std::cout << "  Set DEFAULT core for the input shells (same for protons and neutrons)?" << std::endl;
        std::cout << "(type 1 by default or 0 to do manually) : ";
        std::cin >> default_core;*/
        std::cout << "  A default inert core of " << core_number(shell_identification ( SHELLS, SET, Number_of_distint_elements) - 1)
                << " protons and " <<  core_number(shell_identification ( SHELLS, SET, Number_of_distint_elements) - 1)
                <<" neutrons  for the input shells is settled." << std::endl;
        std::cout << std::endl;
        std::cout << "  Apply density effect:"<< std::endl;
        std::cout<< "       * Don't apply                                                     (0)"<< std::endl;
        std::cout<< "       * Apply just for two body matrix elements                         (1)"<< std::endl;
        std::cout <<"       * Apply for single particle energies and two body matrix elements (2)  : ";
        std::cin >> density_just_for_2B_ME;

        if(default_core){
            //core_number_of_particles = core_number(core_shell(SET,Number_of_distint_elements));
            core_number_of_particles = core_number(shell_identification ( SHELLS, SET, Number_of_distint_elements) - 1);
            // presence of filling __ int __ int __ double (parameters) __ apply density dependence (double)

            if(density_just_for_2B_ME == 1){
                f_salida << "    1    "<< core_number_of_particles <<"    "<< core_number_of_particles;
                f_salida << "    0.333333    0.000000" <<std::endl;
            }
            else if (density_just_for_2B_ME == 2){
                f_salida << "    2    "<< core_number_of_particles <<"    "<< core_number_of_particles;
                f_salida << "    0.333333    0.000000" <<std::endl;
            }
            else{
                f_salida << "    0    "<< core_number_of_particles <<"    "<< core_number_of_particles;
                f_salida << "    0.333333    0.000000" <<std::endl;
            }

        }
        else{
            if(density_just_for_2B_ME == 1){
                f_salida << "    1    ";
            }
            else if(density_just_for_2B_ME == 2){
                f_salida << "    2    ";
            }
            else{
                f_salida << "    0    ";
            }
            // Manual input
            std::cout<< "   >> Enter core number of particles for protons:  ";
            std::cin >> core_number_of_particles;
            f_salida <<  core_number_of_particles;
            std::cout<< "   >> Enter core number of particles for neutrons: ";
            std::cin >> core_number_of_particles;
            f_salida << "    " << core_number_of_particles << "    0.333333    0.000000" <<std::endl;
        }
        
        
	// Almost all the interactions employ the oscillator length b, add A in order to evalate it. 
	double b_param = b_lenght(A);// Default
	bool bool_b_length_manual = false;
	double r_square_m_e;
	bool add_r_square = false;
	
	if((M_E_seletion != "MSDI") || (M_E_seletion == "SDI")){
	    std::cout << std::endl;
	    std::cout << " Oscillator length b is required: " << std::endl;
	    std::cout << "    Type (0) to set with A mass number by default or (1) to introduce directly the b parameter: ";
	    std::cin  >> bool_b_length_manual;
	    if(!bool_b_length_manual){
		std::cout << " 	  >> Insert mass number of the nucleus (A):  ";
		std::cin  >> A;
		b_param = b_lenght(A);
		std::cout << "	 b(" << A <<") = " << b_param <<" fm"<<std::endl;
	    }
	    else{
		std::cout << " 	  >> Then, Set oscillator length in fm:   ";
		std::cin >>b_param;
		
		std::cout << " 	  >> Set mass number A:             ";
		std::cin >>A;
	    }
	    
	}
	

        /// Perform of all combinations of two particles (with repetition)
        QN_1body_jj a,b,c,d;
        QN_2body_jj_Coupling Bra, Ket;
        Fraction j_bra_min, j_ket_min,j_bra_max, j_ket_max, J;
        //Fraction J_min, J_max;

        int i_d_min;
        double numeric;


        /// Selection of "Number_of_distint_elements" states among 2 2-body states, without repetition
        /** COMMENT
         * state |ab> = phase* |ba> and (by a phase) < ab V cd> = < ba cd> = < ab dc> = < ba dc> = < cd ab> = ...
         *
         *    First, we establish the order to proceed (arbitrarily, the SET order)
         * so, Quantum Numbers of the first are always equal or great than the second
         *      a) Condition:( i_b >= i_a && i_d >= i_c)
         * to not repeat the numbers by changing the bra by the ket, bra elements precede the ket ones
         *      b) Condition:(i_c >= i_a)
         *
         *      Finally, when i_a =  i_c, the "i_d = i_c to NDE" definition doubles several evaluated values
         * f.e (ia ib ic id): 0001 0100 or 1113 1311. So here again, if we choose (i_d >= i_b) we omit the doubled
         * to respect Condition a), there is:
         *      c) Condition: (i_d = max(i_b,i_c) if (i_a == i_c))
         */
        for(int i_a=0; i_a < Number_of_distint_elements; i_a++){

            a = SET[i_a];
            for(int i_b=i_a; i_b < Number_of_distint_elements; i_b++){

                b = SET[i_b];
                Bra = {a.n,a.l,a.j, b.n,b.l,b.j, 0, 0, 0, 0};
                // default J,T

                // momentum addition limits
                j_bra_min = abs(a.j - b.j);
                j_bra_max = a.j + b.j;

                for(int i_c=i_a; i_c < Number_of_distint_elements; i_c++){

                    c = SET[i_c];
                    // Selection of the last index.
                    if(i_a == i_c){
                        i_d_min = max(i_c, i_b);
                    }
                    else{
                        i_d_min = i_c;
                    }
                    for(int i_d = i_d_min; i_d < Number_of_distint_elements; i_d++){
                        //std::cout << i_a << i_b << " " << i_c << i_d << std::endl;
                        d = SET[i_d];
                        Ket = {c.n,c.l,c.j, d.n,d.l,d.j, 0,0,0,0};

                        // ket moment addition limits
                        j_ket_min = abs(c.j - d.j);
                        j_ket_max = c.j + d.j;

                        if((j_ket_min <= j_bra_max) && (j_bra_min <= j_ket_max)){

                            // Print the numbers
                            f_salida << " 0 1 " << name_index_antoine(a) << " " << name_index_antoine(b) << " "
                                                << name_index_antoine(c) << " " << name_index_antoine(d) << " ";
                            f_salida << max(j_bra_min,j_ket_min)<< " " << min(j_bra_max,j_ket_max);
                            f_salida << std::endl;
                            /*
                            // Print the numbers
                            f_salida << " 0 1 " << index_antoine(a) << " " << index_antoine(b) << " "
                                                << index_antoine(c) << " " << index_antoine(d) << " ";
                            f_salida << max(j_bra_min,j_ket_min)<< " " << min(j_bra_max,j_ket_max);
                            f_salida << std::endl;

*/
                            for(int T = 0; T <= 1; T++){

                                for(J = max(j_bra_min,j_ket_min); J <= min(j_bra_max,j_ket_max); J += 1){
                                    Bra.J = J;
                                    Ket.J = J;
                                    Bra.T = T;
                                    Ket.T = T;
				      
				    //introduce the r_square m.e
				    if(add_r_square){
					r_square_m_e = -0.5 * h_bar_omega(A)* M_E_Brink_Boeker(Bra,Ket,4,b_param);
				    }

                                    /// Do the selection
				    
                                    //switch(M_E_seletion){
                                    //    case "SDI" :
                                    if(M_E_seletion == "SDI")
                                    {
                                        // Print the matrix elements
                                        numeric = M_E_SDI(Bra, Ket, spin_exchange);// bool isospin_exchange
                                        if(numeric < 0){
                                            if (fabs(numeric) < 1e-7){
                                                // it is zero, display as positive
                                                f_salida << "   0.00000" ;
                                            }
                                            else{
                                                f_salida << "  " << std::fixed << std::setprecision(5) << numeric;
                                            }

                                        }
                                        else{
                                            if (fabs(numeric) < 1e-7){
                                                // it is zero, display as positive
                                                f_salida << "   0.00000" ;
                                            }
                                            else{
                                                f_salida << "   " << std::fixed << std::setprecision(5) << numeric;
                                            }
                                        }
                                    }

                                    //case "MSDI" :
                                    else if(M_E_seletion == "MSDI")
                                    {
                                        // Print the matrix elements
                                        numeric = M_E_MSDI(Bra, Ket);// bool isospin_exchange
                                        if(numeric < 0){
                                            if (fabs(numeric) < 1e-7){
                                                // it is zero, display as positive
                                                f_salida << "   0.00000" ;
                                            }
                                            else{
                                                f_salida << "  " << std::fixed << std::setprecision(5) << numeric;
                                            }

                                        }
                                        else{
                                            if (fabs(numeric) < 1e-7){
                                                // it is zero, display as positive
                                                f_salida << "   0.00000" ;
                                            }
                                            else{
                                                f_salida << "   " << std::fixed << std::setprecision(5) << numeric;
                                            }

                                        }
                                    }

                                    // case "MULTI" :
                                    else if(M_E_seletion == "MULTI")
                                    {
                                        // Print the matrix elements
                                        numeric = M_E_Spin_Isospin_Multipolar(Bra, Ket, order,b_param);// bool isospin_exchange
                                        if(numeric < 0){
                                            if (fabs(numeric) < 1e-7){
                                                // it is zero, display as positive
                                                f_salida << "   0.00000" ;
                                                //f_salida << "   " << std::fixed << std::setprecision(5) << abs(numeric);
                                            }
                                            else{
                                                f_salida << "  " << std::fixed << std::setprecision(5) << numeric;
                                            }

                                        }
                                        else{
                                            if (fabs(numeric) < 1e-7){
                                                // it is zero, display as positive
                                                //f_salida << "   " << std::fixed << std::setprecision(5) << abs(numeric);
                                                f_salida << "   0.00000" ;
                                            }
                                            else{
                                                f_salida << "   " << std::fixed << std::setprecision(5) << numeric;
                                            }

                                        }
                                    }

                                    // Central Coupling
                                    else if (M_E_seletion == "CENTRAL")
                                    {
                                        // Print the matrix elements
                                        numeric = M_E_Central_jj(Bra, Ket, 0, b_param);//
                                        if(numeric < 0){
                                            if (fabs(numeric) < 1e-7){
                                                // it is zero, display as positive
                                                f_salida << "   0.00000" ;
                                                //f_salida << "   " << std::fixed << std::setprecision(5) << abs(numeric);
                                            }
                                            else{
                                                f_salida << "  " << std::fixed << std::setprecision(5) << numeric;
                                            }

                                        }
                                        else{
                                            if (fabs(numeric) < 1e-7){
                                                // it is zero, display as positive
                                                //f_salida << "   " << std::fixed << std::setprecision(5) << abs(numeric);
                                                f_salida << "   0.00000" ;
                                            }
                                            else{
                                                f_salida << "   " << std::fixed << std::setprecision(5) << numeric;
                                            }

                                        }
                                    }
                                    // Spin Orbit Coupling
                                    else if (M_E_seletion == "SO")
                                    {
                                        // Print the matrix elements
                                        numeric = M_E_Spin_Orbit(Bra, Ket, 0, b_param);//
                                        if(numeric < 0){
                                            if (fabs(numeric) < 1e-7){
                                                // it is zero, display as positive
                                                f_salida << "   0.00000" ;
                                                //f_salida << "   " << std::fixed << std::setprecision(5) << abs(numeric);
                                            }
                                            else{
                                                f_salida << "  " << std::fixed << std::setprecision(5) << numeric;
                                            }

                                        }
                                        else{
                                            if (fabs(numeric) < 1e-7){
                                                // it is zero, display as positive
                                                //f_salida << "   " << std::fixed << std::setprecision(5) << abs(numeric);
                                                f_salida << "   0.00000" ;
                                            }
                                            else{
                                                f_salida << "   " << std::fixed << std::setprecision(5) << numeric;
                                            }

                                        }
                                    }

                                    // Tensor Coupling
                                    else if(M_E_seletion == "TENSOR")
                                    {
                                        // Print the matrix elements
                                        numeric = M_E_Tensor(Bra, Ket, 0, b_param);//
                                        if(numeric < 0){
                                            if (fabs(numeric) < 1e-7){
                                                // it is zero, display as positive
                                                f_salida << "   0.00000" ;
                                                //f_salida << "   " << std::fixed << std::setprecision(5) << abs(numeric);
                                            }
                                            else{
                                                f_salida << "  " << std::fixed << std::setprecision(5) << numeric;
                                            }

                                        }
                                        else{
                                            if (fabs(numeric) < 1e-7){
                                                // it is zero, display as positive
                                                //f_salida << "   " << std::fixed << std::setprecision(5) << abs(numeric);
                                                f_salida << "   0.00000" ;
                                            }
                                            else{
                                                f_salida << "   " << std::fixed << std::setprecision(5) << numeric;
                                            }

                                        }
                                    }

                                    //case "BB" :
                                    else if(M_E_seletion == "BB")
                                    {
                                        // Print the matrix elements
					std::cout << name_index_antoine(Bra) <<" | "<< name_index_antoine(Ket) <<" JT:"<<J<<" "<<T<<std::endl;
                                        numeric = M_E_Brink_Boeker(Bra, Ket, option_BB, b_param);// bool isospin_exchange
					
					if(add_r_square){numeric += r_square_m_e;}
					
                                        if(numeric < 0){
                                            if (fabs(numeric) < 1e-7){
                                                // it is zero, display as positive
                                                f_salida << "   0.00000" ;
                                            }
                                            else{
                                                f_salida << "  " << std::fixed << std::setprecision(5) << numeric;
                                            }

                                        }
                                        else{
                                            if (fabs(numeric) < 1e-7){
                                                // it is zero, display as positive
                                                f_salida << "   0.00000" ;
                                            }
                                            else{
                                                f_salida << "   " << std::fixed << std::setprecision(5) << numeric;
                                            }

                                        }
                                    }

                                    // more cases
                                }
                                f_salida << std::endl;
                            }
                        }
                        else{} // Impossible addition

                    }
                }
            }
        }

        delete [] SET;
        delete [] SHELLS;

        f_salida.close();


        // Repetition Loop
        repeat = false;
	std::cout << std::endl;
        std::cout << "-----------------------------------------------------------------------" << std::endl;
        std::cout << " The Matrix Elements have been done and file is written," << std::endl;
        std::cout << " Do you want to compute another shells or interaction ? (type 1 for yes) ";
        std::cin  >> repeat;
        std::cout << std::endl;
    }while(repeat);

}
