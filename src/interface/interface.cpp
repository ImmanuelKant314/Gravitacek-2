#include "interface/interface.hpp"
#include "interface/usefullfunctions.hpp"
#include "gravitacek2/geomotion/spacetimes.hpp"
#include "gravitacek2/integrator/integrator.hpp"
#include "gravitacek2/integrator/odesystems.hpp"
#include "gravitacek2/chaos/linearized_evolution.hpp"

#include <stdexcept>
#include <iostream>
#include <algorithm>
#include <fstream>
#include <array>
#include <cmath>

#include <gsl/gsl_linalg.h>
#include <gsl/gsl_eigen.h>
#include <gsl/gsl_complex.h>
#include <gsl/gsl_complex_math.h>

std::string Interface::substitute(std::string text)
{
    for (int i = 0; i<macros.size(); i++)
    {
        if (text.empty())
            return text;

        std::string macro = macros[i];
        std::string value = values[i];

        size_t start_pos = 0;
        while ((start_pos = text.find(macro, start_pos)) != std::string::npos) 
        {
            text = text.replace(start_pos, macro.length(), value);
        }
    }
    return text;
}

bool Interface::macro_name_valid(std::string name)
{
    // check first character
    char c = name[0];
    if (!('A' <= c && 'Z' >= c))
        return false;
    for (auto& c : name)
    {
        if(!(('A' <= c && 'Z' >=c) || (c == '_') || ('0' <= c || '9' >= c)))
            return false;
    }
    return true;
}

std::string Interface::strip(std::string text)
{
    int start, end;
    int i;
    for (i = 0; i < text.length(); i++)
    {
        if (text[i] != ' ')
        {
            start = i;
            break;
        }
    }

    if (i == text.length())
        return "";

    for (i = text.length()-1; i>=0; i--)
    {
        if (text[i] != ' ')
        {
            end = i;
            break;
        }
    }
    return text.substr(start, end-start+1);
}

void Interface::find_command_name(std::string text, std::string &command, std::string &rest)
{
    for (int i = 0; i < text.length(); i++)
    {
        if(text[i] == ' ' || text[i] == '(')
        {
            command = text.substr(0, i);
            rest = text.substr(i);
            return;
        }
    }
    command = text;
    rest = "";
}

bool Interface::try_apply_operators(std::string text)
{
    std::string command, rest;
    this->find_command_name(strip(text), command, rest);
    if (command == "def")
    {
        this->define_macro(rest);
        return true;
    }
    else if(command == "delete")
    {
        this->delete_macro(rest);
        return true;
    }
    else if(command == "print")
    {
        this->print_macro(rest);
        return true;
    }
    else if(strip(text) == "printm")
    {
        this->print_all_macros();
        return true;
    }
    else if(command == "help")
    {
        this->help(rest);
        return true;
    }

    return false;
}

void Interface::define_macro(std::string text)
{
    text = strip(text);
    // extract name
    int i;
    for (i = 0; i < text.length(); i++)
    {
        if (text[i] == ' ')
        {
            std::string name = text.substr(0, i);
            // check validity on the name
            if (!macro_name_valid(name))
                throw std::invalid_argument("invalid macro name");
            std::string value = strip(text.substr(i+1));
            value = substitute(value);
            for (int j = 0; j < macros.size(); j++)
            {
                if(macros[j] == name)
                {
                    values[j] = value;
                    return;
                }
            }
            macros.push_back(name);
            values.push_back(value);
            return;
        }
    }
    throw std::invalid_argument("macro is empty");
}

void Interface::delete_macro(std::string text)
{
    text = strip(text);
    for (int i = 0; i < macros.size(); i++)
    {
        if(macros[i] == text)
        {
            macros.erase(macros.begin()+i);
            values.erase(values.begin()+i);
            return;
        }
    }
    throw std::invalid_argument("macro with this name can not be deleted");

    // if found, delete it (write macro deleted)
}

void Interface::print_macro(std::string text)
{
    text = strip(text);
    for (int i = 0; i < macros.size(); i++)
    {
        if (macros[i] == text)
        {
            std::cout << values[i] << std::endl;
            return;
        }
    }
    throw std::invalid_argument("no macro with this name");
}

void Interface::print_all_macros()
{
    for (int i = 0; i < macros.size(); i++)
    {
        std::cout << macros[i] << ":" << values[i] << std::endl;
    }
}

void Interface::help(std::string text)
{
    text = this->strip(text);
    if (text.length() == 0)
    {
        int n = this->help_name.size();
        if (n > 0)
            std::cout << this->help_name[0];
        for (int i = 1; i<n; i++)
            std::cout << ", " << this->help_name[i];
        std::cout << std::endl;
    }
    else
    {
        int i;
        for (i = 0; i < this->help_name.size(); i++)
        {
            if (text == this->help_name[i])
            {
                std::cout << this->help_text[i] << std::endl;
                break;
            }
        }
        if (i == this->help_name.size())
            throw std::runtime_error("help for " + text + " was not found");
    }
}

std::shared_ptr<gr2::Weyl> Interface::create_weyl_spacetime(std::string text)
{
    text = strip(text);
    std::string spacetime_name, args_text;
    this->find_command_name(text, spacetime_name, args_text);
    auto args = this->find_function_arguments(args_text);
    std::shared_ptr<gr2::Weyl> spacetime;

    if (spacetime_name == "CombinedWeyl")
    {
        std::vector<std::shared_ptr<gr2::Weyl>> sources = {};
        for (auto &arg : args)
        {
            sources.push_back(this->create_weyl_spacetime(arg));
        }
        spacetime = std::make_shared<gr2::CombinedWeyl>(sources);
    }
    else if (spacetime_name == "WeylSchwarzschild")
    {
        if (args.size() != 1)
            throw std::invalid_argument("invalid number of arguments for WeylSchwarzschild");
        spacetime = std::make_shared<gr2::WeylSchwarzschild>(std::stold(args[0]));
    }
    else if (spacetime_name == "BachWeylRing")
    {   
        if (args.size() != 2)
            throw std::invalid_argument("invalid number of arguments for BachWeylRing");
        spacetime = std::make_shared<gr2::BachWeylRing>(std::stold(args[0]), std::stold(args[1]));
    }
    else if (spacetime_name == "InvertedKuzminToomreDisk")
    {
        if (args.size() != 3)
            throw std::invalid_argument("invalid number of arguments for InvertedKuzminToomreDisk");
        spacetime = std::make_shared<gr2::InvertedKuzminToomreDisk>(std::stoi(args[0]), std::stold(args[1]), std::stold(args[2]));
    }
    else if (spacetime_name == "InvertedMorganMorganDisk")
    {
        if (args.size() != 3)
            throw std::invalid_argument("invalid number of arguments for InvertedMorganMorganDisk");
        spacetime = std::make_shared<gr2::InvertedMorganMorganDisk>(std::stoi(args[0]), std::stold(args[1]), std::stold(args[2]));
    }
    else
    {
        throw std::invalid_argument("spacetime with name " + spacetime_name + " does not exist");
    }
    return spacetime;
}

std::vector<std::string> Interface::find_function_arguments(std::string text)
{
    int i, counter = 0;
    bool start = false;
    int start_index, end_index;
    int word_start, word_end;
    std::vector<std::string> result;
    
    for (i = 0; i < text.length(); i++)
    {
        // counting parenthesis
        if (text[i]=='(')
            counter++;
        else if (text[i]==')')
            counter--;

        // controling parenthesis
        if (counter == 1 && !start)
        {
            start = true;
            start_index=i;
            word_start = i+1;
        }
        else if (start && counter == 0)
        {
            end_index = i;
            break;
        }
        else if (counter < 0)
        {
            throw("invalid usage of parenthesis");
        }

        // getting arguments
        if (start && (counter==1 && text[i] == ','))
        {
            result.push_back(strip(text.substr(word_start, i-word_start)));
            word_start=i+1;
        }
    }
    if (start)
        result.push_back(strip(text.substr(word_start, i-word_start)));
    if (!start)
        throw std::invalid_argument("no argument to function were given");
    return result;
}

bool Interface::try_apply_function(std::string text)
{
    std::string name, rest;
    find_command_name(text, name, rest);
    if (name == "split_args")
    {
        this->split_args(rest);
        return true;
    }
    else if (name == "draw_potential_1D")
    {
        this->draw_potential_1D(rest);
        return true;
    }
    else if (name == "draw_lambda_1D")
    {
        this->draw_lambda_1D(rest);
        return true;
    }
    else if (name == "solve_ode_system")
    {
        // this->solve_ode_system(rest);
        return true;
    }
    else if (name == "local_expansions_weyl")
    {
        this->local_expansions_weyl(rest);
        return true;
    }
    else if (name == "norm_growth_weyl")
    {
        this->norm_growth_weyl(rest);
        return true;
    }
    else if (name == "rest_norm2_weyl")
    {
        this->rest_norm2_weyl(rest);
        return true;
    }
    else if (name == "poincare_border_weyl")
    {
        this->poincare_border_weyl(rest);
        return true;
    }
    else if (name == "poincare_section_weyl")
    {
        this->poincare_section_weyl(rest);
        return true;
    }
    else if (name == "numerical_expansions_weyl")
    {
        this->numerical_expansions_weyl(rest);
        return true;
    }
    else if (name == "trajectory_weyl")
    {
        this->trajectory_weyl(rest);
        return true;
    }
    return false;
}

void Interface::split_args(std::string text)
{
    auto args = find_function_arguments(text);
    for (int i = 0; i<args.size(); i++)
        std::cout << "arg" << i << ": " << args[i] << std::endl;
}

void Interface::draw_potential_1D(std::string text)
{
    //weylspacetime, coordiante, initial_values, min, max, num, file

    // get arguments
    auto args = find_function_arguments(text);
    int number_of_arguments = 7;
    if (args.size() > number_of_arguments)
        throw std::invalid_argument("too much arguments for draw_potential_1D");
    else if (args.size() < number_of_arguments)
        throw std::invalid_argument("too little arguments for draw_potential_1D");

    std::shared_ptr<gr2::Weyl> spacetime = nullptr;
    std::ofstream file;

    try
    {
        // save arguments
        spacetime = this->create_weyl_spacetime(args[0]);
        int coordinate = std::stoi(args[1]);
        auto coordinates = find_function_arguments(args[2]);
        gr2::real min_val = std::stold(args[3]);
        if (coordinates.size() != 4)
            throw std::invalid_argument("invalid number of coordinates");
        gr2::real max_val = std::stold(args[4]);
        int num = std::stoi(args[5]);
        std::string file_name = args[6];

        // prepare calculation
        gr2::real y[4];
        for (int i = 0; i <4; i++)
        {
            y[i] = std::stold(coordinates[i]);
        }

        // calculation
        file.open(file_name);
        for (int i = 0; i< num; i++)
        {
            y[coordinate] = (max_val-min_val)/(num-1)*i + min_val;
            spacetime->calculate_nu(y);
            file << y[coordinate] << ";" << spacetime->get_nu() << std::endl;
        }
    }
    catch(const std::exception& e)
    {
        file.close();
        throw e;
    }
    //TODO: try catch, spacetime is dynamicaly allocated, file has to be closed
    // procede calculation
}

void Interface::draw_lambda_1D(std::string text)
{
    //weylspacetime, coordiante, min, max, num, file

    // get arguments
    auto args = find_function_arguments(text);
    int number_of_arguments = 7;
    if (args.size() > number_of_arguments)
        throw std::invalid_argument("too much arguments for draw_potential_1D");
    else if (args.size() < number_of_arguments)
        throw std::invalid_argument("too little arguments for draw_potential_1D");

    std::shared_ptr<gr2::Weyl> spacetime = nullptr;
    std::ofstream file;

    try
    {
        // save arguments
        spacetime = this->create_weyl_spacetime(args[0]);
        int coordinate = std::stoi(args[1]);
        auto coordinates = find_function_arguments(args[2]);
        gr2::real min_val = std::stold(args[3]);
        if (coordinates.size() != 4)
            throw std::invalid_argument("invalid number of coordinates");
        gr2::real max_val = std::stold(args[4]);
        int num = std::stoi(args[5]);
        std::string file_name = args[6];

        // prepare calculation
        gr2::real y[4];
        for (int i = 0; i <4; i++)
        {
            y[i] = std::stold(coordinates[i]);
        }

        // calculation
        file.open(file_name);
        for (int i = 0; i< num; i++)
        {
            y[coordinate] = (max_val-min_val)/(num-1)*i + min_val;
            spacetime->calculate_lambda_init(y);
            file << y[coordinate] << ";" << spacetime->get_lambda() << std::endl;
        }
    }
    catch(const std::exception& e)
    {
        file.close();
        throw e;
    }
}

std::shared_ptr<gr2::OdeSystem> Interface::create_ode_system(std::string text)
{
    text = strip(text);
    std::string ode_name, args_text;
    this->find_command_name(text, ode_name, args_text);
    auto args = this->find_function_arguments(args_text);
    std::shared_ptr<gr2::OdeSystem> ode;

    if (ode_name == "DampedHarmonicOscillator")
    {
        if (args.size() != 2)
            throw std::invalid_argument("invalid number of arguments for DampedHarmonicOscillator");
        ode = std::make_shared<gr2::DampedHarmonicOscillator>(std::stold(args[0]), std::stold(args[1]));
    }
    else
    {
        throw std::invalid_argument("spacetime with name " + ode_name + "does not exists");
    }
    return ode;
}

// void Interface::solve_ode_system(std::string text)
// {
//     auto args = find_function_arguments(text);
//     int number_of_arguments = 7;
//     if (args.size() < number_of_arguments)
//         throw std::invalid_argument("too little arguments for solve_ode_system");
//     else if (args.size() > number_of_arguments)
//         throw std::invalid_argument("too much arguments for solve_ode_system");
    
//     std::shared_ptr<gr2::OdeSystem> ode = nullptr;
//     std::ofstream file;

//     try
//     {
//         /* code */
//         ode = this->create_ode_system(args[0]);
//         auto initial_conditions = find_function_arguments(args[1]);
//         if (initial_conditions.size() != ode->get_n())
//             throw std::invalid_argument("invalid number of initial_value_conditions");
//         gr2::real t_start = std::stold(args[2]);
//         gr2::real t_end = std::stold(args[3]);
//         gr2::real delta_t = std::stold(args[4]);
//         std::string method = args[5];
//         std::string file_name = args[6];

//         // prepare initial conditions
//         gr2::real *y_initial = new gr2::real [ode->get_n()];
//         for (int i = 0; i < ode->get_n(); i++)
//             y_initial[i] = std::stold(initial_conditions[i]);

//         // calculation
//         gr2::Integrator integrator(*ode, method);
//         DataRecord recorder(ode->get_n());
//         integrator.add_event(&recorder);
//         integrator.integrate(y_initial, t_start, t_end, delta_t);

//         // saving data
//         file.open(file_name);
//         for (auto record:recorder.data)
//         {
//             for (auto d : record)
//                 file << d << ";";
//             file << std::endl;
//         }
//     }
//     catch(const std::exception& e)
//     {
//         file.close();
//         throw e;
//     }
// }

void Interface::local_expansions_weyl(std::string text)
{
    // Initialize calculation
    auto args = find_function_arguments(text);
    int number_of_arguments = 7;
    if (args.size() < number_of_arguments)
        throw std::invalid_argument("too little arguments for local_expansions_Weyl");
    else if (args.size() > number_of_arguments)
        throw std::invalid_argument("too much arguments for local_expansions_Weyl");

    std::shared_ptr<gr2::Weyl> spt = this->create_weyl_spacetime(args[0]);

    gr2::real E = std::stold(args[1]);
    gr2::real L = std::stold(args[2]);
    auto range_rho = find_function_arguments(args[3]);
    if (range_rho.size() != 3)
        throw std::invalid_argument("incorent number of arguments for range in rho");
    gr2::real rho_min = std::stold(range_rho[0]);
    gr2::real rho_max = std::stold(range_rho[1]);
    int n_rho = std::stoi(range_rho[2]);
    gr2::real delta_rho = (rho_max-rho_min)/(n_rho-1);

    auto range_z = find_function_arguments(args[4]);
    if (range_z.size() != 3)
        throw std::invalid_argument("incorent number of arguments for range in z");
    gr2::real z_min = std::stold(range_z[0]);
    gr2::real z_max = std::stold(range_z[1]);
    int n_z = std::stoi(range_z[2]);
    gr2::real delta_z = (z_max-z_min)/(n_z-1);

    int n_angles = std::stoi(args[5]);
    gr2::real delta_angle = 2*gr2::pi/n_angles;
    std::string file_name = args[6];

    std::ofstream file;
    gr2::real y[9]={};
    gsl_matrix* H;
    gsl_vector_complex *eig_vals;
    gsl_eigen_nonsymm_workspace *work_space;

    // Procede in calculation
    try
    {
        // open file
        file.open(file_name);

        if (!file.is_open())
            throw std::runtime_error("file " + file_name + "could not be opened");

        // prepare calculation of eigenvalues
        eig_vals = gsl_vector_complex_alloc(8);
        work_space = gsl_eigen_nonsymm_alloc(8);

        // calculate eigenvalues
        for (int i = 0; i < n_rho; i++)
            for (int j = 0; j < n_z; j++)
            {
                gr2::real rho = rho_min + i*delta_rho;
                gr2::real z = z_min + j*delta_z;
                y[gr2::Weyl::RHO] = rho;
                y[gr2::Weyl::Z] = z;

                // calculate lambda 
                spt->calculate_lambda_init(y);
                y[gr2::Weyl::LAMBDA] = spt->get_lambda();

                // calculate ut (from E)
                spt->calculate_metric(y);
                y[gr2::Weyl::UT] = -E/spt->get_metric()[gr2::Weyl::T][gr2::Weyl::T];

                // calculate uphi (from L)
                y[gr2::Weyl::UPHI] = L/spt->get_metric()[gr2::Weyl::PHI][gr2::Weyl::PHI];

                // calculate size of rest velocity
                gr2::real norm2 = (-1 + y[gr2::Weyl::UT]*E - y[gr2::Weyl::UPHI]*L);
                if (norm2 < 0)
                    continue;
                gr2::real norm2_c = norm2/spt->get_metric()[gr2::Weyl::RHO][gr2::Weyl::RHO];
                gr2::real norm_c = sqrtl(norm2_c);
            
                gr2::real method_value = 0;
                for (int k = 0; k < n_angles; k++)
                {
                    gr2::real value = 0;
                    gr2::real angle = k*delta_angle;

                    // calculate rest of velocity
                    y[gr2::Weyl::URHO] = norm_c*cosl(angle);
                    y[gr2::Weyl::UZ] = norm_c*sinl(angle);

                    // calculate matrix H
                    H = gr2::time_corrected_matrix_H(spt.get(), y);

                    // calculate eigen values
                    gsl_eigen_nonsymm(H, eig_vals, work_space);
                    
                    // find the greatest real part
                    for (int i = 0; i < 8; i++)
                        value = std::max((double)value, GSL_REAL(gsl_vector_complex_get(eig_vals, i)));
                    method_value = std::max(value, method_value);

                    // free the matrix
                    gsl_matrix_free(H);
                }
                // save values to the file
                file << i << ";" << j << ";" << rho << ";" << z << ";" << method_value << "\n";
            }

            // free the work space
            gsl_vector_complex_free(eig_vals);
            gsl_eigen_nonsymm_free(work_space);

            // close file
            file.close();
    }
    catch(const std::exception& e)
    {
        gsl_matrix_free(H);
        gsl_vector_complex_free(eig_vals);
        gsl_eigen_nonsymm_free(work_space);
        file.close();
        throw e;
    }
};

void Interface::norm_growth_weyl(std::string text)
{
    // Initialize calculation
    auto args = find_function_arguments(text);
    int number_of_arguments = 7;
    if (args.size() < number_of_arguments)
        throw std::invalid_argument("too little arguments for local_expansions_Weyl");
    else if (args.size() > number_of_arguments)
        throw std::invalid_argument("too much arguments for local_expansions_Weyl");

    std::shared_ptr<gr2::Weyl> spt = this->create_weyl_spacetime(args[0]);

    gr2::real E = std::stold(args[1]);
    gr2::real L = std::stold(args[2]);
    auto range_rho = find_function_arguments(args[3]);
    if (range_rho.size() != 3)
        throw std::invalid_argument("incorent number of arguments for range in rho");
    gr2::real rho_min = std::stold(range_rho[0]);
    gr2::real rho_max = std::stold(range_rho[1]);
    int n_rho = std::stoi(range_rho[2]);
    gr2::real delta_rho = (rho_max-rho_min)/(n_rho-1);

    auto range_z = find_function_arguments(args[4]);
    if (range_z.size() != 3)
        throw std::invalid_argument("incorent number of arguments for range in z");
    gr2::real z_min = std::stold(range_z[0]);
    gr2::real z_max = std::stold(range_z[1]);
    int n_z = std::stoi(range_z[2]);
    gr2::real delta_z = (z_max-z_min)/(n_z-1);

    int n_angles = std::stoi(args[5]);
    gr2::real delta_angle = 2*gr2::pi/n_angles;
    std::string file_name = args[6];

    std::ofstream file;
    gr2::real y[9]={};
    gsl_matrix *H;
    gsl_vector *eig_vals;
    gsl_eigen_symm_workspace *work_space;

    // Procede in calculation
    try
    {
        // open file
        file.open(file_name);

        if (!file.is_open())
            throw std::runtime_error("file " + file_name + "could not be opened");

        // prepare calculation of eigenvalues
        eig_vals = gsl_vector_alloc(8);
        work_space = gsl_eigen_symm_alloc(8);

        // calculate eigenvalues
        for (int i = 0; i < n_rho; i++)
            for (int j = 0; j < n_z; j++)
            {
                gr2::real rho = rho_min + i*delta_rho;
                gr2::real z = z_min + j*delta_z;
                y[gr2::Weyl::RHO] = rho;
                y[gr2::Weyl::Z] = z;

                // calculate lambda 
                spt->calculate_lambda_init(y);
                y[gr2::Weyl::LAMBDA] = spt->get_lambda();

                // calculate ut (from E)
                spt->calculate_metric(y);
                y[gr2::Weyl::UT] = -E/spt->get_metric()[gr2::Weyl::T][gr2::Weyl::T];

                // calculate uphi (from L)
                y[gr2::Weyl::UPHI] = L/spt->get_metric()[gr2::Weyl::PHI][gr2::Weyl::PHI];

                // calculate size of rest velocity
                gr2::real norm2 = (-1 + y[gr2::Weyl::UT]*E - y[gr2::Weyl::UPHI]*L);
                if (norm2 < 0)
                    continue;
                gr2::real norm2_c = norm2/spt->get_metric()[gr2::Weyl::RHO][gr2::Weyl::RHO];
                gr2::real norm_c = sqrtl(norm2_c);
            
                gr2::real method_value = 0;
                for (int k = 0; k < n_angles; k++)
                {
                    gr2::real value = 0;
                    gr2::real angle = k*delta_angle;

                    // calculate rest of velocity
                    y[gr2::Weyl::URHO] = norm_c*cosl(angle);
                    y[gr2::Weyl::UZ] = norm_c*sinl(angle);

                    // calculate matrix H
                    H = gr2::time_corrected_matrix_H(spt.get(), y);

                    // calculate eigen values
                    gsl_eigen_symm(H, eig_vals, work_space);
                    
                    // find the greatest real part
                    for (int i = 0; i < 8; i++)
                        value = std::max((double)value, gsl_vector_get(eig_vals, i));
                    method_value = std::max(value, method_value);

                    // free the matrix
                    gsl_matrix_free(H);
                }
                // save values to the file
                file << i << ";" << j << ";" << rho << ";" << z << ";" << method_value << "\n";
            }

            // free the work space
            gsl_vector_free(eig_vals);
            gsl_eigen_symm_free(work_space);

            // close file
            file.close();
    }
    catch(const std::exception& e)
    {
        gsl_matrix_free(H);
        gsl_vector_free(eig_vals);
        gsl_eigen_symm_free(work_space);
        file.close();
        throw e;
    }
};

void Interface::rest_norm2_weyl(std::string text)
{
    // Initialize calculation
    auto args = find_function_arguments(text);
    int number_of_arguments = 6;
    if (args.size() < number_of_arguments)
        throw std::invalid_argument("too little arguments for local_expansions_Weyl");
    else if (args.size() > number_of_arguments)
        throw std::invalid_argument("too much arguments for local_expansions_Weyl");

    std::shared_ptr<gr2::Weyl> spt = this->create_weyl_spacetime(args[0]);

    gr2::real E = std::stold(args[1]);
    gr2::real L = std::stold(args[2]);
    auto range_rho = find_function_arguments(args[3]);
    if (range_rho.size() != 3)
        throw std::invalid_argument("incorent number of arguments for range in rho");
    gr2::real rho_min = std::stold(range_rho[0]);
    gr2::real rho_max = std::stold(range_rho[1]);
    int n_rho = std::stoi(range_rho[2]);
    gr2::real delta_rho = (rho_max-rho_min)/(n_rho-1);

    auto range_z = find_function_arguments(args[4]);
    if (range_z.size() != 3)
        throw std::invalid_argument("incorent number of arguments for range in z");
    gr2::real z_min = std::stold(range_z[0]);
    gr2::real z_max = std::stold(range_z[1]);
    int n_z = std::stoi(range_z[2]);
    gr2::real delta_z = (z_max-z_min)/(n_z-1);
    
    std::string file_name = args[5];

    std::ofstream file;
    gr2::real y[9]={};

    // Procede in calculation
    try
    {
        // open file
        file.open(file_name);

        if (!file.is_open())
            throw std::runtime_error("file " + file_name + "could not be opened");

        // calculate norms
        for (int i = 0; i < n_rho; i++)
            for (int j = 0; j < n_z; j++)
            {
                gr2::real rho = rho_min + i*delta_rho;
                gr2::real z = z_min + j*delta_z;
                y[gr2::Weyl::RHO] = rho;
                y[gr2::Weyl::Z] = z;

                // calculate lambda 
                spt->calculate_lambda_init(y);
                y[gr2::Weyl::LAMBDA] = spt->get_lambda();

                // calculate ut (from E)
                spt->calculate_metric(y);
                y[gr2::Weyl::UT] = -E/spt->get_metric()[gr2::Weyl::T][gr2::Weyl::T];

                // calculate uphi (from L)
                y[gr2::Weyl::UPHI] = L/spt->get_metric()[gr2::Weyl::PHI][gr2::Weyl::PHI];

                // calculate size of rest velocity
                gr2::real norm2 = (-1 + y[gr2::Weyl::UT]*E - y[gr2::Weyl::UPHI]*L);
            
                // save values to the file
                file << i << ";" << j << ";" << rho << ";" << z << ";" << norm2 << "\n";
            }

            // close file
            file.close();
    }
    catch(const std::exception& e)
    {
        file.close();
        throw e;
    }
}

void Interface::poincare_border_weyl(std::string text)
{
    // Initialize calculation
    auto args = find_function_arguments(text);
    int number_of_arguments = 5;
    if (args.size() < number_of_arguments)
        throw std::invalid_argument("too little arguments for local_expansions_Weyl");
    else if (args.size() > number_of_arguments)
        throw std::invalid_argument("too much arguments for local_expansions_Weyl");

    std::shared_ptr<gr2::Weyl> spt = this->create_weyl_spacetime(args[0]);

    gr2::real E = std::stold(args[1]);
    gr2::real L = std::stold(args[2]);
    auto range_rho = find_function_arguments(args[3]);
    if (range_rho.size() != 3)
        throw std::invalid_argument("incorent number of arguments for range in rho");
    gr2::real rho_min = std::stold(range_rho[0]);
    gr2::real rho_max = std::stold(range_rho[1]);
    int n_rho = std::stoi(range_rho[2]);
    gr2::real delta_rho = (rho_max-rho_min)/(n_rho-1);

    std::string file_name = args[4];

    std::ofstream file;
    gr2::real y[9]={};

    // Procede in calculation
    try
    {
        // open file
        file.open(file_name);

        if (!file.is_open())
            throw std::runtime_error("file " + file_name + "could not be opened");

        // calculate norms
        for (int i = 0; i < n_rho; i++)
        {
            gr2::real rho = rho_min + i*delta_rho;
            gr2::real z = 1e-5;
            y[gr2::Weyl::RHO] = rho;
            y[gr2::Weyl::Z] = z;

            // calculate lambda 
            spt->calculate_lambda_init(y);
            y[gr2::Weyl::LAMBDA] = spt->get_lambda();

            // calculate ut (from E)
            spt->calculate_metric(y);
            y[gr2::Weyl::UT] = -E/spt->get_metric()[gr2::Weyl::T][gr2::Weyl::T];

            // calculate uphi (from L)
            y[gr2::Weyl::UPHI] = L/spt->get_metric()[gr2::Weyl::PHI][gr2::Weyl::PHI];

            // calculate size of rest velocity
            gr2::real norm2 = (-1 + y[gr2::Weyl::UT]*E - y[gr2::Weyl::UPHI]*L);
            gr2::real urho;
            if (norm2<0)
                urho = 0;
            else
                urho = sqrtl(norm2/spt->get_metric()[gr2::Weyl::RHO][gr2::Weyl::RHO]);
        
            // save values to the file
            file << i << ";" << rho << ";" << urho << "\n";
        }
        // close file
        file.close();
    }
    catch(const std::exception& e)
    {
        file.close();
        throw e;
    }
}

void Interface::poincare_section_weyl(std::string text)
{
    // Initialize calculation
    auto args = find_function_arguments(text);
    int number_of_arguments = 7;
    if (args.size() < number_of_arguments)
        throw std::invalid_argument("too little arguments for local_expansions_Weyl");
    else if (args.size() > number_of_arguments)
        throw std::invalid_argument("too much arguments for local_expansions_Weyl");

    std::shared_ptr<gr2::Weyl> spt = this->create_weyl_spacetime(args[0]);

    gr2::real E = std::stold(args[1]);
    gr2::real L = std::stold(args[2]);
    auto range_rho = find_function_arguments(args[3]);
    if (range_rho.size() != 3)
        throw std::invalid_argument("incorent number of arguments for range in rho");
    gr2::real rho_min = std::stold(range_rho[0]);
    gr2::real rho_max = std::stold(range_rho[1]);
    int n_rho = std::stoi(range_rho[2]);
    gr2::real delta_rho = (rho_max-rho_min)/(n_rho-1);

    int angles = std::stoi(args[4]);
    gr2::real delta_angle = gr2::pi_4/angles;
    gr2::real t_max = std::stoll(args[5]);
    std::string file_name = args[6];

    std::ofstream file;
    gr2::real y[9]={};

    // Procede in calculation
    try
    {
        // open file
        file.open(file_name);

        if (!file.is_open())
            throw std::runtime_error("file " + file_name + "could not be opened");

        gr2::Integrator integrator(spt, "DoPr853", 1e-16, 1e-16, false);
        auto too_close = std::make_shared<StopBeforeBlackHole>(0.4);
        integrator.add_event(too_close);
        auto errorE_too_high = std::make_shared<StopTooHighErrorE>(spt,E,1e-10);
        integrator.add_event(errorE_too_high);
        auto errorL_too_high = std::make_shared<StopTooHighErrorL>(spt,L,1e-10);
        integrator.add_event(errorL_too_high);
        auto stop_on_disk = std::make_shared<StopOnDisk>(spt, 1e-4, true);
        integrator.add_event(stop_on_disk);

        // calculate norms
        for (int i = 0; i < n_rho; i++)
        {
            gr2::real rho = rho_min + i*delta_rho;
            gr2::real z = 1e-3;
            y[gr2::Weyl::RHO] = rho;
            y[gr2::Weyl::Z] = z;

            // calculate lambda 
            spt->calculate_lambda_init(y);
            y[gr2::Weyl::LAMBDA] = spt->get_lambda();

            // calculate ut (from E)
            spt->calculate_metric(y);
            y[gr2::Weyl::UT] = -E/spt->get_metric()[gr2::Weyl::T][gr2::Weyl::T];

            // calculate uphi (from L)
            y[gr2::Weyl::UPHI] = L/spt->get_metric()[gr2::Weyl::PHI][gr2::Weyl::PHI];

            // calculate size of rest velocity
            gr2::real norm2 = (-1 + y[gr2::Weyl::UT]*E - y[gr2::Weyl::UPHI]*L);
            if (norm2 < 0)
                continue;
            gr2::real norm = sqrtl(norm2/spt->get_metric()[gr2::Weyl::RHO][gr2::Weyl::RHO]);

            std::cout << "rho" << y[gr2::Weyl::RHO] << std::endl;
            for (int j = 0; j < angles; j++)
            {
                // calculate initial conditions
                gr2::real angle = j*delta_angle;
                y[gr2::Weyl::URHO] = norm*sinl(angle);
                y[gr2::Weyl::UZ] = norm*cosl(angle);

                // calculate poincare section
                try
                {
                    /* code */
                    integrator.integrate(y, 0, t_max, 0.2);
                }
                catch(const std::exception& e)
                {
                    std::cerr << e.what() << '\n';
                }
                
                // save data
                for (auto& d : stop_on_disk->data)
                    file << d[0] << ";" << d[1] << "\n";

                // delete data
                stop_on_disk->data.clear();
                file.flush();
            }
        }
        // close file
        file.close();
    }
    catch(const std::exception& e)
    {
        file.close();
        throw e;
    }
}

void Interface::numerical_expansions_weyl(std::string text)
{
    // Initialize calculation
    auto args = find_function_arguments(text);
    int number_of_arguments = 10;
    if (args.size() < number_of_arguments)
        throw std::invalid_argument("too little arguments for local_expansions_Weyl");
    else if (args.size() > number_of_arguments)
        throw std::invalid_argument("too much arguments for local_expansions_Weyl");

    std::shared_ptr<gr2::Weyl> spt = this->create_weyl_spacetime(args[0]);
    std::shared_ptr<gr2::OdeSystem> ode = std::make_shared<gr2::CombinedOdeSystem>(std::vector<std::shared_ptr<gr2::OdeSystem>>{spt, spt});

    gr2::real E = std::stold(args[1]);
    gr2::real L = std::stold(args[2]);
    auto range_rho = find_function_arguments(args[3]);
    if (range_rho.size() != 3)
        throw std::invalid_argument("incorent number of arguments for range in rho");
    gr2::real rho_min = std::stold(range_rho[0]);
    gr2::real rho_max = std::stold(range_rho[1]);
    int n_rho = std::stoi(range_rho[2]);
    gr2::real delta_rho = (rho_max-rho_min)/(n_rho-1);

    auto range_z = find_function_arguments(args[4]);
    if (range_z.size() != 3)
        throw std::invalid_argument("incorent number of arguments for range in z");
    gr2::real z_min = std::stold(range_z[0]);
    gr2::real z_max = std::stold(range_z[1]);
    int n_z = std::stoi(range_z[2]);
    gr2::real delta_z = (z_max-z_min)/(n_z-1);
    
    gr2::real rho_start = std::stold(args[5]);
    gr2::real u_rho_frac = std::stold(args[6]);

    gr2::real t_max = std::stold(args[7]);
    std::string file_name = args[8];

    gr2::real eps_pos = 1e-8;

    std::ofstream file;
    gr2::real y[19]={};
    
    // Procede in calculation
    try
    {
        // open file
        file.open(file_name);

        if (!file.is_open())
            throw std::runtime_error("file " + file_name + "could not be opened");

        gr2::Integrator integrator(ode, "DoPr853", 1e-16, 1e-16, true);
        auto too_close = std::make_shared<StopBeforeBlackHole>(0.4);
        integrator.add_event(too_close);
        auto errorE_too_high = std::make_shared<StopTooHighErrorE>(spt,E,1e-10);
        integrator.add_event(errorE_too_high);
        auto errorL_too_high = std::make_shared<StopTooHighErrorL>(spt,L,1e-10);
        integrator.add_event(errorL_too_high);
        auto stop_on_disk = std::make_shared<StopOnDiskTwoParticles>(spt, 1e-4);
        integrator.add_event(stop_on_disk);
        auto renormalization = std::make_shared<RenormalizationOfSecondParticleWeyl>(spt, 1e-5);
        auto num_expansions = std::make_shared<NumericalExpansions>(spt, rho_min, rho_max, n_rho, z_min, z_max, n_z, &(renormalization->log_norm));
        integrator.add_event(num_expansions);
        integrator.add_event(renormalization);

        // ========== initial conditions for the first particle ==========
        gr2::real rho = rho_start;
        gr2::real z = 1e-3;
        y[gr2::Weyl::RHO] = rho;
        y[gr2::Weyl::Z] = z;

        // calculate lambda 
        spt->calculate_lambda_init(y);
        y[gr2::Weyl::LAMBDA] = spt->get_lambda();

        // calculate ut (from E)
        spt->calculate_metric(y);
        y[gr2::Weyl::UT] = -E/spt->get_metric()[gr2::Weyl::T][gr2::Weyl::T];

        // calculate uphi (from L)
        y[gr2::Weyl::UPHI] = L/spt->get_metric()[gr2::Weyl::PHI][gr2::Weyl::PHI];

        // calculate size of rest velocity
        gr2::real norm2 = (-1 + y[gr2::Weyl::UT]*E - y[gr2::Weyl::UPHI]*L);
        if (norm2 < 0)
            return;
        gr2::real norm = sqrtl(norm2/spt->get_metric()[gr2::Weyl::RHO][gr2::Weyl::RHO]);

        // calculate velocity
        y[gr2::Weyl::URHO] = norm*u_rho_frac;
        y[gr2::Weyl::UZ] = norm*sqrtl(1-u_rho_frac*u_rho_frac);

        // ========== initial conditions for second particle ==========
        gr2::real* y_ = y + 9;

        // positions
        for (int j = 0; j < 4; j++)
            y_[j] = y[j] + eps_pos;

        // calculate lambda 
        spt->calculate_lambda_init(y_);
        y_[gr2::Weyl::LAMBDA] = spt->get_lambda();

        // calculate ut (from E)
        spt->calculate_metric(y_);
        y_[gr2::Weyl::UT] = -E/spt->get_metric()[gr2::Weyl::T][gr2::Weyl::T];

        // calculate uphi (from L)
        y_[gr2::Weyl::UPHI] = L/spt->get_metric()[gr2::Weyl::PHI][gr2::Weyl::PHI];

        // calculate size of rest velocity
        norm2 = (-1 + y_[gr2::Weyl::UT]*E - y_[gr2::Weyl::UPHI]*L);
        if (norm2 < 0)
            return;
        norm = sqrtl(norm2/spt->get_metric()[gr2::Weyl::RHO][gr2::Weyl::RHO]);

        // calculate velocity
        y_[gr2::Weyl::URHO] = norm*u_rho_frac;
        y_[gr2::Weyl::UZ] = norm*sqrtl(1-u_rho_frac*u_rho_frac);

        // init norm
        y[18] = 0;

        // calculate numerical expansions
        try
        {
            /* code */
            std::cout << "Integrujeme" << std::endl;
            integrator.integrate(y, 0, t_max, 0.2);
            std::cout << "Dointegrovano" << std::endl;
        }
        catch(const std::exception& e)
        {
            std::cerr << e.what() << '\n';
        }
        
        // TODO: save data
        std::cout << "Jdeme ukladat" << std::endl;
        for (int i = 0; i<n_rho; i++)
            for (int j = 0; j<n_z; j++)
            {
                file << i << ";" << j << ";" << rho_min + i*delta_rho << ";" << rho_max + j*delta_z << ";" << num_expansions->data[i][j] << std::endl;
            }

        file.flush();
        // close file
        file.close();
    }
    catch(const std::exception& e)
    {
        file.close();
        throw e;
    }
}

void Interface::trajectory_weyl(std::string text)
{
    // Initialize calculation
    auto args = find_function_arguments(text);
    int number_of_arguments = 8;
    if (args.size() < number_of_arguments)
        throw std::invalid_argument("too little arguments for local_expansions_Weyl");
    else if (args.size() > number_of_arguments)
        throw std::invalid_argument("too much arguments for local_expansions_Weyl");

    std::shared_ptr<gr2::Weyl> spt = this->create_weyl_spacetime(args[0]);

    gr2::real E = std::stold(args[1]);
    gr2::real L = std::stold(args[2]);
    
    gr2::real rho_start = std::stold(args[3]);
    gr2::real u_rho_frac = std::stold(args[4]);

    gr2::real t_max = std::stold(args[5]);
    gr2::real dt = std::stold(args[6]);
    std::string file_name = args[7];

    std::ofstream file;
    gr2::real y[9]={};
    
    // Procede in calculation
    try
    {
        // open file
        file.open(file_name);

        if (!file.is_open())
            throw std::runtime_error("file " + file_name + "could not be opened");

        gr2::Integrator integrator(spt, "DoPr853", 1e-16, 1e-16, true);
        auto data_monitor = std::make_shared<ConstantStepDataMonitoring>(0, dt);
        integrator.add_event(data_monitor);
        auto too_close = std::make_shared<StopBeforeBlackHole>(0.4);
        integrator.add_event(too_close);
        auto errorE_too_high = std::make_shared<StopTooHighErrorE>(spt,E,1e-10);
        integrator.add_event(errorE_too_high);
        auto errorL_too_high = std::make_shared<StopTooHighErrorL>(spt,L,1e-10);
        integrator.add_event(errorL_too_high);
        auto stop_on_disk = std::make_shared<StopOnDisk>(spt, 1e-4, false);
        integrator.add_event(stop_on_disk);

        // ========== initial conditions for the first particle ==========
        gr2::real rho = rho_start;
        gr2::real z = 1e-3;
        y[gr2::Weyl::RHO] = rho;
        y[gr2::Weyl::Z] = z;

        // calculate lambda 
        spt->calculate_lambda_init(y);
        y[gr2::Weyl::LAMBDA] = spt->get_lambda();

        // calculate ut (from E)
        spt->calculate_metric(y);
        y[gr2::Weyl::UT] = -E/spt->get_metric()[gr2::Weyl::T][gr2::Weyl::T];

        // calculate uphi (from L)
        y[gr2::Weyl::UPHI] = L/spt->get_metric()[gr2::Weyl::PHI][gr2::Weyl::PHI];

        // calculate size of rest velocity
        gr2::real norm2 = (-1 + y[gr2::Weyl::UT]*E - y[gr2::Weyl::UPHI]*L);
        if (norm2 < 0)
            return;
        gr2::real norm = sqrtl(norm2/spt->get_metric()[gr2::Weyl::RHO][gr2::Weyl::RHO]);

        // calculate velocity
        y[gr2::Weyl::URHO] = norm*u_rho_frac;
        y[gr2::Weyl::UZ] = norm*sqrtl(1-u_rho_frac*u_rho_frac);

        // calculate trajectory
        try
        {
            /* code */
            std::cout << "Integrujeme" << std::endl;
            integrator.integrate(y, 0, t_max, 0.2);
            std::cout << "Konec integrace" << std::endl;
        }
        catch(const std::exception& e)
        {
            std::cerr << e.what() << '\n';
        }
        
        // TODO: save data
        std::cout << "Jdeme ukladat" << std::endl;
        for (auto &d:data_monitor->data)
        {
            file << d[0];
            for (int i = 1; i < 9; i++)
                file << ";" << d[i];
            file << std::endl;
        }

        file.flush();
        // close file
        file.close();
    }
    catch(const std::exception& e)
    {
        file.close();
        throw e;
    }
}

Interface::Interface():macros(), values(), help_name(), help_text()
{
    // load help
    std::ifstream file;
    std::string filename = "./data/help.txt";

    file.open(filename);

    if (!file.is_open())
        throw std::runtime_error("Program was not able to load help");

    std::string line, help_name_string, help_text_string;

    while (getline(file, help_name_string))
    {
        while (getline(file, line))
        {
            if (line[0] == '=')
            {
                this->help_name.push_back(this->strip(help_name_string));
                this->help_text.push_back(help_text_string);
                help_text_string = "";
                break;
            }
            help_text_string += help_text_string==""?line:"\n" + line;
        }
    }
}

bool Interface::command(std::string text)
{
    // skip line
    if (text == "")
        return true;

    // end
    text = strip(text);
    if (text == "end" || text == "END" || text == "exit")
        return false;

    // try to apply def, delete
    if(try_apply_operators(text))
        return true;

    std::string new_text = substitute(text);

    // try to apply functions
    if(try_apply_function(new_text))
        return true;

    throw std::invalid_argument("command could not be recognized");
}