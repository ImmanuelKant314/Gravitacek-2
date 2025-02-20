#include "interface/interface.hpp"
#include "gravitacek2/geomotion/spacetimes.hpp"
#include "gravitacek2/integrator/integrator.hpp"
#include "gravitacek2/integrator/odesystems.hpp"
#include <stdexcept>
#include <iostream>
#include <algorithm>
#include <fstream>
#include <array>

class DataRecord : public gr2::Event
{
protected:
    int n;

public:
    std::vector<std::vector<gr2::real>> data;

    DataRecord(int n) : gr2::Event(gr2::EventType::data)
    {
        this->n = n;
    }

    virtual gr2::real value(const gr2::real &t, const gr2::real y[], const gr2::real dydt[]) override
    {
        return 0;
    }

    virtual void apply(gr2::real &t, gr2::real y[], gr2::real dydt[])
    {
        std::vector<gr2::real> record;
        record.push_back(t);
        for (int i=0; i < n; i++)
            record.push_back(y[i]);
        this->data.push_back(record);
    }
};

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

gr2::Weyl* Interface::create_weyl_spacetime(std::string text)
{
    text = strip(text);
    std::string spacetime_name, args_text;
    this->find_command_name(text, spacetime_name, args_text);
    auto args = this->find_function_arguments(args_text);
    gr2::Weyl* spacetime;

    if (spacetime_name == "WeylSchwarzschild")
    {
        if (args.size() != 1)
            throw std::invalid_argument("invalid number of arguments for WeylSchwarzschild");
        spacetime = new gr2::WeylSchwarzschild(std::stold(args[0]));
    }
    else if (spacetime_name == "BachWeylRing")
    {   
        if (args.size() != 2)
            throw std::invalid_argument("invalid number of arguments for BachWeylRing");
        spacetime = new gr2::BachWeylRing(std::stold(args[0]), std::stold(args[1]));
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
        this->solve_ode_system(rest);
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

    gr2::Weyl* spacetime = nullptr;
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
        delete[] spacetime;
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

    gr2::Weyl* spacetime = nullptr;
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
        delete[] spacetime;
        file.close();
        throw e;
    }
}

gr2::OdeSystem* Interface::create_ode_system(std::string text)
{
    text = strip(text);
    std::string ode_name, args_text;
    this->find_command_name(text, ode_name, args_text);
    auto args = this->find_function_arguments(args_text);
    gr2::OdeSystem* ode;

    if (ode_name == "DampedHarmonicOscillator")
    {
        if (args.size() != 2)
            throw std::invalid_argument("invalid number of arguments for DampedHarmonicOscillator");
        ode = new gr2::DampedHarmonicOscillator(std::stold(args[0]), std::stold(args[1]));
    }
    else
    {
        throw std::invalid_argument("spacetime with name " + ode_name + "does not exists");
    }
    return ode;
}

void Interface::solve_ode_system(std::string text)
{
    auto args = find_function_arguments(text);
    int number_of_arguments = 7;
    if (args.size() < number_of_arguments)
        throw std::invalid_argument("too little arguments for solve_ode_system");
    else if (args.size() > number_of_arguments)
        throw std::invalid_argument("too much arguments for solve_ode_system");
    
    gr2::OdeSystem* ode = nullptr;
    std::ofstream file;

    try
    {
        /* code */
        ode = this->create_ode_system(args[0]);
        auto initial_conditions = find_function_arguments(args[1]);
        if (initial_conditions.size() != ode->get_n())
            throw std::invalid_argument("invalid number of initial_value_conditions");
        gr2::real t_start = std::stold(args[2]);
        gr2::real t_end = std::stold(args[3]);
        gr2::real delta_t = std::stold(args[4]);
        std::string method = args[5];
        std::string file_name = args[6];

        // prepare initial conditions
        gr2::real *y_initial = new gr2::real [ode->get_n()];
        for (int i = 0; i < ode->get_n(); i++)
            y_initial[i] = std::stold(initial_conditions[i]);

        // calculation
        gr2::Integrator integrator(*ode, method);
        DataRecord recorder(ode->get_n());
        integrator.add_event(&recorder);
        integrator.integrate(y_initial, t_start, t_end, delta_t);

        // saving data
        file.open(file_name);
        for (auto record:recorder.data)
        {
            for (auto d : record)
                file << d << ";";
            file << std::endl;
        }
    }
    catch(const std::exception& e)
    {
        delete[] ode;
        file.close();
        throw e;
    }
    
}

Interface::Interface():macros(), values()
{

}

bool Interface::command(std::string text)
{
    // skip line
    if (text == "")
        return true;

    // end
    text = strip(text);
    if (text == "end" || text == "END")
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