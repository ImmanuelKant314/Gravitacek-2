#include "interface/interface.hpp"
#include <stdexcept>
#include <iostream>
#include <algorithm>

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

bool Interface::macro_name_char(char c)
{
    if (('A' <= c && 'Z' >= c) || c == '_')
        return true;
    return false;
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

std::string Interface::find_function_arguments(std::string text)
{
    int counter = 0;
    bool start = false;
    int start_index, end_index;
    for (int i = 0; i < text.length(); i++)
    {
        if (text[i]=='(')
            counter ++;
        else if (text[i]==')')
            counter--;

        if (counter == 1 && !start)
        {
            start = true;
            start_index=i;
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
    }
    return text.substr(start_index, end_index-start_index);
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
        else if (!macro_name_char(text[i]))
            throw std::invalid_argument("invalid macro name");
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

bool Interface::try_apply_function(std::string text)
{
    return false;
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