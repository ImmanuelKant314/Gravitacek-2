#include <iostream>
#include <chrono>

#include "gravitacek2/setup.hpp"
#include "gravitacek2/integrator/odesystem.hpp"
#include "interface/interface.hpp"

int main()
{
    std::cout << "Gravitacek2" << std::endl;
    Interface intf = Interface();
    std::string text;
    bool test = true;

    while(test)
    {
        std::getline(std::cin, text);
        try
        { 
            auto start = std::chrono::high_resolution_clock::now();
            test = intf.command(text);
            auto end = std::chrono::high_resolution_clock::now();
            std::cout << "Time of execution: " << std::chrono::duration_cast<std::chrono::duration<double>>(end-start).count() << " s" << std::endl;
        }
        catch(const std::exception& e)
        {
            std::cerr << e.what() << '\n';
        }
        
    }
}