#include <iostream>
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
            test = intf.command(text);
        }
        catch(const std::exception& e)
        {
            std::cerr << e.what() << '\n';
        }
        
    }
}