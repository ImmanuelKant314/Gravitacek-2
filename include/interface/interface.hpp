#pragma once

#include <string>
#include <vector>

#include "gravitacek2/setup.hpp"
#include "gravitacek2/geomotion/geomotion.hpp"
#include "gravitacek2/geomotion/weyl.hpp"

class Interface
{
protected:
    // ==================== Macros ==================== 
    std::vector<std::string> macros;    //!<vector of macro names
    std::vector<std::string> values;    //!<vector of macro numbers

    /**
     * @brief Substitute text using macros.
     * 
     * @param text text for substitution
     * @return substituted string
     */
    std::string substitute(std::string text);

    /**
     * @brief Check validity of characted in macro name.
     * 
     * @param c characted in macro name
     * @return true if character is valid
     * @return false if character is not valid
     */
    bool macro_name_char(char c);

    /**
     * @brief Delete spaces on the begining and end.
     * 
     * @param text modified text
     * @return text without spaces on on begining and end
     */
    std::string strip(std::string text);

    // finding
    /**
     * @brief Izolate command from text.
     * 
     * @param text analyzed text
     * @param command extracted command
     * @param rest remaining text (without command)
     */
    void find_command_name(std::string text, std::string &command, std::string &rest);

    // operators
    /**
     * @brief Try to apply operators.
     * 
     * @param text analyzed text
     * @return true if operator was applied
     * @return false if operator was not applied
     */
    bool try_apply_operators(std::string text);

    /**
     * @brief Define macro.
     * 
     * @param text macro name with macro value separated by space
     */
    void define_macro(std::string text);

    /**
     * @brief Delete macro.
     * 
     * @param text name of macro
     */
    void delete_macro(std::string text);

    /**
     * @brief Print value of macro.
     * 
     * @param text name of macro
     */
    void print_macro(std::string text);

    // ==================== Functions ==================== 

    /**
     * @brief Find arguments in text.
     * 
     * @param text string with arguments
     * @return vector of arguments
     */
    std::vector<std::string> find_function_arguments(std::string text);

    /**
     * @brief Try applying function.
     * 
     * @param text string with possible function
     * @return true if function was applied
     * @return false  if function was not applied
     */
    bool try_apply_function(std::string text);

    /**
     * @brief Write given arguments.
     * 
     * @param text string with arguments
     */
    void split_args(std::string text);

    /**
     * @brief Create spacetime with given parameters.
     * 
     * @param text definition of spacetime
     * @return pointer to spacetime object
     */
    gr2::GeoMotion* create_spacetime(std::string text);

    /**
     * @brief Create Weyl spacetime with given parameters.
     * 
     * @param text definition of Weyl spacetime
     * @return pointer to spacetime object
     */
    gr2::Weyl* create_weyl_spacetime(std::string text);

    /**
     * @brief Calculate potential with given parameters.
     * 
     * @param text parameters for function
     */
    void draw_potential_1D(std::string text);

    /**
     * @brief Calculate metric function lambda with given parameters.
     * 
     * @param text parameters for function
     */
    void draw_lambda_1D(std::string text);

public:
    Interface();

    /**
     * @brief Execute command
     * 
     * @param text command
     * @return true if application should not stop
     * @return false if application should stop
     */
    bool command(std::string text);
};