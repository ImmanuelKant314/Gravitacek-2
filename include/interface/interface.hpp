#pragma once

#include <string>
#include <vector>
#include <memory>

#include "gravitacek2/setup.hpp"
#include "gravitacek2/geomotion/geomotion.hpp"
#include "gravitacek2/geomotion/weyl.hpp"

class Interface
{
protected:
    // ==================== Macros ==================== 
    std::vector<std::string> macros;    //!<vector of macro names
    std::vector<std::string> values;    //!<vector of macro numbers

    // ==================== Help ==================== 
    std::vector<std::string> help_name; //!<vector of help names
    std::vector<std::string> help_text; //!<vector of texts for help

    /**
     * @brief Substitute text using macros.
     * 
     * @param text text for substitution
     * @return substituted string
     */
    std::string substitute(std::string text);

    /**
     * @brief Check validity of name of macro.
     * 
     * @param name potential name for macro
     * @return true if name is valid
     * @return false if name is not valid
     */
    bool macro_name_valid(std::string name);

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

    /**
     * @brief Print all saved macros.
     * 
     */
    void print_all_macros();

    /**
     * @brief Print help for functions.
     * 
     * @param text parameter for help
     */
    void help(std::string text);

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
     * @return false if function was not applied
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
    std::shared_ptr<gr2::GeoMotion> create_spacetime(std::string text);

    /**
     * @brief Create Weyl spacetime with given parameters.
     * 
     * @param text definition of Weyl spacetime
     * @return pointer to spacetime object
     */
    std::shared_ptr<gr2::Weyl> create_weyl_spacetime(std::string text);

    /**
     * @brief Create OdeSystem with given parameters.
     * 
     * Argument should be in form:
     * `ode_system_name(ode_system_parameters)`
     * 
     * @param text definition of OdeSystem
     * @return pointer to OdeSystem object
     */
    std::shared_ptr<gr2::OdeSystem> create_ode_system(std::string text);

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

    /**
     * @brief Solve given differential equation.
     * 
     * Argument should be in form:
     * (ode_system_name(ode_system_params),initial_conditions,t_start,t_end,delta_t,method,file)
     * 
     * @param text arguments for solve_ode
     */
    void solve_ode_system(std::string text);

    /**
     * @brief Calculate values of local expansion for Weyl spacetime.
     * 
     * Argument should be in form:
     * (weyl_spacetime(weyl_spacetimes_params),E,L,(rho_min,rho_max,n_rho),(z_min,z_max,n_z),angles,file)
     * 
     * @param text arguments for local_expansions
     */
    void local_expansions_Weyl(std::string text);
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