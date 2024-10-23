#include <string>
#include <vector>

class Interface
{
protected:
    // macros
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

    std::string find_function_arguments(std::string text);

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

    bool try_apply_function(std::string text);
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