#include "braiding_option.h"

#ifndef NEW_BRAIDING
#define NEW_BRAIDING

namespace braiding {
using cgarside::EndLine;
using cgarside::ind_cout;
using cgarside::IndentedOStream;
using cgarside::sint16;

/**
 * @brief Exception raised when asking for help.
 *
 * An exception that is raised when asking for help (e.g., to ask what a legal
 * factor is).
 */
struct HelpAskedFor {};

/**
 * @brief Exception raised when asking to return to option selection.
 *
 * An exception that is raised when asking to abort the current task and return
 * to option selection.
 */
struct InterruptAskedFor {};

enum class Option {
    LCF,
    RCF,
    LGCD,
    RGCD,
    LLCM,
    RLCM,
    SSS,
    USS,
    SCS,
    Conjugacy,
    Header,
    Quit,
#if BRAIDING_CLASS == 0
    ThurstonType,
#endif
};

/**
 * @brief Reads a line from `is` as a braid.
 *
 * Reads a line from `is`, parses it as a braid and sets `b` to that braid using
 * `b.OfString`.
 *
 * @param b The braid that is to be set to what is read from `is`.
 * @param is The input stream we read from. Its default value is `std::cin`.
 * @exception `cgarside::InvalidStringError`: may be raised by `OfString`.
 * @exception `HelpAskedFor`: raised when `"?"` is entered.
 * @exception `InterruptAskedFor`: raised when `"q"` or `"Q"` is entered.
 */
void read_braid(Braid &b, std::istream &is = std::cin);

/**
 * @brief Reads a line from `is` as a parameter.
 *
 * Reads a line from `is`, parses it as a parameter and returns it.
 *
 * @param is The input stream we read from. Its default value is `std::cin`.
 * @exception `cgarside::InvalidStringError`: may be raised by
 * `parameter_of_string`.
 * @exception `HelpAskedFor`: raised when `"?"` is entered.
 * @exception `InterruptAskedFor`: raised when `"q"` or `"Q"` is entered.
 */
Braid::ParameterType read_braid_parameter(std::istream &is = std::cin);

/**
 * @brief Prints a text that explains what a correct braid input is.
 *
 * Prints a text that explains what a correct braid input is. Typically used
 * when catching `HelpAskedFor`.
 */
void explain_braid_input();

/**
 * @brief Prints a text that explains what a correct braid parameter input is.
 *
 * Prints a text that explains what a correct braid parameter input is.
 * Typically used when catching `HelpAskedFor`.
 */
void explain_braid_parameter_input();

/**
 * @brief Prompts the user to enter a braid.
 *
 * Prompts the user to enter a braid in standard input,
 * then parses the input and sets `b` to that braid using
 * `read_braid`.
 *
 * Keeps asking for braids until a valid one is entered.
 *
 * @param b The braid that is to be set to what is received.
 * @exception `InterruptAskedFor`: raised when `"q"` or `"Q"` is entered.
 */
void prompt_braid(Braid &b);

/**
 * @brief Prompts the user to enter a braid parameter.
 *
 * Prompts the user to enter a braid parameter in standard input,
 * then parses the input and returns the parameters.
 *
 * Keeps asking for parameters until a valid one is entered.
 *
 * @exception `InterruptAskedFor`: raised when `"q"` or `"Q"` is entered.
 */
Braid::ParameterType prompt_braid_parameter();

/**
 * @brief Prints a line.
 *
 * Prints a line made of 64 `'─'` in `os`.
 *
 * @param os The output stream it is printed in.
 */
void print_line(IndentedOStream &os = ind_cout);

/**
 * @brief Prints Braiding's header.
 *
 * Prints Braiding 2.0's beautiful header (with real braid Unicode art!) in
 * `os`.
 *
 * Go admire it in braiding.cpp.
 *
 * @param os The output stream it is printed in.
 */
void print_header(IndentedOStream &os = ind_cout);

/**
 * @brief Prints Braiding's options.
 *
 * Prints the different options in `os`.
 *
 * @param os The output stream it is printed in.
 */
void print_options(IndentedOStream &os = ind_cout);

/**
 * @brief Prompts the user to enter an option.
 *
 * Prompts the user to enter an option in standard input,
 * then parses the input return the corresponding option.
 *
 * Keep asking for options until a valid one is entered.
 */
Option prompt_option();

/**
 * @brief Handles the lcf case.
 *
 * Prompts the user to enter a braid, put it in LCF, then prints it in standard
 * output.
 *
 * @exception `InterruptAskedFor`: raised when `"q"` or `"Q"` is entered on an
 * input.
 */
void lcf_case();

void rcf_case();

void left_gcd_case();

void right_gcd_case();

void left_lcm_case();

void right_lcm_case();

void sss_case();

void uss_case();

void scs_case();

void conjugacy_case();

#if BRAIDING_CLASS == 0

void thurston_type_case();

#endif

} // namespace braiding

#endif