/**
 * @file braiding.h
 * @author Matteo Wei (matteo.wei@ens.psl.eu)
 * @brief Header file for most of braiding.
 * @version 0.1
 * @date 2024-07-28
 *
 * @copyright Copyright (C) 2024. Distributed under the GNU General Public
 * License, version 3.
 *
 */

/*
 * GarCide Copyright (C) 2024 Matteo Wei.
 *
 * This program is free software: you can redistribute it and/or modify
 * it under the terms of the GNU General Public License as published by
 * the Free Software Foundation, either version 3 of the License, or
 * (at your option) any later version.
 *
 * This program is distributed in the hope that it will be useful,
 * but WITHOUT ANY WARRANTY; without even the implied warranty of
 * MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE. See the
 * GNU General Public License in LICENSE for more details.
 *
 * You should have received a copy of the GNU General Public License
 * along with this program.  If not, see <http://www.gnu.org/licenses/>.
 */

#include "braiding_option.h"

#ifndef BRAIDING
#define BRAIDING

/**
 * @brief Namespace for everything related to the _Braiding_ executable.
 */
namespace braiding {
using garcide::EndLine;
using garcide::ind_cout;
using garcide::IndentedOStream;
using garcide::i16;

/**
 * @brief Exception raised to ask for help.
 */
struct HelpAskedFor {};

/**
 * @brief Exception raised to return to option selection.
 */
struct InterruptAskedFor {};

/**
 * @brief An `enum` that enumerates the options for _Braiding_'s menu.
 *
 */
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
    Centralizer,
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
 * Then sets `b` to that braid using `b.of_string`.
 *
 * @param b The braid that is to be set to what is read from `is`.
 * @param is The input stream we read from. Its default value is `std::cin`.
 * @exception `garcide::InvalidStringError`: may be raised by `of_string`.
 * @exception `HelpAskedFor`: raised when `"?"` is entered.
 * @exception InterruptAskedFor Raised when `"q"` or `"Q"` is entered.
 */
void read_braid(Braid &b, std::istream &is = std::cin);

/**
 * @brief Reads a line from `is` as a parameter and returns it.
 *
 * @param is The input stream we read from. Its default value is `std::cin`.
 * @exception `garcide::InvalidStringError`: may be raised by
 * `parameter_of_string`.
 * @exception `HelpAskedFor`: raised when `"?"` is entered.
 * @exception InterruptAskedFor Raised when `"q"` or `"Q"` is entered.
 */
Braid::Parameter read_braid_parameter(std::istream &is = std::cin);

/**
 * @brief Prints a text that explains what a correct braid input is.
 *
 * Typically used when catching `HelpAskedFor`.
 */
void explain_braid_input();

/**
 * @brief Prints a text that explains what a correct braid parameter input is.
 *
 * Typically used when catching `HelpAskedFor`.
 */
void explain_braid_parameter_input();

/**
 * @brief Prompts the user to enter a braid.
 *
 * Then parses the input and sets `b` to that braid using
 * `read_braid`.
 *
 * Keeps asking for braids until a valid one is entered.
 *
 * @param b The braid that is to be set to what is received.
 * @exception InterruptAskedFor Raised when `"q"` or `"Q"` is entered.
 */
void prompt_braid(Braid &b);

/**
 * @brief Prompts the user to enter a braid parameter.
 *
 * Then parses the input and returns the parameters.
 *
 * Keeps asking for parameters until a valid one is entered.
 *
  * @exception InterruptAskedFor Raised when `"q"` or `"Q"` is entered.
 */
Braid::Parameter prompt_braid_parameter();

/**
 * @brief Prints a line in output stream `os`.
 *
 * It is made of 64 `'`─`'`.
 *
 * @param os The output stream it is printed in.
 */
void print_line(IndentedOStream &os = ind_cout);

/**
 * @brief Prints _Braiding_'s header.
 *
 * Prints _Braiding_'s beautiful header (with real braid Unicode art!) in
 * `os`.
 *
 * Go admire it in braiding.cpp.
 *
 * @param os The output stream it is printed in.
 */
void print_header(IndentedOStream &os = ind_cout);

/**
 * @brief Prints _Braiding_'s options in output stream `os`.
 *
 * @param os The output stream it is printed in.
 */
void print_options(IndentedOStream &os = ind_cout);

/**
 * @brief Prompts the user to enter an option.
 *
 * Then parses the input return the corresponding option.
 *
 * Keep asking for options until a valid one is entered.
 */
Option prompt_option();

/**
 * @brief Handles the LCF case.
 *
 * Prompts the user to enter a braid, puts it in LCF, then prints it in standard
 * output.
 *
 * @exception InterruptAskedFor Raised when `"q"` or `"Q"` is entered.
 */
void lcf_case();

/**
 * @brief Handles the RCF case.
 *
 * Prompts the user to enter a braid, puts it in RCF, then prints it in standard
 * output.
 *
 * @exception InterruptAskedFor Raised when `"q"` or `"Q"` is entered.
 */
void rcf_case();

/**
 * @brief Handles the left GCD case.
 *
 * Prompts the user to enter two braids, computes their left GCD, then prints it
 * in standard output.
 *
 * @exception InterruptAskedFor Raised when `"q"` or `"Q"` is entered.
 */
void left_gcd_case();

/**
 * @brief Handles the right GCD case.
 *
 * Prompts the user to enter two braids, computes their right GCD, then prints
 * it in standard output.
 *
 * @exception InterruptAskedFor Raised when `"q"` or `"Q"` is entered.
 */
void right_gcd_case();

/**
 * @brief Handles the left LCM case.
 *
 * Prompts the user to enter two braids, computes their left LCM, then prints it
 * in standard output.
 *
 * @exception InterruptAskedFor Raised when `"q"` or `"Q"` is entered.
 */
void left_lcm_case();

/**
 * @brief Handles the right LCM case.
 *
 * Prompts the user to enter two braids, computes their right LCM, then prints
 * it in standard output.
 *
 * @exception InterruptAskedFor Raised when `"q"` or `"Q"` is entered.
 */
void right_lcm_case();

/**
 * @brief Handles the \f$\mathrm{SSS}\f$ case.
 *
 * Prompts the user to enter a braid, computes its super summit set, then prints
 * it in standard output.
 *
 * @exception InterruptAskedFor Raised when `"q"` or `"Q"` is entered.
 */
void sss_case();

/**
 * @brief Handles the \f$\mathrm{USS}\f$ case.
 *
 * Prompts the user to enter a braid, computes its ultra summit set, then prints
 * it in standard output.
 *
 * @exception InterruptAskedFor Raised when `"q"` or `"Q"` is entered.
 */
void uss_case();

/**
 * @brief Handles the \f$\mathrm{SCS}\f$ case.
 *
 * Prompts the user to enter a braid, computes its sliding circuits set, then
 * prints it in standard output.
 *
 * @exception InterruptAskedFor Raised when `"q"` or `"Q"` is entered.
 */
void scs_case();

/**
 * @brief Handles the centralizer case.
 *
 * Prompts the user to enter a braid, computes its centralizer, then prints
 * it in standard output.
 *
 * @exception InterruptAskedFor Raised when `"q"` or `"Q"` is entered.
 */
void centralizer_case();

/**
 * @brief Handles the conjugacy check case.
 *
 * Prompts the user to enter two braids, then checks if they are conjugates, and
 * prints the result. If they are conjugate, also prints a conjugator in
 * standard output.
 *
 * @exception InterruptAskedFor Raised when `"q"` or `"Q"` is entered.
 */
void conjugacy_case();

#if BRAIDING_CLASS == 0

/**
 * @brief Handles the Thurston type case.
 *
 * Prompts the user to enter an Artin braid, then computes its Thurston type and
 * prints it in standard output.
 *
 * Only defined if preprocessor variable `BRAIDING_CLASS` is set to `0`.
 * 
 * @exception InterruptAskedFor Raised when `"q"` or `"Q"` is entered.
 */
void thurston_type_case();

#endif

} // namespace braiding

#endif