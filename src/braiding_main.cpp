/**
 * @file braiding_main.cpp
 * @author Matteo Wei (matteo.wei@ens.psl.eu)
 * @brief Main file for braiding executable.
 * @version 1.0.0
 * @date 2024-08-31
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

#include "braiding/braiding.hpp"

using namespace braiding;

/**
 * @brief Main function.
 * 
 * @return `0`.
 */
int main() {
    ind_cout << EndLine();
    print_header();
    print_options();
    while (true) {
        Option option = prompt_option();
        ind_cout << EndLine();
        switch (option) {
        case Option::LCF: {
            try {
                lcf_case();
            } catch (InterruptAskedFor) {
                ind_cout << EndLine();
            }
            break;
        }
        case Option::RCF: {
            try {
                rcf_case();
            } catch (InterruptAskedFor) {
                ind_cout << EndLine();
            }
            break;
        }
        case Option::LGCD: {
            try {
                left_gcd_case();
            } catch (InterruptAskedFor) {
                ind_cout << EndLine();
            }
            break;
        }
        case Option::RGCD: {
            try {
                right_gcd_case();
            } catch (InterruptAskedFor) {
                ind_cout << EndLine();
            }
            break;
        }
        case Option::LLCM: {
            try {
                left_lcm_case();
            } catch (InterruptAskedFor) {
                ind_cout << EndLine();
            }
            break;
        }
        case Option::RLCM: {
            try {
                right_lcm_case();
            } catch (InterruptAskedFor) {
                ind_cout << EndLine();
            }
            break;
        }
        case Option::SSS: {
            try {
                sss_case();
            } catch (InterruptAskedFor) {
                ind_cout << EndLine();
            }
            break;
        }
        case Option::USS: {
            try {
                uss_case();
            } catch (InterruptAskedFor) {
                ind_cout << EndLine();
            }
            break;
        }
        case Option::SCS: {
            try {
                scs_case();
            } catch (InterruptAskedFor) {
                ind_cout << EndLine();
            }
            break;
        }
        case Option::Centralizer: {
            try {
                centralizer_case();
            } catch (InterruptAskedFor) {
                ind_cout << EndLine();
            }
            break;
        }
        case Option::Conjugacy: {
            try {
                conjugacy_case();
            } catch (InterruptAskedFor) {
                ind_cout << EndLine();
            }
            break;
        }
#if BRAIDING_CLASS == 0
        case Option::ThurstonType: {
            try {
                thurston_type_case();
            } catch (InterruptAskedFor) {
                ind_cout << EndLine();
            }
            break;
        }
#endif
        case Option::Header: {
            print_header();
            break;
        }
        case Option::Garside: {
            explain_garside_structure();
            break;
        }
        case Option::Quit: {
            ind_cout << "Leaving Braiding." << EndLine(1);
            return 0;
        }
        }
    }
}