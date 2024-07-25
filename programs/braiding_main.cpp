#include "braiding.h"

using namespace braiding;
using cgarside::EndLine;
using cgarside::ind_cout;
using cgarside::IndentedOStream;
using cgarside::sint16;

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
        case Option::Quit: {
            ind_cout << "Leaving Braiding." << EndLine(1);
            return 0;
        }
        }
    }
}