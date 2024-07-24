#include "braiding.h"
#include "garsidesc.h"
#include "garsideuss.h"

namespace braiding {

void read_braid(Braid &b, std::istream &is) {
    std::string str;
    std::getline(is, str);
    if (std::regex_match(str, std::regex{"[\\s\\t]*\\?[\\s\\t]*"})) {
        throw HelpAskedFor();
    }
    if (std::regex_match(str, std::regex{"[\\s\\t]*[qQ][\\s\\t]*"})) {
        throw InterruptAskedFor();
    }
    b.OfString(str);
}

Braid::ParameterType read_braid_parameter(std::istream &is) {
    std::string str;
    std::getline(is, str);
    if (std::regex_match(str, std::regex{"[\\s\\t]*\\?[\\s\\t]*"})) {
        throw HelpAskedFor();
    }
    if (std::regex_match(str, std::regex{"[\\s\\t]*[qQ][\\s\\t]*"})) {
        throw InterruptAskedFor();
    }
    return Braid::parameter_of_string(str);
}

void explain_braid_input() {
    // If for some reason you are using a braid implementation
    // that is not a cgarside::Braid<F>, then you may want to add additionnal
    // processor directives in there.
    ind_cout
        << "A braid is entered as a sequence of factors raised to powers."
        << EndLine()
        << "This sequence is interpreted as the product of these factors."
        << EndLine() << "It does not need to be weighted in any manner."
        << EndLine(1)
        << "It should match regexp (W | '.')* (F (W ^ W Z)? (W | '.')*)*."
        << EndLine() << "Where W = (' ' | '\\t')* matches whitespaces,"
        << EndLine()
        << "Z = Z = '-'? (['1' - '9'] ['0' - '9']* | '0') matches integers,"
        << EndLine() << "And F matches (a subset) of factors." << EndLine(1);
#if BRAIDING_CLASS == 0

    ind_cout
        << "For Artin braids, factors must match regexp Z | 'D'." << EndLine()
        << "An integer i represents classic Artin generator σ_i." << EndLine()
        << "It must therefore be in [1, n[, with n the number of strands."
        << EndLine() << "'D' represents the half-twist, Δ." << EndLine(1)
        << "For example, with n = 3," << EndLine()
        << " \"1 ^ 1 2 1\",  \"2. 1 . 2\" or \"D\" Are three ways to enter Δ."
        << EndLine(1);

#endif
}

void explain_braid_parameter_input() {
#if BRAIDING_CLASS == 0

    ind_cout << "Enter the number of strands." << EndLine(1);

#endif
}

void prompt_braid(Braid &b) {
    ind_cout << "Enter a braid (? for help, q to abort): " << EndLine(1);
    while (true) {
        ind_cout << ">>> ";
        try {
            read_braid(b);
            ind_cout << EndLine();
            return;
        } catch (cgarside::InvalidStringError inval) {
            ind_cout << EndLine() << "This is not a valid braid!" << EndLine(1)
                     << inval.error_source << EndLine(1)
                     << "Please try again (? for help, q to abort):"
                     << EndLine(1);
        } catch (HelpAskedFor _) {
            ind_cout << EndLine();
            explain_braid_input();
            ind_cout << "Please try again (? for help, q to abort):"
                     << EndLine(1);
        }
    }
}

Braid::ParameterType prompt_braid_parameter() {
    ind_cout << "Enter the parameter (? for help, q to abort): " << EndLine(1);
    while (true) {
        ind_cout << ">>> ";
        try {
            Braid::ParameterType p = read_braid_parameter();
            ind_cout << EndLine();
            return p;
        } catch (cgarside::InvalidStringError inval) {
            ind_cout << EndLine() << "This is not a valid braid parameter!"
                     << EndLine(1) << inval.error_source << EndLine(1)
                     << "Please try again (? for help, q to abort):"
                     << EndLine(1);
        } catch (HelpAskedFor _) {
            ind_cout << EndLine();
            explain_braid_parameter_input();
            ind_cout << "Please try again (? for help, q to abort):"
                     << EndLine(1);
        }
    }
}

void print_header(IndentedOStream &os) {
    os << "┌──────────────────┬────────────────────────┬──────────────────┐"
       << EndLine()
       << "│──────────────────│  This is Braiding 2.0  │──────────────────│"
       << EndLine()
       << "├──────┬───────────┴────────────────────────┴───────────┬──────┤"
       << EndLine()
       << "││││││││           Copyright and authors here           ││││││││"
       << EndLine()
       << "││││││││   Braiding comes with ABSOLUTELY NO WARRANTY   ││││││││"
       << EndLine()
       << "││││││││             This is free software              ││││││││"
       << EndLine()
       << "││││││││   See GNU General Public License in GPL.txt    ││││││││"
       << EndLine()
       << "├──────┴──┬──────────────────────────────────────────┬──┴──────┤"
       << EndLine()
       << "│││││││││││     _____________________ __________     │││││││││││"
       << EndLine()
       << "│││││││││││     _____________ _______/_ ________     │││││││││││"
       << EndLine()
       << "│││││││││││     _______ _____/_ _______/_ ______     │││││││││││"
       << EndLine()
       << "│││││││││││     ___ ___/_ _____/_ _______/_ ____     │││││││││││"
       << EndLine()
       << "│││││││││││     _ _/_ ___/_ _____/_ _______/_ __     │││││││││││"
       << EndLine()
       << "│││││││││││     _/___/_____/_______/_________/__     │││││││││││"
       << EndLine()
       << "│││││││││││                                          │││││││││││"
       << EndLine()
       << "└─────────┴──────────────────────────────────────────┴─────────┘"
       << EndLine(1);
}

void print_options(IndentedOStream &os) {
#if BRAIDING_CLASS == 0
    os << "Using Artin Braids." << EndLine(1);
#endif
    os << "l:      Left Normal Form        r:      Right Normal Form       "
       << EndLine(1)
       << "^l:     Left GCD                ^r:     Right GCD               "
       << EndLine(1)
       << "vl:     Left LCM                vr:     Right LCM               "
       << EndLine(1)
       << "sss:    Super Summit Set        uss:    Ultra Summit Set        "
       << EndLine(1)
       << "scs:    Sliding Circuits Set    c:      Conjugacy Test          "
       << EndLine(1)
#if BRAIDING_CLASS == 0
       << "t:      Thurston Type           "
#endif
       << "q:      Quit" << EndLine(1);
}

Option prompt_option() {
    while (true) {
        print_line();
        ind_cout << EndLine();
        ind_cout << "Choose an option (? for help):" << EndLine(1) << ">>> ";
        std::string str;
        std::getline(std::cin, str);
        if (std::regex_match(str, std::regex{"[\\s\\t]*\\?[\\s\\t]*"})) {
            ind_cout << EndLine();
            print_options();
        } else if (std::regex_match(str,
                                    std::regex{"[\\s\\t]*[lL][\\s\\t]*"})) {
            return Option::LCF;
        } else if (std::regex_match(str,
                                    std::regex{"[\\s\\t]*[rR][\\s\\t]*"})) {
            return Option::RCF;
        } else if (std::regex_match(str,
                                    std::regex{"[\\s\\t]*\\^[lL][\\s\\t]*"})) {
            return Option::LGCD;
        } else if (std::regex_match(str,
                                    std::regex{"[\\s\\t]*\\^[rR][\\s\\t]*"})) {
            return Option::RGCD;
        } else if (std::regex_match(str,
                                    std::regex{"[\\s\\t]*v[lL][\\s\\t]*"})) {
            return Option::LLCM;
        } else if (std::regex_match(str,
                                    std::regex{"[\\s\\t]*v[rR][\\s\\t]*"})) {
            return Option::RLCM;
        } else if (std::regex_match(
                       str, std::regex{"[\\s\\t]*[sS][sS][sS][\\s\\t]*"})) {
            return Option::SSS;
        } else if (std::regex_match(
                       str, std::regex{"[\\s\\t]*[uU][sS][sS][\\s\\t]*"})) {
            return Option::USS;
        } else if (std::regex_match(
                       str, std::regex{"[\\s\\t]*[sS][cC][sS][\\s\\t]*"})) {
            return Option::SCS;
        } else if (std::regex_match(str,
                                    std::regex{"[\\s\\t]*[cC][\\s\\t]*"})) {
            return Option::Conjugacy;
        } else if (std::regex_match(str,
                                    std::regex{"[\\s\\t]*[qQ][\\s\\t]*"})) {
            return Option::Quit;
#if BRAIDING_CLASS == 0
        } else if (std::regex_match(str,
                                    std::regex{"[\\s\\t]*[tT][\\s\\t]*"})) {
            return Option::ThurstonType;
#endif
        } else {
            ind_cout << EndLine() << "Not a valid option!" << EndLine(1);
        }
    }
}

void print_line(IndentedOStream &os) {
    os << "────────────────────────────────────────────────────────────────"
       << EndLine();
}

void lcf_case() {
    Braid b(prompt_braid_parameter());
    prompt_braid(b);
    ind_cout << EndLine() << "Its left normal form is:" << EndLine(1) << b
             << EndLine(1);
}

void rcf_case() {
    Braid b(prompt_braid_parameter());
    prompt_braid(b);
    b.MakeRCFFromLCF();
    ind_cout << EndLine() << "Its right normal form is:" << EndLine(1);
    b.print_rcf();
    ind_cout << EndLine(1);
}

void left_gcd_case() {
    Braid::ParameterType p = prompt_braid_parameter();
    Braid b(p), c(p);
    prompt_braid(b);
    prompt_braid(c);
    ind_cout << EndLine() << "Their left gcd is:" << EndLine(1) << (b ^ c)
             << EndLine(1);
}

void right_gcd_case() {
    Braid::ParameterType p = prompt_braid_parameter();
    Braid b(p), c(p);
    prompt_braid(b);
    prompt_braid(c);
    ind_cout << EndLine() << "Their right gcd is:" << EndLine(1)
             << b.RightMeet(c) << EndLine(1);
}

void left_lcm_case() {
    Braid::ParameterType p = prompt_braid_parameter();
    Braid b(p), c(p);
    prompt_braid(b);
    prompt_braid(c);
    ind_cout << EndLine() << "Their left lcm is:" << EndLine(1) << b.LeftJoin(c)
             << EndLine(1);
}

void right_lcm_case() {
    Braid::ParameterType p = prompt_braid_parameter();
    Braid b(p), c(p);
    prompt_braid(b);
    prompt_braid(c);
    ind_cout << EndLine() << "Their right lcm is:" << EndLine(1)
             << b.RightJoin(c) << EndLine(1);
}

void sss_case() {
    Braid b(prompt_braid_parameter());
    prompt_braid(b);
    ind_cout << EndLine() << "Its supper summit set is:" << EndLine(1);
    cgarside::super_summit::SSS(b).print(ind_cout);
    ind_cout << EndLine(1);
}

void uss_case() {
    Braid b(prompt_braid_parameter());
    prompt_braid(b);
    ind_cout << EndLine() << "Its ultra summit set is:" << EndLine(1);
    cgarside::ultra_summit::USS(b).print(ind_cout);
    ind_cout << EndLine(1);
}

void scs_case() {
    Braid b(prompt_braid_parameter());
    prompt_braid(b);
    ind_cout << EndLine() << "Its sliding circuit set is:" << EndLine(1);
    cgarside::sliding_circuit::SCS(b).print(ind_cout);
    ind_cout << EndLine(1);
}

void conjugacy_case() {
    Braid::ParameterType p = prompt_braid_parameter();
    Braid b(p), c(p), conj(p);
    prompt_braid(b);
    prompt_braid(c);
    if (cgarside::sliding_circuit::AreConjugate(b, c, conj)) {
        ind_cout << EndLine() << "They are conjugates." << EndLine()
                 << "A conjugating element is:" << EndLine() << conj
                 << EndLine(1);
    } else {
        ind_cout << EndLine() << "They are not conjugates.";
    }
}

#if BRAIDING_CLASS == 0

void thurston_type_case() {
    Braid b(prompt_braid_parameter());
    prompt_braid(b);
    ind_cout << EndLine() << "Its Thurston type is " << cgarside::artin::thurston_type(b) << "." << EndLine(1);
}

#endif

} // namespace braiding