/**
 * @file braiding.cpp
 * @author Matteo Wei (matteo.wei@ens.psl.eu)
 * @brief Implementation file for most of braiding.
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

#include "braiding/braiding.h"
#include "garcide/centralizer.h"
#include "garcide/sliding_circuits.h"

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
    b.of_string(str);
}

Braid::Parameter read_braid_parameter(std::istream &is) {
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
        << "For Artin braids, F = Z | 'D'." << EndLine()
        << "The integer \"i\" represents classic Artin generator σ_i."
        << EndLine()
        << "It must therefore be in [1, n[, with n the number of strands."
        << EndLine() << "\"D\" represents the half-twist, Δ." << EndLine(1)
        << "For example, with n = 3," << EndLine()
        << "\"1 ^ 1 2 1\",  \"2. 1 . 2\" or \"D\" are three ways to enter Δ."
        << EndLine(1);

#elif BRAIDING_CLASS == 1

    ind_cout << "For Artin dual braids, F = '(' Z ','? Z ')' | 'D'."
             << EndLine()
             << "The couple of integers \"(i, j)\" (optional comma)"
             << EndLine() << "stands for Birman-Ko-Lee's a_(i,j)." << EndLine()
             << "i and j must therefore be distinct and in [1, n]," << EndLine()
             << "where n is the number of strands." << EndLine()
             << "\"D\" represents Birman-Ko-Lee's cyclic permutation δ."
             << EndLine(1) << "For example, with n = 3," << EndLine()
             << "\"(3 2) ^1 (1, 2)\" or \"D\" are two ways to enter δ."
             << EndLine(1);

#elif BRAIDING_CLASS == 2

    ind_cout << "For B-series dual braids, F = '(' Z ','? Z ')' | Z | 'D'."
             << EndLine()
             << "The couple of integers \"(i, j)\" (with an optional comma)"
             << EndLine() << "stands for what Bessis calls a short generator."
             << EndLine()
             << "It behaves (in the Coxeter group) as a double transposition,"
             << EndLine()
             << "Swapping i and j, and i + n and j + n (indexes are mod 2n),"
             << EndLine()
             << "where n is the parameter such that we are working in A(B_n)."
             << EndLine() << "The integer i stands for a long generator."
             << EndLine()
             << "In the Coxeter group, it is the transposition (i i + n)."
             << EndLine()
             << "\"D\" represents the Δ element for the Garside structure,"
             << EndLine() << "which is a cyclic permutation." << EndLine(1)
             << "For example, with n = 3," << EndLine()
             << "\"(3 2) ^ 1 (1, 2) 1\" or \"D\" are two ways to enter Δ."
             << EndLine(1);

#elif BRAIDING_CLASS == 3

    ind_cout
        << "For I-series dual braids, F = ('s' '_'?)? Z | 'D'." << EndLine()
        << "\"s_k\" (whith optional s_) represents the reflection that sends 1"
        << EndLine() << "on ζ^k, where ζ = exp(i τ / n)." << EndLine()
        << "\"D\" represents the rotation Δ." << EndLine(1)
        << "For example, with n = 3," << EndLine()
        << "\"s0 2\",  \"1 . s_0\" or \"D ^ 1\" are three ways to enter Δ."
        << EndLine(1);

#elif BRAIDING_CLASS == 4

    ind_cout << "For the dual structure on B(e, e, n + 1)," << EndLine()
             << "F = '(' Z ','? Z ')' | Z | 'D'." << EndLine()
             << "The couple of integers \"(i, j)\" (with an optional comma)"
             << EndLine()
             << "stands for Bessis-Corran's symmetric generator a_i,j,"
             << EndLine()
             << "whence |i - j| must be at most n - 1 (taking indexes mod en)."
             << EndLine()
             << "The integer \"i\" stands for assymetric generator a_i."
             << EndLine()
             << "\"D\" represents the Δ element for the Garside structure."
             << EndLine(1) << "For example, with n = 3, e = 3," << EndLine()
             << "\"(10 9) ^ 1 (9, 8) 8\" or \"D\" are two ways to enter Δ."
             << EndLine(1);

#elif BRAIDING_CLASS == 5

    ind_cout
        << "For the semi-classic structure on B(e, e, n)," << EndLine()
        << "F = 's' '_'? Z | 't' '_'? Z | 'D'." << EndLine()
        << "\"s_i\" (optional underscore) stands for Corran-Picantin's s_i,"
        << EndLine() << "so i must belong to [3, n]." << EndLine()
        << "\"t_i\" (optional underscore) stands for Corran-Picantin's t_i."
        << EndLine()
        << "'D' represents the Δ element for the Garside structure."
        << EndLine(1) << "For example, with n = 3, e = 3," << EndLine()
        << " \"t1 ^ 1 t_0 s_3 t1 t0 s3\" or \"D\" are two ways to enter Δ."
        << EndLine(1);

#elif BRAIDING_CLASS == 6

    ind_cout
        << "For Artin braids, F = Z | 'D'." << EndLine()
        << "\"e_i\" (optional \"e_\") represents base vector e_i."
        << EndLine()
        << "i must therefore be in [0, n[, with n the dimension."
        << EndLine() << "\"D\" represents vector Δ = (1, ..., 1)." << EndLine(1)
        << "For example, with n = 3," << EndLine()
        << "\"1 ^ 1 0 2\",  \"0. 1 . 2\" or \"D\" are three ways to enter Δ."
        << EndLine(1);

#else

    ind_cout << "Enter a braid (no more details for that group)." << EndLine(1);

#endif
}

void explain_braid_parameter_input() {
#if (BRAIDING_CLASS == 0 || BRAIDING_CLASS == 1)

    ind_cout << "Enter the number of strands (an integer)." << EndLine(1);

#elif (BRAIDING_CLASS == 2 || BRAIDING_CLASS == 3)

    ind_cout << "Enter the group parameter (an integer)." << EndLine(1);

#elif (BRAIDING_CLASS == 4)

    ind_cout << "Enter a tuple '(' Z ','? Z ')' of integers." << EndLine()
             << "((e, n) for B(e, e, n + 1).)" << EndLine(1);

#elif (BRAIDING_CLASS == 5)

    ind_cout << "Enter a tuple '(' Z ','? Z ')' of integers." << EndLine()
             << "((e, n) for B(e, e, n).)" << EndLine(1);

#elif (BRAIDING_CLASS == 6)

    ind_cout << "Enter the dimension (an integer)." << EndLine(1);

#else

    ind_cout << "Enter the group parameter (no more details for that group)." << EndLine(1);

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
        } catch (garcide::InvalidStringError inval) {
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

Braid::Parameter prompt_braid_parameter() {
    ind_cout << "Enter the parameter (? for help, q to abort): " << EndLine(1);
    while (true) {
        ind_cout << ">>> ";
        try {
            Braid::Parameter p = read_braid_parameter();
            ind_cout << EndLine();
            return p;
        } catch (garcide::InvalidStringError inval) {
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
       << "│──────────────────│    This is Braiding    │──────────────────│"
       << EndLine()
       << "│──────────────────│      Pre-release       │──────────────────│"
       << EndLine()
       << "├───┬──────────────┴────────────────────────┴──────────────┬───┤"
       << EndLine()
       << "│││││            Copyright (C) 2024 Matteo Wei.            │││││"
       << EndLine()
       << "│││││                   Based on Braiding                  │││││"
       << EndLine()
       << "│││││       Copyright (C) 2004 Juan Gonzalez-Meneses.      │││││"
       << EndLine()
       << "├───┴──┬────────────────────────────────────────────────┬──┴───┤"
       << EndLine()
       << "││││││││  Braiding comes with ABSOLUTELY NO WARRANTY; ││││││││"
       << EndLine()
       << "││││││││   this is free software, and you are welcome   ││││││││"
       << EndLine()
       << "││││││││  to redistribute it under certain conditions.  ││││││││"
       << EndLine()
       << "││││││││  See GNU General Public License in LICENCE.md. ││││││││"
       << EndLine()
       << "├──────┴──┬──────────────────────────────────────────┬──┴──────┤"
#if BRAIDING_CLASS == 0
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
#elif BRAIDING_CLASS == 1
       << EndLine()
       << "│││││││││││               _ __________               │││││││││││"
       << EndLine()
       << "│││││││││││               _/_ ________               │││││││││││"
       << EndLine()
       << "│││││││││││               ___/_ ______               │││││││││││"
       << EndLine()
       << "│││││││││││               _____/_ ____               │││││││││││"
       << EndLine()
       << "│││││││││││               _______/_ __               │││││││││││"
       << EndLine()
       << "│││││││││││               _________/__               │││││││││││"
       << EndLine()
       << "│││││││││││                                          │││││││││││"
#elif BRAIDING_CLASS == 2
       << EndLine()
       << "│││││││││││                 +------+                 │││││││││││"
       << EndLine()
       << "│││││││││││                / \\      \\                │││││││││││"
       << EndLine()
       << "│││││││││││               /   \\      \\               │││││││││││"
       << EndLine()
       << "│││││││││││              +     +------+              │││││││││││"
       << EndLine()
       << "│││││││││││               \\   /      /               │││││││││││"
       << EndLine()
       << "│││││││││││                \\ /      /                │││││││││││"
       << EndLine()
       << "│││││││││││                 +------+                 │││││││││││"
#elif BRAIDING_CLASS == 3
       << EndLine()
       << "│││││││││││                 +------+                 │││││││││││"
       << EndLine()
       << "│││││││││││                /        \\                │││││││││││"
       << EndLine()
       << "│││││││││││               /          \\               │││││││││││"
       << EndLine()
       << "│││││││││││              +            +              │││││││││││"
       << EndLine()
       << "│││││││││││               \\          /               │││││││││││"
       << EndLine()
       << "│││││││││││                \\        /                │││││││││││"
       << EndLine()
       << "│││││││││││                 +------+                 │││││││││││"
#elif BRAIDING_CLASS == 4
       << EndLine()
       << "│││││││││││                 <------>                 │││││││││││"
       << EndLine()
       << "│││││││││││                                          │││││││││││"
       << EndLine()
       << "│││││││││││                                          │││││││││││"
       << EndLine()
       << "│││││││││││              <-----∧      ∧              │││││││││││"
       << EndLine()
       << "│││││││││││               \\   /      /               │││││││││││"
       << EndLine()
       << "│││││││││││                \\ /      /                │││││││││││"
       << EndLine()
       << "│││││││││││                 ∨      ∨                 │││││││││││"
#elif BRAIDING_CLASS == 5
       << EndLine()
       << "│││││││││││                _                         │││││││││││"
       << EndLine()
       << "│││││││││││              | j 0 0 0 0 |               │││││││││││"
       << EndLine()
       << "│││││││││││              | 0 j 0 0 0 |               │││││││││││"
       << EndLine()
       << "│││││││││││              | 0 0 j 0 0 |               │││││││││││"
       << EndLine()
       << "│││││││││││              | 0 0 0 j 0 |               │││││││││││"
       << EndLine()
       << "│││││││││││              | 0 0 0 0 j |               │││││││││││"
       << EndLine()
       << "│││││││││││                                          │││││││││││"
#elif (BRAIDING_CLASS == 6)
       << EndLine()
       << "│││││││││││   |    |    |    |    |    |    |    |   │││││││││││"
       << EndLine()
       << "│││││││││││---+----+----+----+----+----+----+----+---│││││││││││"
       << EndLine()
       << "│││││││││││   |    |    |    |    |    |    |    |   │││││││││││"
       << EndLine()
       << "│││││││││││---+----+----+----+----+----+----+----+---│││││││││││"
       << EndLine()
       << "│││││││││││   |    |    |    |    |    |    |    |   │││││││││││"
       << EndLine()
       << "│││││││││││---+----+----+----+----+----+----+----+---│││││││││││"
       << EndLine()
       << "│││││││││││   |    |    |    |    |    |    |    |   │││││││││││"
#else
       << EndLine()
       << "│││││││││││                                          │││││││││││"
       << EndLine()
       << "│││││││││││                                          │││││││││││"
       << EndLine()
       << "│││││││││││                                          │││││││││││"
       << EndLine()
       << "│││││││││││                                          │││││││││││"
       << EndLine()
       << "│││││││││││                                          │││││││││││"
       << EndLine()
       << "│││││││││││                                          │││││││││││"
       << EndLine()
       << "│││││││││││                                          │││││││││││"
#endif
       << EndLine()
       << "└─────────┴──────────────────────────────────────────┴─────────┘"
       << EndLine(1);
}

void print_options(IndentedOStream &os) {
#if BRAIDING_CLASS == 0
    os << "Using Garside's classic structure for Artin braids." << EndLine(1);
#elif BRAIDING_CLASS == 1
    os << "Using Birman-Ko-Lee's dual structure for Artin braid groups."
       << EndLine(1);
#elif BRAIDING_CLASS == 2
    os << "Using the dual structure for B-series Artin groups." << EndLine(1);
#elif BRAIDING_CLASS == 3
    os << "Using the dual structure for I-series Artin groups." << EndLine(1);
#elif BRAIDING_CLASS == 4
    os << "Using Bessis-Corran's dual structure for B(e,e,n) complex"
       << EndLine() << "reflection braid groups." << EndLine(1);
#elif BRAIDING_CLASS == 5
    os << "Using Corran-Picantin's semi-classic structure for B(e,e,n)"
       << EndLine() << "complex reflection braid groups." << EndLine(1);
#endif
    os << "l:      Left Normal Form        r:      Right Normal Form       "
       << EndLine(1)
       << "^l:     Left GCD                ^r:     Right GCD               "
       << EndLine(1)
       << "vl:     Left LCM                vr:     Right LCM               "
       << EndLine(1)
       << "sss:    Super Summit Set        uss:    Ultra Summit Set        "
       << EndLine(1)
       << "scs:    Sliding Circuits Set    ctr:    Centralizer             "
       << EndLine(1) << "c:      Conjugacy Test          "
#if BRAIDING_CLASS == 0
       << "t:      Thurston Type           " << EndLine(1)
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
        } else if (std::regex_match(
                       str, std::regex{"[\\s\\t]*[cC][tT][rR][\\s\\t]*"})) {
            return Option::Centralizer;
        } else if (std::regex_match(str,
                                    std::regex{"[\\s\\t]*[cC][\\s\\t]*"})) {
            return Option::Conjugacy;
        } else if (std::regex_match(str,
                                    std::regex{"[\\s\\t]*[qQ][\\s\\t]*"})) {
            return Option::Quit;
        } else if (std::regex_match(str,
                                    std::regex{"[\\s\\t]*[hH][\\s\\t]*"})) {
            return Option::Header;
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
    b.lcf_to_rcf();
    ind_cout << EndLine() << "Its right normal form is:" << EndLine(1);
    b.print_rcf();
    ind_cout << EndLine(1);
}

void left_gcd_case() {
    Braid::Parameter p = prompt_braid_parameter();
    Braid b(p), c(p);
    prompt_braid(b);
    prompt_braid(c);
    ind_cout << EndLine() << "Their left gcd is:" << EndLine(1) << (b ^ c)
             << EndLine(1);
}

void right_gcd_case() {
    Braid::Parameter p = prompt_braid_parameter();
    Braid b(p), c(p);
    prompt_braid(b);
    prompt_braid(c);
    ind_cout << EndLine() << "Their right gcd is:" << EndLine(1)
             << b.right_meet(c) << EndLine(1);
}

void left_lcm_case() {
    Braid::Parameter p = prompt_braid_parameter();
    Braid b(p), c(p);
    prompt_braid(b);
    prompt_braid(c);
    ind_cout << EndLine() << "Their left lcm is:" << EndLine(1) << b.left_join(c)
             << EndLine(1);
}

void right_lcm_case() {
    Braid::Parameter p = prompt_braid_parameter();
    Braid b(p), c(p);
    prompt_braid(b);
    prompt_braid(c);
    ind_cout << EndLine() << "Their right lcm is:" << EndLine(1)
             << b.right_join(c) << EndLine(1);
}

void sss_case() {
    Braid b(prompt_braid_parameter());
    prompt_braid(b);
    ind_cout << EndLine();
    garcide::super_summit::super_summit_set(b).print(ind_cout);
    ind_cout << EndLine(1);
}

void uss_case() {
    Braid b(prompt_braid_parameter());
    prompt_braid(b);
    ind_cout << EndLine();
    garcide::ultra_summit::ultra_summit_set(b).print();
}

void scs_case() {
    Braid b(prompt_braid_parameter());
    prompt_braid(b);
    ind_cout << EndLine();
    garcide::sliding_circuits::sliding_circuits_set(b).print();
}

void centralizer_case() {
    Braid b(prompt_braid_parameter());
    prompt_braid(b);
    ind_cout << EndLine();
    garcide::centralizer::centralizer(b).print();
}

void conjugacy_case() {
    Braid::Parameter p = prompt_braid_parameter();
    Braid b(p), c(p), conj(p);
    prompt_braid(b);
    prompt_braid(c);
    if (garcide::sliding_circuits::are_conjugate(b, c, conj)) {
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
    ind_cout << EndLine() << "Its Thurston type is "
             << garcide::artin::thurston_type(b) << "." << EndLine(1);
}

#endif

} // namespace braiding