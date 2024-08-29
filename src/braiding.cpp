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

#include "braiding/braiding.hpp"
#include "garcide/centralizer.hpp"
#include "garcide/sliding_circuits.hpp"

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

    ind_cout << "For Artin braids, F = ('s' '_'?)? Z | 'D'." << EndLine()
             << "\"s_i\" (optional s_) represents classic Artin generator σ_i."
             << EndLine()
             << "i must therefore be in [1, n[, with n the number of strands."
             << EndLine() << "\"D\" represents the half-twist, Δ_n."
             << EndLine(1) << "For example, with n = 3," << EndLine()
             << "\"1 ^ 1 s2 1\",  \"s_2. s1 . 2\" or \"D\" are three ways to "
                "enter Δ_n."
             << EndLine(1);

#elif BRAIDING_CLASS == 1

    ind_cout << "For Artin dual braids, F = ('a' '_'?)? '(' Z ','? Z ')' | 'D'."
             << EndLine()
             << "\"a_(i,j)\" (optional a_ and \",\") stands for "
                "Birman-Ko-Lee's a_i,j."
             << EndLine() << "i and j must therefore be distinct and in [1, n],"
             << EndLine() << "where n is the number of strands." << EndLine()
             << "\"D\" represents Birman-Ko-Lee's cyclic permutation δ_n."
             << EndLine(1) << "For example, with n = 3," << EndLine()
             << "\"a_(3 2) ^1 (1, 2)\" or \"D\" are two ways to enter δ_n."
             << EndLine(1);

#elif BRAIDING_CLASS == 2

    ind_cout
        << "For B-series dual braids," << EndLine()
        << "F = ('s' '_'?)? '(' Z ','? Z ')' | ('l' '_'?)? Z | 'D'."
        << EndLine()
        << "\"s_(i,j)\" (optional \"s_\" and \",\") stands for what Bessis' "
           "short generator."
        << EndLine()
        << "It behaves (in the Coxeter group) as a double transposition,"
        << EndLine()
        << "Swapping i and j, and i + n and j + n (indexes are mod 2n),"
        << EndLine()
        << "where n is the parameter such that we are working in A(B_n)."
        << EndLine() << "\"l_i\" (optional l_) stands for a long generator."
        << EndLine()
        << "In the Coxeter group, it is the transposition (i i + n)."
        << EndLine()
        << "\"D\" represents the Garside element, which is a cyclic "
           "permutation."
        << EndLine(1) << "For example, with n = 3," << EndLine()
        << "\"(3 2) ^ 1 s(1, 2) l1\" or \"D\" are two ways to enter Δ."
        << EndLine(1);

#elif BRAIDING_CLASS == 3

    ind_cout << "For I-series dual braids, F = ('s' '_'?)? Z | 'D'."
             << EndLine()
             << "\"s_k\" (optional \"s_\") represents the reflection that "
                "sends 1"
             << EndLine() << "on ζ_n^k, where ζ_n = exp(i τ / n)." << EndLine()
             << "\"D\" represents the rotation Δ." << EndLine(1)
             << "For example, with n = 3," << EndLine()
             << "\"s0 2\",  \"1 . s_0\" or \"D ^ 1\" are three ways to enter Δ."
             << EndLine(1);

#elif BRAIDING_CLASS == 4

    ind_cout
        << "For the dual structure on B(e, e, n + 1)," << EndLine()
        << "F = '(' Z ','? Z ')' | Z | 'D'." << EndLine()
        << "\"a_(i, j)\" (optional \"s_\" and \",\") stands for Bessis-Corran's"
        << EndLine() << "symmetric generator a_i,j," << EndLine()
        << "whence |i - j| must be at most n - 1 (taking indexes mod en)."
        << EndLine()
        << "\"a_i\" (optional \"s_\") represents assymetric generator a_i."
        << EndLine() << "\"D\" represents the Garside element." << EndLine(1)
        << "For example, with n = 3, e = 3," << EndLine()
        << "\"a(10 9) ^ 1 (9, 8) 8\" or \"D\" are two ways to enter Δ."
        << EndLine(1);

#elif BRAIDING_CLASS == 5

    ind_cout << "For the semi-classic structure on B(e, e, n)," << EndLine()
             << "F = 's' '_'? Z | 't' '_'? Z | 'D'." << EndLine()
             << "\"s_i\" (optional \"_\") stands for Corran-Picantin's s_i,"
             << EndLine() << "so i must belong to [3, n]." << EndLine()
             << "\"t_i\" (optional \"_\") stands for Corran-Picantin's t_i."
             << EndLine()
             << "'D' represents the Δ element for the Garside structure."
             << EndLine(1) << "For example, with n = 3, e = 3," << EndLine()
             << " \"t1 ^ 1 t_0 s_3 t1 t0 s3\" or \"D\" are two ways to enter Δ."
             << EndLine(1);

#elif BRAIDING_CLASS == 6

    ind_cout
        << "For Artin braids, F = Z | 'D'." << EndLine()
        << "\"e_i\" (optional \"e_\") represents base vector e_i." << EndLine()
        << "i must therefore be in [0, n[, with n the dimension." << EndLine()
        << "\"D\" represents vector Δ = (1, ..., 1)." << EndLine(1)
        << "For example, with n = 3," << EndLine()
        << "\"1 ^ 1 0 e_2\",  \"0. e1 . 2\" or \"D\" are three ways to enter Δ."
        << EndLine(1);

#else

    ind_cout << "Enter a braid (no more details for that group)." << EndLine(1);

#endif
}

void explain_garside_structure() {

#if BRAIDING_CLASS == 0

    ind_cout << "In the classic Garside structure of braid group B_n there is a"
             << EndLine()
             << "one-to-one mapping between canonical factors and permutations."
             << EndLine(1)
             << "The atoms are the Artin generators σ_i, i ∈ [1, n - 1]."
             << EndLine()
             << "As a permutation, σ_i is the transposition (i i+1), and it is"
             << EndLine() << "printed as si." << EndLine(1)
             << "Factors are printed as words in the σ_i." << EndLine(1)
             << "The Garside element Δ_n corresponds to the permutation that"
             << EndLine()
             << "sends i to n - i + 1 (Garside's half-twist). It has length"
             << EndLine() << "n (n - 1) / 2 as a product in the atoms."
             << EndLine(1);

#elif BRAIDING_CLASS == 1

    ind_cout
        << "In the dual Garside structure of braid group B_n, canonical"
        << EndLine()
        << "factors correspond to a subset of the n-th symmetric group."
        << EndLine()
        << "These are the permutations whose disjoint cycle decomposition"
        << EndLine() << "is composed of decreasing cycles, that give rise to a"
        << EndLine() << "non-crossing partition on [1,n]." << EndLine(1)
        << "Thus there is a one-to-one correspondance between factors and"
        << EndLine()
        << "non-crossing partitions. The lattice structure (on both sides)"
        << EndLine() << "is the one induced by that mapping." << EndLine(1)
        << "The atoms are the Birman-Ko-Lee generators a_i,j, where"
        << EndLine() << "(i, j) ∈ [1, n] and i > j." << EndLine()
        << "As a permutation, a_i,j is the transposition (i j), and it is"
        << EndLine()
        << "printed as a(i, j). The corresponding partition is the one"
        << EndLine() << "whose cells are {i, j}, and the {k}, for k != i, j."
        << EndLine(1)
        << "Factors are printed as words in the a_i,j, in a way that makes"
        << EndLine()
        << "it easy to parse the cycle decomposition: for instance, cycle"
        << EndLine() << "(5 2 1) would be printed as a(5, 2) a(2, 1)."
        << EndLine(1)
        << "The Garside element δ_n corresponds to the permutation that"
        << EndLine()
        << "sends i to i + 1 (mod n). It has length n - 1 as a product in"
        << EndLine()
        << "the atoms, and corresponds to the partition with only one cell."
        << EndLine(1)
        << "See Birman, Ko, Lee, A New Approach to the Word and Conjugacy"
        << EndLine() << "Problems in the Braid Groups, 1997, arXiv:math/9712211"
        << EndLine() << "[math.GT]." << EndLine(1);

#elif BRAIDING_CLASS == 2

    ind_cout
        << "There is a morphism from B-series Artin Group A_n(B) to B_2n"
        << EndLine()
        << "that induces an injective Garside group morphism between their"
        << EndLine() << "dual Garside structures." << EndLine(1)
        << "Factors correspond to the permutations whose disjoint cycle"
        << EndLine()
        << "decomposition is composed of decreasing cycles, that give rise"
        << EndLine()
        << "to a non-crossing partition on [1,2n] that is stable by the"
        << EndLine() << "permutation that sends i to i + n (mod 2n)."
        << EndLine(1)
        << "Thus there is a one-to-one correspondance between factors and"
        << EndLine()
        << "those non-crossing partitions. The lattice structure (on both"
        << EndLine() << "sides) is the one induced by that mapping."
        << EndLine(1)
        << "The atoms are the Bessis short and long generators, s_i,j and"
        << EndLine() << "l_i, with i, j ∈ [1, 2n], and i and j distinct mod n."
        << EndLine(1) << "As a permutation, s_i,j is the double transposition"
        << EndLine() << "(i j) (i+n j+n), and it is printed as s(i, j)."
        << EndLine() << "The corresponding partition is the one whose cells are"
        << EndLine()
        << "{i, j}, {i+n, j+n} and the {k}, for k != i, j, i+n, j+n."
        << EndLine(1)
        << "l_i is the transposition (i i+n), and is printed as li."
        << EndLine() << "The corresponding partition is the one" << EndLine()
        << "whose cells are {i, i+n} and the {k}, for k != i, i+n" << EndLine(1)
        << "Factors are printed as words in the s_i,j and l_k, in a way that"
        << EndLine()
        << "makes it easy to parse the cycle decomposition: for instance,"
        << EndLine()
        << "double cycle (5 2 1) (10 7 6) would be printed as s(5, 2)"
        << EndLine() << "s(2, 1)." << EndLine(1)
        << "The Garside element Δ corresponds to the permutation that sends"
        << EndLine()
        << "i to i + 1 (mod 2n). It has length n as a product in the atoms,"
        << EndLine() << "and corresponds to the partition with only one cell."
        << EndLine(1)
        << "See Bessis, The Dual Braid Monoid, 2001, arXiv:math/0101158"
        << EndLine() << "[math.GR]." << EndLine(1);

#elif BRAIDING_CLASS == 3

    ind_cout
        << "Canonical factors for the dual Garside structure of I-series"
        << EndLine()
        << "Artin group A_n(I) are in bijection with a subset of dihedral"
        << EndLine()
        << "group D_2n. More specifically, the identity, the reflections"
        << EndLine()
        << "(which are the atoms), and a rotation (the Garside element)."
        << EndLine(1)
        << "The divisibility lattice has height 2, thus it is particularily"
        << EndLine()
        << "trivial: the identity is the min, the rotattion the max, and"
        << EndLine() << "the reflexions make up the middle layer." << EndLine(1)
        << "The reflexion that sends 0 to exp(k/n iτ) is printed as sk."
        << EndLine(1);

#elif BRAIDING_CLASS == 6

    ind_cout
        << "In the Garside structure for Z^n, canonical factors are the"
        << EndLine() << "vectors whose coordinates are all in {0, 1}."
        << EndLine(1)
        << "The divisibility order is the natural componentwise order."
        << EndLine(1)
        << "The Garside element is vector (1, ..., 1), and the atoms are the"
        << EndLine() << "base vectors, printed as ei." << EndLine(1);

#else

    ind_cout << "No more details." << EndLine(1);

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

    ind_cout << "Enter the group parameter (no more details for that group)."
             << EndLine(1);

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
       << "│││││       Copyright (C) 2004 Juan González-Meneses.      │││││"
       << EndLine()
       << "├───┴──┬────────────────────────────────────────────────┬──┴───┤"
       << EndLine()
       << "││││││││   Braiding comes with ABSOLUTELY NO WARRANTY;  ││││││││"
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

#elif BRAIDING_CLASS == 6

    os << "Using the Garside structure for euclidean lattice Z^n.";

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
        ind_cout << "Choose an option (? for help, gar for a description of the"
                 << EndLine() << "Garside structure):" << EndLine(1) << ">>> ";
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
        } else if (std::regex_match(
                       str, std::regex{"[\\s\\t]*[gG][aA][rR][\\s\\t]*"})) {
            return Option::Garside;
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
    ind_cout << EndLine() << "Their left lcm is:" << EndLine(1)
             << b.left_join(c) << EndLine(1);
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