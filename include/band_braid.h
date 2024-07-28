#ifndef BAND
#define BAND

#include "cgarside.h"

#ifdef USE_CLN

#include <cln/integer.h>

#endif

namespace cgarside::band {

class Underlying {

  protected:
    sint16 PresentationParameter;

    std::vector<sint16> permutation_table;

  public:
    /**
     * @brief Maximum braid index.
     *
     * The greatest index that may be used for braids.
     *
     * It is used because we use `thread_local` objects to avoid some
     * allocations, and their size must be known at compile time.
     *
     * Having too big `thread_local` objects might cause some issue with thread
     * spawning.
     */
    static const sint16 MaxBraidIndex = 256;

    typedef sint16 Parameter;

    static Parameter parameter_of_string(const std::string &str);

    Parameter get_parameter() const;

    sint16 lattice_height() const;

    // Constructor
    Underlying(sint16 n);

    sint16 at(sint16 i) const { return permutation_table[i]; }
    sint16 &at(sint16 i) { return permutation_table[i]; }

    /**
     * @brief Extraction from string.
     *
     * Reads the string `str`, starting at position `pos`, and sets `this` to
     * the corresponding atom.
     *
     * Letting `W = (\s | \t)*` be the language of whitespaces and `Z = -? ([1 -
     * 9] [0 - 9]* | 0)` be the language of integers, accepted strings are those
     * represented by the regular expression `\(W Z W,? W Z W\)`, under the
     * additional hypothesis that the two integers lie in [`1`, `Parameter`] and
     * are not equal.
     *
     * @param str The string to extract from.
     * @param pos The position to start from.
     * @exception `InvalidStringError`: Thrown when there is no subword starting
     * from `pos` that matches `\(W Z W,? W Z W\)`, or if there is one, if
     * either integer does not belong to [`1`, `Parameter`], or both are equal.
     */
    void of_string(const std::string &str, size_t &pos);

    /**
     * @brief Prints internal representation to `os`.
     *
     * Prints the factor's `permutation_table` to `os`, typically for debugging
     * purposes.
     *
     * @param os The output stream it prints to.
     */
    void debug(IndentedOStream &os) const;

    void AssignDCDT(sint16 *x) const;

    void OfDCDT(const sint16 *x);

    /**
     * @brief Prints the factor to `os`.
     *
     * Prints the factor to `os` as a product of atoms.
     *
     * @param os The output stream it prints to.
     */
    void print(IndentedOStream &os) const;

    // Set to the identity element (here the identity).
    void identity();

    // Set to delta.
    void delta();

    Underlying left_meet(const Underlying &b) const;

    Underlying right_meet(const Underlying &b) const;

    // Equality check.
    // We check wether the underlying permutation table are (pointwise) equal.
    bool compare(const Underlying &b) const;

    // Computes the factor corresponding to the inverse permutation.
    // Used to simplify complement operation.
    Underlying Inverse() const;

    // product under the hypothesis that it is still simple.
    Underlying product(const Underlying &b) const;

    // Under the assumption a <= b, a.left_complement(b) computes
    // The factor c such that ac = b.
    Underlying left_complement(const Underlying &b) const;

    Underlying right_complement(const Underlying &b) const;

    // Generate a random factor.
    void randomize();

    // List of atoms.
    std::vector<Underlying> atoms() const;

    // Conjugate by delta^k.
    // Used to speed up calculations compared to the default implementation.
    void delta_conjugate_mut(sint16 k);

    size_t hash() const;

    void of_ballot_sequence(const sint8 *s);
};

typedef FactorTemplate<Underlying> Factor;

typedef BraidTemplate<Factor> Braid;

#ifdef USE_CLN

void ballot_sequence(sint16 n, cln::cl_I k, sint8 *s);

#endif

} // namespace cgarside::band

#endif