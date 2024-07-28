#ifndef OCTAHEDRAL
#define OCTAHEDRAL

#include "cgarside.h"

namespace cgarside::octahedral {

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

    /**
     * @brief Extraction from string.
     *
     * Reads the string `str`, starting at position `pos`, and sets `this` to
     * the corresponding atom.
     *
     * Letting `W = (\s | \t)*` be the language of whitespaces and `Z = -? ([1 -
     * 9] [0 - 9]* | 0)` be the language of integers, accepted strings are those
     * represented by regular expression `\(W Z W,? W Z W\) | Z`, under the
     * additional hypothesis that in the first case the two integers are not
     * equal mod `PresentationParameter'.
     *
     * The first case stands for short generators (double, symmetric,
     * transpositions), and the second one for long generators (transposition of
     * antipodals points).
     *
     * @param str The string to extract from.
     * @param pos The position to start from.
     * @exception `InvalidStringError`: Thrown when there is no subword starting
     * from `pos` that matches `\(W Z W,? W Z W\) | Z`, or if there is one, it
     * matches `\(W Z W,? W Z W\)`, and both integers are equal mod
     * `PresentationParameter'.
     */
    void of_string(const std::string &str, size_t &pos);

    void debug(IndentedOStream &os) const;

    void AssignDCDT(sint16 *x) const;

    void OfDCDT(const sint16 *x);

    // print to os. Be wary, as it side-effects!
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

    std::size_t hash() const {
        std::size_t h = 0;
        for (sint16 i = 1; i <= 2 * get_parameter(); i++) {
            h = h * 31 + permutation_table[i];
        }
        return h;
    }
};

typedef FactorTemplate<Underlying> Factor;

typedef BraidTemplate<Factor> Braid;

} // namespace CGarside

#endif