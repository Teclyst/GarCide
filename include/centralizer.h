#ifndef CENTRALIZER
#define CENTRALIZER

#include "garsideuss.h"

namespace cgarside::centralizer {

template <class B> class Centralizer {
  private:
    std::vector<B> generators;

  public:
    using Iterator = typename std::vector<B>::iterator;
    using ConstIterator = typename std::vector<B>::const_iterator;

    inline Iterator begin() const { return generators.begin(); }

    inline ConstIterator cbegin() const { return generators.begin(); }

    inline Iterator end() const { return generators.end(); }

    inline ConstIterator cend() const { return generators.end(); }

    inline void insert(const B &f) { generators.push_back(B); }
};

template <class F>
Centralizer<BraidTemplate<F>>
centralizer(const ultra_summit::UltraSummitSet<BraidTemplate<F>> &uss,
            const std::vector<F> &mins, const std::vector<sint16> &prev) {
    BraidTemplate<F> b = uss.at(0, 0);
    typename F::ParameterType p = b.GetParameter();
    Centralizer<BraidTemplate<F>> centralizer;
}

} // namespace cgarside::centralizer

#endif