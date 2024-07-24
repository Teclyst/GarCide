#ifndef CENTRALIZER
#define CENTRALIZER

#include "garsideuss.h"

template <class F> class Centralizer {
  private:
    std::vector<F> generators;

  public:
    using Iterator = typename std::vector<F>::iterator;
    using ConstIterator = typename std::vector<F>::const_iterator;

    inline Iterator begin() const {
        return generators.begin();
    }

    inline ConstIterator cbegin() const {
        return generators.begin();
    }

    inline Iterator end() const {
        return generators.end();
    }

    inline ConstIterator cend() const {
        return generators.end();
    }

    inline void insert(const F &f) {
        generators.push_back(f);
    }
};

#endif