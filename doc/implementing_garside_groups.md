# How to implement a Garside group with _GarCide_

## The main idea

The idea behind GarCide is that algorithms for Garside groups are fairly generic, with only base operations (product and lattice operations) depending on the specific group.

Because of that, it is possible to abstract away a lot of implementation details, and write a lot of code at a high level of abstraction, while only requiring that a small amount of methods be implemented.

## The way things work in _GarCide_

in C++, this sort of things can be done using _templates_. Code is written using (abstract) class arguments, and when an instantiation of the object, for given (concrete) class arguments is given, the code is copied, substituting every occurence of the abstract argument by its concrete value.

This is exactly what is done in _GarCide_. Most of the code is written in terms of two template classes, defined in `garcide.hpp`:

```cpp
template<class U>
FactorTemplate;

template<class F>
BraidTemplate;
```

`FactorTemplate` expects a class argument `U` representing the underlying objects for canonical factors. `U` must implement a certain number of methods (detailed later) for compilation not to fail, but nothing is required about the actual data structures used. `FactorTemplate<U>` will then provide a full suite of functions and operators that may be expected for a Garside groups canonical factors (see documentation).

`BraidTemplate` expects a class argument `F` representing canonical factors. `F` must implement at least all of `FactorTemplate`'s methods, so in most cases `F` will be a `FactorTemplate<U>`, for some `U` (but that is not mandatory). `BraidTemplate<F>` will then provide methods for manipulating Garside groups elements.

The typical way to implement a Garside group is to write code for a class `Underlying`, and then get factor and braid classes from `FactorTemplate` and `BraidTemplate`.

## Implementing a Garside group

### Code Organization and naming conventions

It is advised separate the code in two files (if you are working within the _GarCide_ library). On one hand, declarations should go in a _header_ file (ending with `.hpp`) in `inc/garcide/groups`. On the other hand, definitions (a.k.a. the real code) should go in an _implementation file_ (ending with `.cpp`) located in `lib/garcide/groups`.

It is also advised to work in a sub-namespace of `cgarcide` (_e.g._ `cgarcide::your_garside_group`).

The existing code follows _Rust_ naming conventions:

* `snake_case` for variables, functions, class members, namespaces, and file names.
* `CamelCase` for classes, structs and enums, and everything typing-related.
* Concise `CamelCase` for template parameters (often single-letter).
* `SCREAMING_SNAKE_CASE` for constants, compile options and basically everything that occurs at compile time.

### Implementing the underlying factors

This is where most of the work occurs. You should implement a class that has following signature (in the sense that it should at least have all those members):

```cpp
class Underlying {

  public:
    // A type representing a group parameter,
    // such the number of strands for standard braids.
    using Parameter = SomeType;

    // Tries to convert a string to a parameter.
    // Throws garcide::InvalidStringError if input is incorrect.
    static Parameter parameter_of_string(const std::string&);

    // Gets the factor's parameter. 
    Parameter get_parameter() const;

    // Constructor.
    // Constructs a new factor whose parameter is the argument.
    // In general, we do not require that the factor be initialized.
    Underlying(Parameter);

    // Returns the lattice's height.
    garcide::sint16 lattice_height() const;

    // Tries to find a substring matching a factor, 
    // starting at the specified position.
    // If it succeeds, sets the factor to what matches that substring,
    // And increases the position by the substring's length.
    // Throws garcide::InvalidStringError if it fails.
    // Should recognize "D" as a Delta.
    void of_string(const std::string&, size_t&);

    // Prints internal data to an outstream (for debugging purposes).
    void debug(garcide::IndentedOStream&) const;

    // Prints the factor to an outstream.
    void print(garcide::IndentedOStream&) const;

    // Sets the factor to the identity.
    void identity();

    // Sets the factor to delta.
    void delta();

    // Compare the factor to another one (equality test).
    bool compare(const Underlying&) const;

    // Computes the product of the factor with another one,
    // under the assumption that it is smaller than delta.
    // (The factor is right-multiplied by the other.)
    Underlying product(const Underlying&) const;

    // Computes the left complement of the factor to another one,
    // under the assumption that the factor right-divides the other one.
    Underlying left_complement(const Underlying&) const;

    // Computes the right complement of the factor to another one,
    // under the assumption that the factor left-divides the other one.
    Underlying right_complement(const Underlying&) const;

    // Conjugates this factor by a power of Delta,
    // and sets it to that conjugate. 
    void delta_conjugate_mut(sint16);

    // Computes the left meet of two factors.
    Underlying left_meet(const Underlying&) const;

    // Computes the right meet of two factors.
    Underlying right_meet(const Underlying&) const;

    // Randomizes the factor.
    void randomize();

    // Hashes the factor.
    size_t hash() const;
};
```

On top of that,

```cpp
cgarcide::IndentedOStream
&garcide::IndentedOStream::operator<< <SomeType>(const SomeType&);
```

must be implemented. This is in particular the case if the same is true, but for `std::ostream` instead of `cgarcide::IndentedOStream`, and therefore in particular for all base types (integers, booleans, ...).

If `SomeType` is more complicated however, you may have to specialize the template yourself. This should be done in namespace `cgarcide` specifically (otherwise it will not compile, even in one of its strict sub-namespace). See `dual_complex.hpp` and `dual_complex.cpp` for an example.

These class members will be called all the time (and it is very likely that most of execution time will be spent executing these), so it is a good place to optimize.

It may happen that you have two possible data structures, with one being more efficient for group operations and the other for lattice operations (typically for dual Garside structures). More complicated functions tend to use more of the former than of the latter, so it is often a good idea to go for the data structure that is best for group operations.

### Getting factor and braid classes

You can then construct a factor and braid classes with

```cpp
using Factor = cgarside::FactorTemplate<Underlying>;

using Braid = cgarside::BraidTemplate<Factor>;
```

They are aliased for brevity.

### Compiling your garside group implementation as a part of the GarCide library

Skip this step if you do not want to do that (but you will need to if you want to use for _Braiding_).

You simply need to tell _CMake_ that your files should be compiled into the library. Go to `lib/CMakeLists.txt` and add `groups/your_header.hpp` to the `set(HEADERS_LIST ...)` instruction, and `garcide/groups/your_implementation.cpp` to the `add_library(garcide ...)` instruction (of course, change path if you did not add those files there).

### Making your implementation usable for _Braiding_

You will need to make a bunch of other modifications to _Braiding_ files.

* You have to tell _CMake_ that (depending on options) your files may have to be compiled into `braiding.exe`. To do this, go to `src/CMakeLists.txt`, and add another `elsif()` branch before the `endif()`:

    ```cmake
    elseif (${USE_FOR_BRAIDING} STREQUAL "YOUR_GROUP")
        message("Add a confirmation message here.")
        target_compile_definitions(braiding.exe PRIVATE -DBRAIDING_CLASS=your_integer)
    ```

* In `inc/braiding/braiding_option.h`, add another `#elif` branch before the each `#endif`:

    ```cpp
    #elif BRAIDING_CLASS == your_integer

    #include "garcide/groups/your_header.hpp"
    ```

    and

    ```cpp
    #elif BRAIDING_CLASS == your_integer

    using Factor = YourFactor;
    using Braid = YourBraid;
    ```

    This will set the factor and braid classes _Braiding_ uses to yours, if it is compiled with flag `-DBRAIDING_CLASS=your_integer`. This will be the case if `USE_FOR_BRAIDING` is set to `"YOUR_GROUP"` with `cmake -DUSE_FOR_BRAING=YOUR_GROUP ..`.
* In `braiding.cpp`, you should add another `#elif` branch in `explain_braid_input`, where you detail how to enter group elements in the terminal. This is optional, but is good if you are not the only user.

    ```cpp
    #elif BRAIDING_CLASS == your_integer

    ind_cout << ...; // Add an explanation message there.
    ```

* In `braiding.cpp`, same thing for `explain_braid_parameter_input`, where you detail how to enter group parameters in the terminal.
* In `braiding.cpp`, add specific UTF8 art for the header, in `print_header`, under an `#elif` branch. (Obviously optional, yet of utmost importance).
* Edit the __README__. Once again, this is optional but good if you plan to share your version.

### Writing documentation

This project uses _Doxygen_ to generate its documentation. _Doxygen_ is integrated within C++, so you only need to add specific comments over your declarations, then build documentation, and _voilà_.

To look at what these comments should look like, go read _Doxygen_ documentation or look at what is already written. It is also likely that there exist _Doxygen_ extensions for your IDE.

At any rate, here is a reasonable (if somewhat maximal) example of writing documentation with _Doxygen_:

```cpp
/** 
 * @brief A brief description.
 * 
 * A longer description (extends brief).
 * 
 * @param x An argument.
 * @tparam T A template argument.
 * @return The return value.
 * @exception Fail Thrown when failing.
 */
template<class T>
ReturnType foo(T x);
```

## An example: Implementing $\mathbb Z^n$ as a Garside group

The following assumes very limited C++ knowledge, and hence will be very detailled.

### $\mathbb Z^n$'s Garside structure

It is known in folklore that $\mathbb Z^n$ has a Garside structure, with its Garside element being $\Delta=(1,\ldots, 1)$, and the sublattice of canonical factors being the cube $[e, \Delta]=\{0,1\}^n$, with the divisibility order being the usual coordinate-wise order.

$[e,\Delta]$ is a section of the canonical surjection $\mathbb Z^n\to(\mathbb Z/2\mathbb Z)^n$, so it is good idea to represent canonical factors as boolean vectors (as it uses less memory). Product is then coordinate-wise XOR, and meets are coordinate-wise MIN (as $\mathbb Z^n$ is abelian, operations are always the same on both sides).

### Owned data

Factor implementation is split between a header file (`inc/garcide/groups/euclidean_lattice.hpp`) and an implementation file (`lib/garcide/groups/euclidean_lattice.cpp`) (go look at source code!).

To avoid naming conflicts with other group implementations, everything is defined within namespace `garcide::euclidean_lattice`. As it is a sub-namespace of `garcide`, `garcide` objects may be called without the initial `garcide::`. As conflicts are avoided that way, we are free to call our underlying factors `Underlying`, without risks of confusion (outside of namespace `garcide::euclidean_lattice`, it would be referred to as `garcide::euclidean_lattice::Underlying`)

As we are in fact implementing Garside structures for $\mathbb Z^n$, for all $n\geqslant 1$, we need a way to get the particular $n$ for a C++ object representing an element of $\mathbb Z^n$. This is done by associating a `size_t` dimension for each factor. In fact, we are going to need it to need to specify a dimension to construct new factors, and in general, implementations are written for families of groups, so a way to find out the particular group must be given. This is done by adding a `Parameter` type to our class.

So we first define in the header a `Parameter` type, and define it as an alias of `size_t`.

`euclidean_lattice.hpp`:

```cpp
class Underlying {

  public:

    using Parameter = size_t;

}
```

A factor is going to need to store some data: its coordinates (seen as an element of $(\mathbb Z/2\mathbb Z)^n$) in the canonical basis, that we are going to store as `std::vector<bool>`, as it is more memory efficient (well, depending on the standard library implementation) than using integers (even `char`s).

We could also use an array of booleans on the stack (a `bool*`), but owned naked pointers are not advised in modern standards, because their ownership semantics aren't clear, in opposition to standard library containers.

Using naked pointers would force us to add a copy constructor and an assignment operator, while by using `std::vector` we get them by default.

So we add a `coordinates` member to our class.

`euclidean_lattice.hpp`:

```cpp
class Underlying {

    ...

  private:

    std::vector<bool> coordinates;

}
```

`coordinates` is `private` to avoid unchecked write access (this guarantees that an initialized factor is always valid), and because it makes the API cleaner (you can do less things with an `Underlying`).

### Implementing parameters

Now that we have added `coordinates`, we can define how to get its dimension from a factor: its `coordinates` size.

So we add a `get_parameter()` member.

`euclidean_lattice.hpp`:

```cpp
class Underlying {

    ...

  public:

    Parameter get_parameter() const;

}
```

`euclidean_lattice.cpp`:

```cpp
Underlying::Parameter Underlying::get_parameter() const {
    return coordinates.size();
}
```

The code is split between a _declaration_ in `euclidean_lattice.hpp`, and a _definition_ in `euclidean_lattice.cpp`.

The `const` keyword signals that `get_parameter` does not modify the `Underlying` it is called on.

We need to be able to read a dimension from a string (for instance if we want to use it for _Braiding_). It does not really make sense to call this member function from a factor (it does not depend on it). Therefore, we are going to declare that function as `static` to signal that it in fact does not depend on the factor:

`euclidean_lattice.hpp`:

```cpp
class Underlying {

    ...

  public:

    ...
    
    static Parameter parameter_of_string(const std::string&);

}
```

The implementation is rather straightforward so won't be detailled.

On the other side, we also need to be able to print parameters. This is done by overloading operator `<<` for `IndentedOStream` and `Parameter`; but in this case, `Parameter` desugars to base type `size_t`. Therefore, there is a specialization of `<<` for `std::ostream` and `Parameter`, so that is also the case for `IndentedOStream` (because there is a generic specialization), and we don't have to do anything else.

### Adding a constructor

We are going to need a way to construct new factors. This is done by adding a _constructor_.

`euclidean_lattice.hpp`:

```cpp
class Underlying {

    ...

  public:

    Underlying(Parameter n);

}
```

`euclidean_lattice.cpp`:

```cpp
Underlying::Underlying(Underlying::Parameter n) : coordinates(n, false) {}
```

New factors are constructed by creating a new factor and filling it with the right amount of `false`s (we need `coordinates` to be the right length so that `get_parameter` returns the correct dimension).

### Handling input

We need a function that converts strings to factors. This is done with the following function:

`euclidean_lattice.hpp`:

```cpp
class Underlying {

    ...

  public:

    ...
    
    void of_string(const std::string &str, size_t &pos);

}
```

`of_string`'s behaviour is a bit more complicated that `parameter_of_string`'s, because it incorporates a position `pos`.

The reason for this is that the corresponding `BraidTemplate` function tries to extract factors from a string one by one, by calling `of_string`. Therefore `of_string` needs to keep track of the position within the string (in order to know where to start from), and it must update it after it succesfully extracts a factor. Let's take a look at the code.

Also notice that while `parameter_of_string` returns a `Parameter`, `of_string` does not return an `Underlying`, and instead modifies the factor it is called on.

`euclidean_lattice.cpp`:

```cpp
void Underlying::of_string(const std::string &str, size_t &pos) {
    Parameter n = get_parameter();

    std::smatch match;

    if (std::regex_search(
            str.begin() + pos, str.end(), match,
            std::regex{"(?:e[\\s\\t]*_?[\\s\\t]*)?(" + number_regex + ")"},
            std::regex_constants::match_continuous)) {
        sint16 i;
        try {
            i = std::stoi(match[1]);
        } catch (std::out_of_range const &) {
            throw InvalidStringError("Index is too big!\n" + match.str(1) +
                                     " cannot be converted to a C++ integer.");
        }
        pos += match[0].length();
        if ((i >= 0) && (i < int(n))) {
            identity();
            coordinates[i] = true;
        } else {
            throw InvalidStringError(
                "Invalid index for canonical base vector!\n" + match.str(1) +
                " is not in [0, " + std::to_string(n) + "[.");
        }
    } else if (std::regex_search(str.begin() + pos, str.end(), match,
                                 std::regex{"D"},
                                 std::regex_constants::match_continuous)) {
        pos += match[0].length();
        delta();
    } else {
        throw InvalidStringError(std::string(
            "Could not extract a factor from\n\"" + str.substr(pos) +
            "\"!\nA factor should match regex ('e' '_'?)? Z | 'D',\nwhere Z "
            "matches integers and ignoring whitespaces."));
    }
}
```

Factors (and _a fortiori_ braids) are written as products in the atoms, therefore it is enough to only recognize atoms in `of_string` (as even though a non-atomic factor will not be read as one by `of_string`, it will anyway be parsed by the corresponding `BraidTemplate` member function as single-factor braid, which is pretty much the same).

We use regular expressions (regexes from now on) to specify what patterns should be recognized as factors. It is advised to use them because they are a very convenient way to specify reasonably simple patterns, and are typically enough to recognize atoms.

Typically, in this case we want to allow several different ways of entering the same atom $e_i$: it's going to be printed as `ei`, but it could be equally be reasonnable to enter it as shorter `i`, or longer `e i`, `e_i`, `e _ i`, or some wilder alternative with $42$ tabs. (For convenience sake however we should make sure that output is compatible with input, so `ei` should be a valid input).

Now the great thing about regular expressions is that they make it very easy to specify all these kinds of patterns at the same time: it is simply expressed by regular expression

```math
(\texttt e (\backslash\texttt{s} \mid \backslash\texttt{t})^*\texttt \_?(\backslash\texttt{s} \mid \backslash\texttt{t})^*)?([\texttt 1 - \texttt 9][\texttt 0 - \texttt 9]^* \mid \texttt 0).
```

The second part, $[\texttt 1 - \texttt 9][\texttt 0 - \texttt 9]^* \mid \texttt 0$, matches either $\texttt 0$ or a sequence of digits starting by a non-zero one, _i.e._ a positive integer. (Here $[a - b]$ matches the _range_ of characters between $a$ and $b$).

The first is optional (the whole group is tagged with a postfix quotation mark), and is composed of an $\texttt e$, followed by a sequence of whitespaces ($\backslash\texttt{s}$) and tabs ($\backslash\texttt{t}$), with an optional underscore in the middle.

Notice than in all cases, the meaningful part does not change: it's always the second one. This can be expressed in C++ regexes using _capture groups_. `std::regex_search`, when executed with positions, and flag `std::regex_constants::match_continuous`, will try to find a match for the regex starting at the first provided position and ending before the second. When it finds one, it will also extract substrings from all subexpressions enclosed by parenthesis in the big one, and store them in `match`. These parenthesis denote capture groups.

In our case, we only want to capture the integer (the rest is not important). But we also want to group the first part (to apply $?$ on all of it). In C++ this can be done with a non-capture group, that is formed by enclosing a group between `(?:` and `)`. So we can write that regex, in C++, as `(?:e(?:\\s|\\t)*_?(?:\\s|\\t)*)?([1-9][0-9]*|0)`, or as the slightly shorter, and equivalent, `(?:e[\\s\\t]*_?[\\s\\t])?([1-9][0-9]*|0)` (as `[]` denotes a set of characters). We already defined `number_regex` as the string `\\-?[1-9][0-9]*|0` elsewhere, so we shorten it even more to the version that is in the code (althought this means we will have to test for positivity later on).

Notice that we must use double `\` to escape them.

All that's left to do is now to try to call `std::regex_search` to find a factor, extract the submatch that corresponds to an integer, and then convert it to that integer, construct the atom, and increase `pos` by the match's length. But the conversion may fail: the provided integer may be strictly negative, or too big - bigger than the dimension, or even worse, so big that it cannot be represented as a C++ integer. In the latter case in particular, an `std::out_of_range` exception will be thrown, and if never caught this would result in the program (_e.g._, _Braiding_) crashing, which is not good practice.

It is much better to instead catch that exception, and then throw an `InvalidStringError`. That way, the error will be caught at the _Braiding_ strata, and the program will fail in a controlled, non-crashing manner, simply displaying the error message that was sent up before asking the user to try again.

Then we also have to recognize `"D"` as $\Delta$. This is done by adding a second case.

It may happen that the input string is incorrect and does not match either of those cases. In this case, the program should once again fail in a controlled manner, using `InvalidStringError` to send an error message back up.

### Handling output

Output is handled by two different functions:

`euclidean_lattice.hpp`

```cpp
class Underlying {

    ...

  public:

    ...
    
    void print(IndentedOStream &os) const;

    void debug(IndentedOStream &os) const;
}
```

They serve different purposes: `print` is meant to be used as the ouput functions (the one that will print to a terminal), and will write factors in a conventional form (typically as product of atoms), that should be parsed by `of_string` to the original factor. In our case, vector $(1, 0, 1)$ would for instance be printed as `e0 e2`.

On the other hand, `debug` is meant to be used (surprinsingly) as a debug tool, and should print the factor's internal data. So our example $(1, 0, 1)$ would be debugged as:

```txt
{   coordinates:
        [true, false, true]
}
```

### Implementing lattice operations

To finish implementing the lattice structure, we must also add functions to set a factor to both the identity and delta, and two functions for lateral meets.

`euclidean_lattice.hpp`:

```cpp
class Underlying {

    ...

  public:

    ...
    
    void identity();

    void delta();

    Underlying left_meet(const Underlying &b) const;

    Underlying right_meet(const Underlying &b) const; 
}
```

The first two are fairly trivial: just set every coordinate to `false` or `true` respectively.

Meets aren't very difficult either: they are both coordinate-wise AND.

`euclidean_lattice.cpp`:

```cpp
Underlying Underlying::left_meet(const Underlying &b) const {
    Underlying meet(get_parameter());
    for (size_t i = 0; i < get_parameter(); i++) {
        meet.coordinates[i] = coordinates[i] && b.coordinates[i];
    }
    return meet;
}
```

Both meets are in fact the same, as $\mathbb Z^n$ is abelian, so `right_meet` can be defined as a function returning the result of `left_meet`, called on the same arguments; but that would be slightly inefficient, because calling `right_meet` results in an useless call, before the relevant function comes into play. Calls are not free (though not very expensive), because they require writing down an activation table to the stack. Therefore, it would be great if we could avoid that useless call, without having to copy and paste code.

This can be done in C++ by _inlining_ the function. Doing this will replace at compile time every occurrence of a function by its code. That way, the second call can be ignored!

A C++ function can be inlined by adding the `inline` keyword before its declaration, and providing its definition directly in the header file. In the case at hand, it gives the following code:

`euclidean_lattice.hpp`:

```cpp
class Underlying {

    ...

  public:

    ...
    
    inline Underlying right_meet(const Underlying &b) const {
        return left_meet(b);
    }; 
}
```

In general, it is better to inline all small functions (and in fact, basically all functions in this implementation would have been legitimate inlining candidates). It is, however, best to avoid inlining bigger functions. The drawbacks to inlining are that it makes compilation longer (as it is the time at which inlining is performed), may result in longer binary files, and can pollute APIs with implementation details. When inlining big functions, performance may also actually worsen because code locality decreases.

### Comparing factors

We'll also have to implement equality tests, because it is not automatically derived, unlike copy and move semantics.

However, it is quite easy to do (as we just have to compare coordinates, and comparison for `std::vector` already does what we want to), so we just `inline` it once again.

`euclidean_lattice.hpp`:

```cpp
class Underlying {

    ...

  public:

    ...

    inline bool compare(const Underlying &b) const {
        return coordinates == b.coordinates
    };

}
```

### The group structure

Next we have to add the group operations: product and lateral complement operations; these three are actually the same, because $(\mathbb Z/2\mathbb Z)^n$ is abelian and all of its elements are of order $2$. That operation is coordinate wise XOR, so we just implement `product` that way, then define the other two as `product`, inlining them.

There is also conjugation by delta, which does nothing as the group is abelian, and is therefore safe to inline.

`euclidean_lattice.hpp`

```cpp

class Underlying {

    ...

  public:

    ...

    Underlying product(const Underlying &b) const;

    inline Underlying left_complement(const Underlying &b) const {
        return product(b);
    }

    inline Underlying right_complement(const Underlying &b) const {
        return product(b);
    }

    inline void delta_conjugate_mut(__attribute__((unused)) sint16 k) {}
}
```

Notice the `__attribute__((unused))` flag before `k` that tells the compiler that `k` is unused, to avoid warnings.

### The end: randomizing and hashing

`euclidean_lattice.hpp`

```cpp

class Underlying {

    ...

  public:

    ...

    void randomize();

    size_t hash() const;
}
```

Randomizing is pretty simple: just do it coordinate-wise.

For hashing, we want a function that is easy to compute, and collision-resistant; we need hashing for several data structures whose efficiency really depends on `hash` being resistant to collisions. In our case, we actually have an injection of canonical factors to $[\![0, 2^n[\![$, given by interpreting a factor as a polynomial in $2$.

`euclidean_lattice.cpp`:

```cpp
size_t Underlying::hash() const {
    size_t hash = 0;
    for (size_t i = 0; i < get_parameter(); i++) {
        hash <<= 1;
        hash += (at(i) ? 1 : 0);
    }
    return hash;
}
```

In general, an often good way to produce hash functions is to interpret in some way or another an object as a polynomial in some reasonably big prime number (often $31$).
