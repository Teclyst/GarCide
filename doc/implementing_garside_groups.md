# How to implement a Garside group with _GarCide_

## The main idea

The idea behind GarCide is that algorithms for Garside groups are fairly generic, with only base operations (product and lattice operations) depending on the specific group.

Because of that, it is possible to abstract away a lot of implementation details, and write a lot of code at a high level of abstraction, while only requiring that a small amount of methods be implemented.

## The way things work in _GarCide_

in C++, this sort of things can be done using _templates_. Code is written using (abstract) class arguments, and when an instantiation of the object, for given (concrete) class arguments is given, the code is copied, substituting every occurence of the abstract argument by its concrete value.

This is exactly what is done in _GarCide_. Most of the code is written in terms of two template classes, defined in `garcide.h`:

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

It is advised separate the code in two files (if you are working within the _GarCide_ library). On one hand, declarations should go in a _header_ file (ending with `.hpp` or `.h`) in `inc/garcide/groups`. On the other hand, definitions (a.k.a. the real code) should go in an _implementation file_ (ending with `.cpp`) located in `lib/garcide/groups`.

It is also advised to work in a sub-namespace of `cgarcide` (_e.g._ `cgarcide::your_garside_group`).

The existing code follows _Rust_ naming conventions:
*   `snake_case` for variables, functions, class members, namespaces, and file names.
*   `CamelCase` for classes, structs and enums, and everything typing-related.
*   Concise `CamelCase` for template parameters (often single-letter).
*   `SCREAMING_SNAKE_CASE` for constants, compile options and basically everything that occurs at compile time. 

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
    Underlying(Parameter);

    // Returns the lattice's height.
    garcide::sint16 lattice_height() const;

    // Tries to find a substring matching a factor, 
    // starting at the specified position.
    // If it succeeds, sets the factor to what matches that substring,
    // And increases the position by the substring's length.
    // Throws garcide::InvalidStringError if it fails is.
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

If `SomeType` is more complicated however, you may have to specialize the template yourself. This should be done in namespace `cgarcide` specifically (otherwise it will not compile, even in one of its strict sub-namespace). See `dual_complex.h` and `dual_complex.cpp` for an example.

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

Skip this step if you do not want to do that (but you will need to if you want to use for _Braiding 2_).

You simply need to tell _CMake_ that your files should be compiled into the library. Go to `lib/CMakeLists.txt` and add `groups/your_header.hpp` to the `set(HEADERS_LIST ...)` instruction, and `garcide/groups/your_implementation.cpp` to the `add_library(garcide ...)` instruction (of course, change path if you did not add those files there).

### Making your implementation usable for _Braiding 2_

You will need to make a bunch of other modifications to _Braiding 2_ files.

*   You have to tell _CMake_ that (depending on options) your files may have to be compiled into `braiding.exe`. To do this, go to `src/CMakeLists.txt`, and add another `elsif()` branch before the `endif()`:
    ```cmake
    elseif (${USE_FOR_BRAIDING} STREQUAL "YOUR_GROUP")
        message("Add a confirmation message here.")
        target_compile_definitions(braiding.exe PRIVATE -DBRAIDING_CLASS=your_integer)
    ```
*   In `inc/braiding/braiding_option.h`, add another `#elif` branch before the each `#endif`:
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
    This will set the factor and braid classes _Braiding 2_ uses to yours, if it is compiled with flag `-DBRAIDING_CLASS=your_integer`. This will be the case if `USE_FOR_BRAIDING` is set to `"YOUR_GROUP"` with `cmake -DUSE_FOR_BRAING=YOUR_GROUP ..`.
*   In `braiding.cpp`, you should add another `#elif` branch in `explain_braid_input`, where you detail how to enter group elements in the terminal. This is optional, but is good if you are not the only user.
    ```cpp
    #elif BRAIDING_CLASS == your_integer

    ind_cout << ...; // Add an explanation message there.
    ```
*   In `braiding.cpp`, same thing for `explain_braid_parameter_input`, where you detail how to enter group parameters in the terminal.
*   In `braiding.cpp`, add specific UTF8 art for the header, in `print_header`, under an `#elif` branch. (Obviously optional, yet of utmost importance).
*   Edit the __README__. Once again, this is optional but good if you plan to share your version.

### Writing documentation

This project uses _Doxygen_ to generate its documentation. _Doxygen_ is integrated within C++, so you only need to add specific comments over your declarations, then build documentation, and _voilà_.

To look at what these comments should look like, go read _Doxygen_ documentation or look at what is already written. It is also likely that there exist _Doxygen_ extensions for your IDE.

At any rate, here is a reasonable (if somewhat maximal) example of writing documentation with _Doxygen_:

```cpp
/** 
 * @brief A brief description.
 * 
 * A longer description.
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

### Implementing the underlying factors

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

We also need to be able to read a dimension from a string (for instance if we want to use it with _Braiding 2_). It does not really make sense to call this member function from a factor (it does not depend on it). Therefore, we are going to declare that function as `static` to signal that it in fact does not depend on the factor:

`euclidean_lattice.hpp`:
```cpp
class Underlying {

    ...

  public:

    ...
    
    static Parameter parameter_of_string(const std::string&);

}
```

The implementation is pretty straightforward so won't be detailled.