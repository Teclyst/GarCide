As of 2024, code was factored so as to allow easy implementations of other Garside structures, taking advantage of the fact that algorithms for Garside Groups are generic, once a very limited subset of operations are define. 

## How the code is currently structured

Most Functions in the Braiding part are templated, and can be instantiated for classes that implement a certain set of base method needed to make computations in an abstraction of a braid group.

In turn, the _CGarside_ part of the project implement templated classes for (generic) canonical factors an braids that make implementing Garside structures rather simple: basically, once base operations are implemented for a specific Garside structure, most of the rest can be derived from it at a general, abstract level.

Note that by design, we do not assume anything about the way canonical factors are actually represented. This is because base operations fundamentally depend on the group being considered, and finding the kind of representation that will be computationally efficient is something we believe should be left to the user if it's not for a group we provide code for.

## So how do I implement my favourite Garside group?

It's very simple actually!

You'll only need to define a class for the underlying representation of canonical factors, with a total of 12 methods (but some of them should be extremely easy to implement).

Here is the signature needed:

```c++
class Underlying {

  // A copy constructor.
  // (Should have the same name as the class.)
  Underlying(const Underlying &u);

  // A method that initializes based on a string str.
  void OfString(const std::string &str);

  // A method that prints to an output stream os.
  void Print(const std::ostream &os);

  // A method that initializes to the neutral element.
  void Neutral();

  // A method that initializes to Delta.
  void Delta();

  // A method to check equality with another factor.
  bool Compare(const Underlying &v) const;

  // A method that computes the left meet with another factor v.
  Underlying LeftMeet(const Underlying &v) const;

  // A method that computes the right meet with another factor v.
  Underlying RightMeet(const Underlying &v) const;

  // A method that computes the product with another factor v.,
  // assuming that the product is still below Delta.
  Underlying Product(const Underlying &v) const;

  // A method such that, given two factors u and v = wu,
  // u.LeftComplement(v) returns w.
  Underlying LeftComplement(const Underlying &v) const;

  // A method such that, given two factors u and v = uw,
  // u.RightComplement(v) returns w.
  Underlying RightComplement(const Underlying &v) const;

  // A method that returns the list of the generators of the group.
  std::list<Underlying> Atoms() const;

}
```


Realistically you'll also need other constructors, as with only those methods you won't be able to create new factors (only to copy them). These constructors may depend on representation-dependent arguments (_e.g._ the number of strands for braids); that is why they cannot be used in a general framework.

You'll maybe have noticed that `Atoms` does not arguments; this might seem surprising because you'd expect that the group generators may depend on some parameter (say the number of strands). But adding parameters would result in the same kind of problems: the generic part cannot be assumed to know how the Garside structure will be parameterized, so we cannot have non-generic arguments for `Atoms`.

Therefore, to be able to retrieve these informations, they will have to be encapsulated within the underlying representation of factors. Note that this solves our problem: in the case of braids the number of strands can be extracted from the factor the `Atoms` method is called upon, and this can be used to generate the corresponding generators.

Once you have a class (called, say, `Underlying`) that satisfies that signature, you can then define another class that will derive more methods for factors: this is done by using the template class `Factor`; this new class can then be used to derive another, _Braiding_-compatible, class representing the elements of your Garside Group.

You may want to specialize the Factor method `DeltaConjugate` if you have a more efficient way to compute conjugates by powers of Delta than simply iterating conjugations.

```c++
typedef Factor<Underlying> CanonicalFactor;

typedef Braid<CanonicalFactor> Element;
```