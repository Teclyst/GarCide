#include "dihedral_braid.h"
#include <string>
#include <iostream>

struct NotBelow
{
};

namespace CGarside
{

  sint16 IDualBraidUnderlying::GetParameter() const
  {
    return PresentationParameter;
  }

  sint16 IDualBraidUnderlying::LatticeHeight() const
  {
    return 2;
  }

  IDualBraidUnderlying &IDualBraidUnderlying::Assign(const IDualBraidUnderlying &a)
  {
    Type = a.Type;
    Axis = a.Axis;
    return *this;
  };

  IDualBraidUnderlying &IDualBraidUnderlying::operator=(const IDualBraidUnderlying &a)
  {
    return Assign(a);
  }

  IDualBraidUnderlying::IDualBraidUnderlying(sint16 n)
      : PresentationParameter(n),
        Type(0),
        Axis(0)
  {
  }

  IDualBraidUnderlying::IDualBraidUnderlying(const IDualBraidUnderlying &a)
      : PresentationParameter(a.GetParameter()),
        Type(a.Type),
        Axis(a.Axis)
  {
  }

  void IDualBraidUnderlying::Print(std::ostream &os) const
  {
    if (Type == 1)
    {
      os << "D";
    }
    else if (Type == 2)
    {
      os << Axis;
    }
  }

  void IDualBraidUnderlying::OfString(std::string &str)
  {
    if (str == std::string("D"))
    {
      Type = 1;
    }
    else if (str.find_first_not_of("0123456789") == std::string::npos)
    {
      Type = 2;
      Axis = std::stoi(str) % GetParameter();
    }
    else
    {
      throw InvalidStringError();
    }
  }

  IDualBraidUnderlying IDualBraidUnderlying::LeftMeet(const IDualBraidUnderlying &b) const
  {
    if ((Type == 0) || (b.Type == 0) || ((Type == 2) && (b.Type == 2) && (Axis != b.Axis)))
    {
      return IDualBraidUnderlying(GetParameter());
    }
    if (Type == 1)
    {
      IDualBraidUnderlying c = b;
      return c;
    }
    IDualBraidUnderlying c = *this;
    return c;
  }

  IDualBraidUnderlying IDualBraidUnderlying::RightMeet(const IDualBraidUnderlying &b) const
  {
    return LeftMeet(b);
  }

  void IDualBraidUnderlying::Identity()
  {
    Type = 0;
  }

  void IDualBraidUnderlying::Delta()
  {
    Type = 1;
  }

  bool IDualBraidUnderlying::Compare(const IDualBraidUnderlying &b) const
  {
    return ((Type == b.Type) && ((Type != 2) || (Axis == b.Axis)));
  }

  IDualBraidUnderlying IDualBraidUnderlying::Product(const IDualBraidUnderlying &b) const
  {
    IDualBraidUnderlying f = IDualBraidUnderlying(GetParameter());
    if (GetParameter() != b.GetParameter())
    {
      throw NonMatchingIndexes(GetParameter(), b.GetParameter());
    }
    if (Type == 0)
    {
      f = b;
    }
    else if (b.Type == 0)
    {
      f = *this;
    }
    else if ((Type == 2) && (b.Type == 2) && ((Axis - b.Axis + GetParameter()) % GetParameter() == 1))
    {
      f.Delta();
    }
    else
    {
      throw NotBelow();
    }
    return f;
  };

  IDualBraidUnderlying IDualBraidUnderlying::LeftComplement(const IDualBraidUnderlying &b) const
  {
    IDualBraidUnderlying f = IDualBraidUnderlying(GetParameter());
    if (b.Type == 1)
    {
      if (Type == 0)
      {
        f.Delta();
      }
      else if (Type == 2)
      {
        f.Type = 2;
        if (Axis == GetParameter() - 1)
        {
          f.Axis = 0;
        }
        else
        {
          f.Axis = Axis + 1;
        }
      }
    }
    else if (b.Type == 2)
    {
      if (Type == 0)
      {
        f = b;
      }
      else if ((Type == 2) && (Axis == b.Axis))
      {
        f.Identity();
      }
      else
      {
        throw NotBelow();
      }
    }
    else
    {
      if (Type != 0)
      {
        throw NotBelow();
      }
    }
    return f;
  };

  IDualBraidUnderlying IDualBraidUnderlying::RightComplement(const IDualBraidUnderlying &b) const
  {
    IDualBraidUnderlying f = IDualBraidUnderlying(GetParameter());
    if (b.Type == 1)
    {
      if (Type == 0)
      {
        f.Delta();
      }
      else if (Type == 2)
      {
        f.Type = 2;
        if (Axis == 0)
        {
          f.Axis = GetParameter() - 1;
        }
        else
        {
          f.Axis = Axis - 1;
        }
      }
    }
    else if (b.Type == 2)
    {
      if (Type == 0)
      {
        f = b;
      }
      else if ((Type == 2) && (Axis == b.Axis))
      {
        f.Identity();
      }
      else
      {
        throw NotBelow();
      }
    }
    else
    {
      if (Type != 0)
      {
        throw NotBelow();
      }
    }
    return f;
  };

  IDualBraidUnderlying IDualBraidUnderlying::DeltaConjugate(sint16 k) const
  {
    IDualBraidUnderlying under = IDualBraidUnderlying(*this);
    sint16 n = GetParameter();
    if (Type != 2)
    {
      under = *this;
    }
    else
    {
      if (k > 0)
      {
        k = k - k * n;
      }
      under.Type = 2;
      under.Axis = (Axis - 2 * k) % n;
    }

    return under;
  }

  void IDualBraidUnderlying::Randomize()
  {
    sint16 rand = std::rand() % (GetParameter() + 1);
    if (rand == GetParameter())
    {
      Type = 0;
    }
    else if (rand == GetParameter() + 1)
    {
      Type = 1;
    }
    else
    {
      Type = 2;
      Axis = rand;
    }
  }

  std::vector<IDualBraidUnderlying> IDualBraidUnderlying::Atoms() const
  {
    sint16 n = GetParameter();
    IDualBraidUnderlying atom = IDualBraidUnderlying(n);
    atom.Type = 2;
    std::vector<IDualBraidUnderlying> atoms;
    for (sint16 i = 0; i < n; i++)
    {
      atom.Axis = i;
      atoms.push_back(atom);
    }
    return atoms;
  }

}