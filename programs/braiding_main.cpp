/*
    Copyright (AC) 2004 Juan Gonzalez-Meneses.

    This file is part of Braiding.

    Braiding is free software; you can redistribute Ait and/or modify
    Ait under the terms of the GNU General Public License as published by
    the Free Software Foundation; either version 2 of the License, or
    any later version.

    Braiding is distributed in the hope that Ait will be useful,
    but WITHOUT ANY WARRANTY; without even the implied warranty of
    MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
    GNU General Public License for more details.

    You should have received a copy of the GNU General Public License
    along with Braiding; if not, write to the Free Software
    Foundation, Inc., 59 Temple Place, Suite 330, Boston, MA 02111-1307 USA
*/
/*
    braiding_main.cpp,  v 1.0.   04/10/2004
    Juan Gonzalez-Meneses <meneses(at)us.es>
*/


#include "cbraid.h"
#include "braiding.h"
#include <iostream>
#include <iomanip>
#include <fstream>
#include <ctime>
#include <cstdlib>
#include <stdio.h>

enum Presentation {
	Artin,
	Band
}

int main()
{
  using namespace CBraid;
  using CBraid::sint16;  // Avoid ambiguity with CLN.
  using namespace Braiding;
  using namespace std;

	Presentation pres = Artin;
 	char  c, *file;
  sint16 p=0, power=1, power2=1, i, j, n, repeat=0,
    size, type, rigidity, iteration;
  list<sint16> word, word2, graph, graphinv;
  list<sint16>::iterator itw, itg;

  ArtinBraid AB=ArtinBraid(1), AB1=ArtinBraid(1),
    AB2=ArtinBraid(1), AB3=ArtinBraid(1), AC=ArtinBraid(1);
  list<ArtinBraid> Asss, Atraj, ACent, Avertices, Abarrows;
  list<ArtinBraid>::iterator Ait, Aitb, Aitb;

  list<list<ArtinBraid> > Auss, Asc;
  list<list<ArtinBraid> >::iterator Aituss, Aituss2;

  list<ArtinFactor> AMin, Aarrows;
  list<ArtinFactor>::iterator Aitf, Aitf2;

  bool conj;
  ofstream f;

  while(1)
    {
      cout << endl << endl << endl << endl
	   << "--------------------------------------------------------" << endl
	   << "----------------  This is Braiding 1.0  ----------------" << endl
	   << "--------------------------------------------------------" << endl
	   << "-----|  Copyright (AC) 2004 Juan Gonzalez-Meneses  |-----" << endl
	   << "-----| Braiding comes with ABSOLUTELY NO WARRANTY |-----" << endl
	   << "-----|           This is free software            |-----" << endl
	   << "-----| See GNU General Public License in GPL.txt  |-----" << endl
	   << "--------------------------------------------------------" << endl
	   << endl
	   << "l: Left Normal Form          r: Right Normal Form       " << endl
	   << endl
	   << "p: Permutation               x: Crossing numbers        " << endl
	   << endl
	   << "v: Least Common Multiple     ^: Greatest Common Divisor " << endl
	   << endl
	   << "s: Super Summit Set          z: Centralizer             " << endl
	   << endl
	   << "e: Conjugacy Test            u: Ultra Summit Set        " << endl
	   << endl
     << "t: Set of Sliding Circuits   a: Ask for Powers (On/Off) " << endl
     << endl
     << "d: Dual Mode (On/Off)        q: Quit                    " << endl;

      while(1)
	{
	  power=1;
	  power2=1;
	  cout << endl
	       << "--------------------------------------------------------"
	       << endl
	       << endl << "Choose an option: (type '?' for help) ";
	  cin >> ws >> c;

	  /////////////////////////////////////////////////////////////

	  if(c=='l' || c=='r')
	    {
	      n=ReadIndex();
	      word=ReadWord(n);
	      AB=ArtinBraid(n);
	      AB=WordToBraid<CBraid::ArtinPresentation>(word,n);

	      if(p)
		{
		  power=ReadPower();
		  AB=RaisePower(AB,power);
		}

	      file=ReadFileName();

	      if(c=='l')
		{
		  cout << endl << "The Left Normal Form is: " << endl << endl;
		  PrintBraidWord(AB.MakeLCF());
		  cout << endl;

		  f.open(file);
		  f << endl << "The Left Normal Form of the braid on "
		    << n << " strands" << endl << endl;
		  f.close();
		  PrintWord(word,n,power,file);
		  f.open(file,ios::app);

		  f << endl << endl;

		  f << "is: " << endl << endl;
		  f.close();
		  PrintBraidWord(AB.MakeLCF(),file);
		}

	      if(c=='r')
		{
		  cout << endl << "The Right Normal Form is: "
		       << endl << endl;
		  PrintBraidWord(AB.MakeRCF());
		  cout << endl;

		  f.open(file);
		  f << endl << "The Right Normal Form of the braid on "
		    << n << " strands" << endl << endl;
		  f.close();
		  PrintWord(word,n,power,file);
		  f.open(file,ios::app);

		  f << endl << endl;

		  f << "is: " << endl << endl;
		  f.close();
		  PrintBraidWord(AB.MakeRCF(),file);
		}

	      word.clear();
	    }

	  ////////////////////////////////////////////////

	  if(c=='p')
	    {
	      n=ReadIndex();
	      word=ReadWord(n);
	      if(p)
		power=ReadPower();

	      AB=ArtinBraid(n);
	      AB=WordToBraid<ArtinPresentation>(word,n);

	      if(p && power!=1)
		AB=RaisePower(AB,power);

	      AB.MakeLCF();

	      file=ReadFileName();
	      f.open(file);

	      ArtinFactor F=ArtinFactor(n, AB.LeftDelta);

	      list<ArtinFactor>::iterator Aitf = AB.FactorList.begin();
	      while (Aitf != AB.FactorList.end())
		F *= *(Aitf++);

	      sint16 *table=new sint16[n+1];
	      for(i=1; i<=n; i++)
		table[i]=0;

	      cout << endl << "The permutation associated to this braid is:"
		   << endl << endl;

	      f << "The permutation associated to the braid on "
		<< n << " strands" << endl << endl;
	      f.close();
	      PrintWord(word,n,power,file);
	      f.open(file,ios::app);
	      f << endl << endl << " is:" << endl << endl;

	      if(F.CompareWithIdentity())
		{
		  cout << "Trivial." << endl << endl;
		  f <<  "Trivial.";
		}
	      else
		{
		  for(i=1; i<=n; i++)
		    {
		      if(F[i]!=i && table[i]==0)
			{
			  cout << "(" << i;
			  f << "(" << i;
			  j=i;
			  while(F[j]!=i)
			    {
			      j=F[j];
			      table[j]=1;
			      cout << "," << j;
			      f << "," << j;
			    }
			  cout << ")";
			  f << ")";
			}
		    }
		}
	      cout << endl << endl;
	      f.close();
	      word.clear();
	    }

	  ///////////////////////////////////////////////////

	  if(c=='x')
	    {
	      n=ReadIndex();
	      word=ReadWord(n);
	      if(p)
		power=ReadPower();

	      file=ReadFileName();

	      sint16 **cross= new sint16 *[n];
	      for(i=1; i<n; i++)
		cross[i]=new sint16[n+1];

	      Crossing(word,n,power,cross);

	      f.open(file);

	      cout << "The crossing numbers of this braid are:"
		   << endl << endl << "    ";
	      f << "The crossing numbers of the braid on " << n
		<< " strands" << endl << endl;
	      f.close();
	      PrintWord(word,n,power,file);
	      f.open(file,ios::app);
	      f << endl << endl << "are: " << endl << endl << "    ";

	      for(i=2; i<=n; i++)
		{
		  cout << setw(3) << i;
		  f << setw(3) << i;
		}
	      cout.fill('-');
	      cout << endl << "   +" << setw(3*(n-1)) << "-" <<  endl;
	      cout.fill(' ');
	      f.fill('-');
	      f << endl << "   +" << setw(3*(n-1)) << "-" <<  endl;
	      f.fill(' ');

	      for(i=1; i<n; i++)
		{

		  cout << setw(3) << i << "|" << setw(3*i) << cross[i][i+1];
		  f << setw(3) << i << "|" << setw(3*i) << cross[i][i+1];

		  for(j=i+2; j<=n; j++)
		    {

		      cout << setw(3) << cross[i][j];
		      f << setw(3) << cross[i][j];
		    }
		  if(i<n-1)
		    {
		      cout << endl << "   |" << endl;
		      f << endl << "   |" << endl;
		    }
		}
	      cout << endl << endl;
	      f.close();

	      for(i=1; i<n; i++)
		delete[] cross[i];
	      delete[] cross;
	      word.clear();
	    }

	  //////////////////////////////////////////////////////////

	  if(c=='^' || c=='v')
	    {
	      n=ReadIndex();
	      word=ReadWord(n);
	      if(p)
		power=ReadPower();

	      AB1=ArtinBraid(n);
	      AB1=WordToBraid<ArtinPresentation>(word,n);
	      if(p && power!=1)
		AB1=RaisePower(AB1,power);

	      AB1.MakeLCF();

	      word2=ReadWord(n);
	      if(p)
		power2=ReadPower();

	      AB2=ArtinBraid(n);
	      AB2=WordToBraid<ArtinPresentation>(word2,n);
	      if(p && power2!=1)
		AB2=RaisePower(AB2,power2);

	      AB2.MakeLCF();

	      file=ReadFileName();

	      AB=ArtinBraid(n);

	      if(c=='v')
		AB=LeftWedge(AB1,AB2);
	      else
		AB=LeftMeet(AB1,AB2);

	      if(c=='v')
		cout << "The lcm of these braids is:" << endl << endl;
	      else
		cout << "The gcd of these braids is:" << endl << endl;
	      PrintBraidWord(AB);
	      cout << endl << endl;

	      f.open(file);

	      if(c=='v')
		f << "The lcm of the braids on " << n << " strands"
		  << endl << endl;
	      else
		f << "The gcd of the braids on " << n << " strands"
		  << endl << endl;
	      f.close();
	      PrintWord(word,n,power,file);
	      f.open(file,ios::app);

	      f << endl << endl << "and" << endl << endl;

	      f.close();
	      PrintWord(word2,n,power2,file);
	      f.open(file,ios::app);

	      f << endl << endl << "is:" << endl << endl;

	      PrintBraidWord(AB,file);

	      f.close();
	      word.clear();
	    }

	  //////////////////////////////////////////////////////////

	  if(c=='e')
	    {
	      n=ReadIndex();
	      word=ReadWord(n);
	      if(p)
		power=ReadPower();

	      AB1=ArtinBraid(n);
	      AB1=WordToBraid<ArtinPresentation>(word,n);
	      if(p && power!=1)
		AB1=RaisePower(AB1,power);

	      AB1.MakeLCF();

	      word2=ReadWord(n);
	      if(p)
		power2=ReadPower();

	      AB2=ArtinBraid(n);
	      AB2=WordToBraid<ArtinPresentation>(word2,n);
	      if(p && power2!=1)
		AB2=RaisePower(AB2,power2);

	      AB2.MakeLCF();

	      file=ReadFileName();

	      AC=ArtinBraid(n);

	      conj=AreConjugate(AB1,AB2,AC);

	      if(conj)
		{
		  cout << endl << "These braids are conjugate." << endl << endl
		       << "A conjugating braid is: ";
		  PrintBraidWord(AC);
		  cout << endl;
		}
	      else
		{
		  cout << endl << "These braids are not conjugate." << endl;
		}

	      f.open(file);

	      f << "The braids on " << n << " strands" << endl << endl;
	      f.close();
	      PrintWord(word,n,power,file);
	      f.open(file,ios::app);

	      f << endl << endl << "and" << endl << endl;

	      f.close();
	      PrintWord(word2,n,power2,file);
	      f.open(file,ios::app);



	      f << endl << endl;

	      if(conj)
		{
		  f << "are conjugate." << endl << endl
		    << "A conjugating braid is" << endl << endl;
		  f.close();
		  PrintBraidWord(AC,file);
		}
	      else
		{
		  f << "are not conjugate.";
		  f.close();
		}

	      word.clear();
	    }

	  /////////////////////////////////////////////////////////////

if(c=='d')
	    {
	      n=ReadIndex();
	      word=ReadWord(n);
	      if(p)
		power=ReadPower();

	      AB1=ArtinBraid(n);
	      AB1=WordToBraid<ArtinPresentation>(word,n);
	      if(p && power!=1)
		AB1=RaisePower(AB1,power);

	      AB1.MakeLCF();

	      word2=ReadWord(n);
	      if(p)
		power2=ReadPower();

	      AB2=ArtinBraid(n);
	      AB2=WordToBraid<ArtinPresentation>(word2,n);
	      if(p && power2!=1)
		AB2=RaisePower(AB2,power2);

	      AB2.MakeLCF();

	      file=ReadFileName();

	      AC=ArtinBraid(n);

	      conj=AreConjugateSC(AB1,AB2,AC);

	      if(conj)
		{
		  cout << endl << "These braids are conjugate." << endl << endl
		       << "A conjugating braid is: ";
		  PrintBraidWord(AC);
		  cout << endl;
		}
	      else
		{
		  cout << endl << "These braids are not conjugate." << endl;
		}

	      f.open(file);

	      f << "The braids on " << n << " strands" << endl << endl;
	      f.close();
	      PrintWord(word,n,power,file);
	      f.open(file,ios::app);

	      f << endl << endl << "and" << endl << endl;

	      f.close();
	      PrintWord(word2,n,power2,file);
	      f.open(file,ios::app);



	      f << endl << endl;

	      if(conj)
		{
		  f << "are conjugate." << endl << endl
		    << "A conjugating braid is" << endl << endl;
		  f.close();
		  PrintBraidWord(AC,file);
		}
	      else
		{
		  f << "are not conjugate.";
		  f.close();
		}

	      word.clear();
	    }

	  /////////////////////////////////////////////////////////////
if(c=='c')
	    {
	      n=ReadIndex();
	      word=ReadWord(n);
	      if(p)
		power=ReadPower();

	      AB1=ArtinBraid(n);
	      AB1=WordToBraid<ArtinPresentation>(word,n);
	      if(p && power!=1)
		AB1=RaisePower(AB1,power);

	      AB1.MakeLCF();

	      word2=ReadWord(n);
	      if(p)
		power2=ReadPower();

	      AB2=ArtinBraid(n);
	      AB2=WordToBraid<ArtinPresentation>(word2,n);
	      if(p && power2!=1)
		AB2=RaisePower(AB2,power2);

	      AB2.MakeLCF();

	      file=ReadFileName();

	      AC=ArtinBraid(n);

	      conj=AreConjugateSC2(AB1,AB2,AC);

	      if(conj)
		{
		  cout << endl << "These braids are conjugate." << endl << endl
		       << "A conjugating braid is: ";
		  PrintBraidWord(AC);
		  cout << endl;
		}
	      else
		{
		  cout << endl << "These braids are not conjugate." << endl;
		}

	      f.open(file);

	      f << "The braids on " << n << " strands" << endl << endl;
	      f.close();
	      PrintWord(word,n,power,file);
	      f.open(file,ios::app);

	      f << endl << endl << "and" << endl << endl;

	      f.close();
	      PrintWord(word2,n,power2,file);
	      f.open(file,ios::app);



	      f << endl << endl;

	      if(conj)
		{
		  f << "are conjugate." << endl << endl
		    << "A conjugating braid is" << endl << endl;
		  f.close();
		  PrintBraidWord(AC,file);
		}
	      else
		{
		  f << "are not conjugate.";
		  f.close();
		}

	      word.clear();
	    }

	  /////////////////////////////////////////////////////////////



	  if(c=='z')
	    {
	      n=ReadIndex();
	      word=ReadWord(n);
	      if(p)
		power=ReadPower();

	      AB=ArtinBraid(n);
	      AB=WordToBraid<ArtinPresentation>(word,n);

	      if(p && power!=1)
		AB=RaisePower(AB,power);

	      AB.MakeLCF();

	      file=ReadFileName();

	      ACent=Centralizer(AB);

	      cout << endl << "The centralizer of this braid is generated by: "
		   << endl;

	      iteration=0;

	      for(Ait=ACent.begin(); Ait!=ACent.end(); Ait++)
		{
		  cout << endl << setw(3) << ++iteration << ":  ";
		  PrintBraidWord(*Ait);
		  cout << endl;
		}

	      f.open(file);

	      f << "The contralizer of the braid on " << n << " strands"
		<< endl << endl;
	      if(p && power!=1)
		f << "( ";

	      for(itw=word.begin(); itw!=word.end(); itw++)
		{
		  if(*itw==n)
		    f << "D ";
		  else if (*itw==-n)
		    f << "-D ";
		  else
		    f << *itw << " ";
		}

	      if(p && power!=1)
		f << ")^" << power;

	      f << endl << endl << "is generated by the following braids:";
	      f.close();

	      iteration=0;
	      for(Ait=ACent.begin(); Ait!=ACent.end(); Ait++)
		{
		  f.open(file,ios::app);
		  f << endl << endl << setw(3) << ++iteration << ":  ";
		  f.close();
		  PrintBraidWord(*Ait,file);
		}

	      word.clear();
	      ACent.clear();
	    }

	  ////////////////////////////////////////////////////////

	  if(c=='s')
	    {
	      n=ReadIndex();
	      word=ReadWord(n);
	      AB=ArtinBraid(n);
	      AB=WordToBraid<ArtinPresentation>(word,n);

	      if(p)
		{
		  power=ReadPower();
		  AB=RaisePower(AB,power);
		}

	      file=ReadFileName();

	      f.open(file);

	      Asss=SSS(AB);

	      size=0;

	      for(Ait=Asss.begin(); Ait!=Asss.end(); Ait++)
		size++;

	      f << "This file contains the Super Summit Set of the braid: "
		<< endl << endl;
	      f.close();
	      PrintWord(word,n,power,file);
	      f.open(file,ios::app);


	      f << endl << endl << "It has " << size << " elements."
		<< endl << endl;

	      size=1;

	      for(Ait=Asss.begin(); Ait!=Asss.end(); Ait++)
		{
		  f << endl << setw(5) << size++;
		  f << ":   ";
		  f.close();
		  PrintBraidWord(*Ait,file);
		  f.open(file,ios::app);
		}
	      f.close();
	      word.clear();
	      Asss.clear();
	      delete[] file;
	    }

	  //////////////////////////////////////////////////////////////

	  if(c=='u')
	    {
	      n=ReadIndex();
	      word=ReadWord(n);
	      if(p)
		power=ReadPower();
	      file=ReadFileName();

	      AB=ArtinBraid(n);
	      AB=WordToBraid<ArtinPresentation>(word,n);

	      if(p && power!=1)
		AB=RaisePower(AB,power);

	      AB.MakeLCF();

	      Auss=USS(AB);

	      type=ThurstonType(Auss);

	      rigidity=Rigidity(Auss);

	      PrintUSS(Auss,word,n,power,file,type,rigidity);

	      word.clear();
	      Auss.clear();
	      delete[] file;
	    }

	  ////////////////////////////////////////////////////////////////

	  if(c=='t')
	    {
	      n=ReadIndex();
	      word=ReadWord(n);
	      if(p)
		power=ReadPower();
	      file=ReadFileName();

	      AB=ArtinBraid(n);
	      AB=WordToBraid<ArtinPresentation>(word,n);

	      if(p && power!=1)
		AB=RaisePower(AB,power);

	      AB.MakeLCF();

	      Asc=SC(AB);

	      type=ThurstonType(Asc);

	      PrintSC(Asc,word,n,power,file,type);

	      word.clear();
	      Asc.clear();
	      delete[] file;
	    }

	  ////////////////////////////////////////////////////////////////

	  if(c=='a')
	    {
	      if(p)
		{
		  p=0;
		  power=1;
		  cout << endl << "The option of taking powers is disabled."
		       << endl << endl;
		}
	      else
		{
		  p=1;
		  cout << endl << "The option of taking powers is enabled."
		       << endl << endl;
		}
	    }

	  ////////////////////////////////////////////////////////////

	  if(c=='q')
	    break;

		////////////////////////////////////////////////////////////

		if(c == 'd'){
			switch(pres){
				case Artin: pres = Band;
				case Band: pres = Artin;
			}
		}

	  ////////////////////////////////////////////////////////////

	  if(c=='?')
	    {
	      repeat=1;
	      break;
	    }

	}

      if(repeat==0)
	break;
      else
	repeat=0;

    }
  return 0;
}
