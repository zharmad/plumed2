/* +++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
   Copyright (c) 2011-2014 The plumed team
   (see the PEOPLE file at the root of the distribution for a list of names)

   See http://www.plumed-code.org for more information.

   This file is part of plumed, version 2.

   plumed is free software: you can redistribute it and/or modify
   it under the terms of the GNU Lesser General Public License as published by
   the Free Software Foundation, either version 3 of the License, or
   (at your option) any later version.

   plumed is distributed in the hope that it will be useful,
   but WITHOUT ANY WARRANTY; without even the implied warranty of
   MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
   GNU Lesser General Public License for more details.

   You should have received a copy of the GNU Lesser General Public License
   along with plumed.  If not, see <http://www.gnu.org/licenses/>.
+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++ */
#include "Function.h"
#include "ActionRegister.h"

#include <cmath>

using namespace std;

namespace PLMD{
namespace function{

//+PLUMEDOC FUNCTION OUTERPRODUCT 
/*
Calculate the matrix from outer product with a set of arguments

The functional form of this function is
\f[
C[i,j]= v1_i * v1_j   
\f]


\par Examples
\endverbatim
(See also \ref PRINT and \ref DISTANCE).


*/
//+ENDPLUMEDOC


class OuterProduct :
  public Function
{
public:
  OuterProduct(const ActionOptions&);
  vector<std::string> names;	
  vector<Value*> vptr;
  void calculate();
  static void registerKeywords(Keywords& keys);
};


PLUMED_REGISTER_ACTION(OuterProduct,"OUTERPRODUCT")

void OuterProduct::registerKeywords(Keywords& keys){
  Function::registerKeywords(keys);
  keys.use("ARG"); 
}

OuterProduct::OuterProduct(const ActionOptions&ao):
Action(ao),
Function(ao)
{

  checkRead();
  for(unsigned i=0;i<getNumberOfArguments();i++){
  	for(unsigned j=i;j<getNumberOfArguments();j++){
		ostringstream oo;
		oo<<getPntrToArgument(i)->getName()<<"_"<<getPntrToArgument(j)->getName();			
		names.push_back(oo.str());
                addComponent(oo.str());componentIsNotPeriodic(oo.str());		
		vptr.push_back(getPntrToComponent(names.back()));
	}
  }

}

void OuterProduct::calculate(){
  unsigned k;k=0;
  for(unsigned i=0;i<getNumberOfArguments();i++){
  	for(unsigned j=i;j<getNumberOfArguments();j++){
		vptr[k]->set(getArgument(i)*getArgument(j));
		k++;	
	}
  }	
}

}
}


