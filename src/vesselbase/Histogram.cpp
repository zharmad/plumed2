/* +++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
   Copyright (c) 2013,2014 The plumed team
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
#include "VesselRegister.h"
#include "ShortcutVessel.h"

namespace PLMD{
namespace vesselbase{

class Histogram : public ShortcutVessel {
public:
  static void registerKeywords( Keywords& keys );
  static void reserveKeyword( Keywords& keys );
  Histogram( const VesselOptions& da );
};

PLUMED_REGISTER_VESSEL(Histogram,"HISTOGRAM")

void Histogram::registerKeywords( Keywords& keys ){
  ShortcutVessel::registerKeywords( keys );
  keys.add("compulsory","LOWER","the lower boundary on the range of interest");
  keys.add("compulsory","UPPER","the upper boundary on the range of interest");
  keys.add("compulsory","NBINS","The number of equal width bins you want to divide the range into");
  keys.addFlag("NORM",false,"calculate the fraction of values rather than the number"); 
}

void Histogram::reserveKeyword( Keywords& keys ){
  keys.reserve("optional","HISTOGRAM","calculate a discretized histogram of the distribution of values. "
                                      "This shortcut allows you to calculates NBIN quantites like BETWEEN.");
}

Histogram::Histogram( const VesselOptions& da ):
ShortcutVessel(da)
{
  bool norm; parseFlag("NORM",norm); std::string normstr="";
  if(norm) normstr=" NORM";

  unsigned nbins; parse("NBINS",nbins);
  double min; parse("LOWER",min);
  double max; parse("UPPER",max);  

  std::string instr = getAllInput();
  double delr = ( max-min ) / static_cast<double>( nbins );
  for(unsigned i=0;i<nbins;++i){
      std::string lb, ub, bind;
      Tools::convert( min+i*delr, lb );
      Tools::convert( min+(i+1)*delr, ub );
      bind = instr + " " + "LOWER=" + lb + " " + "UPPER=" + ub + normstr;
      addVessel("BETWEEN",bind);
  }
}

}
}
