/* +++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
   Copyright (c) 2014 The plumed team
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
#include "core/ActionRegister.h"
#include "tools/NewHistogramBead.h"
#include "MultiColvarFilter.h"

//+PLUMEDOC MCOLVARF MFILTER_BETWEEN
/*
This action can be used to filter the distribution of colvar values in a multicolvar 
so that one can compute the mean and so on for only those multicolvars within a certain range.

\par Examples

*/
//+ENDPLUMEDOC

namespace PLMD {
namespace multicolvar {

class FilterBetween : public MultiColvarFilter {
private:
  double lowb, highb;
  NewHistogramBead hb;
public:
  static void registerKeywords( Keywords& keys );
  FilterBetween(const ActionOptions& ao);
  double applyFilter( const double& val, double& df ) const ;
}; 

PLUMED_REGISTER_ACTION(FilterBetween,"MFILTER_BETWEEN")

void FilterBetween::registerKeywords( Keywords& keys ){
  MultiColvarFilter::registerKeywords( keys );
  keys.add("compulsory","LOWER","the lower boundary for the range of interest");
  keys.add("compulsory","UPPER","the upper boundary for the range of interest");
  keys.add("compulsory","SWITCH","the switching function to use around the region of interest. "
                                 "For more information see \\ref histogrambead");
}

FilterBetween::FilterBetween(const ActionOptions& ao):
Action(ao),
MultiColvarFilter(ao)
{
  // Read in the switching function
  std::string sinput, errors; 
  parse("SWITCH",sinput); hb.set(sinput,errors);
  if( errors.length()>0 ) error("problem reading SWITCH keyword : " + errors);
  parse("LOWER",lowb); parse("UPPER",highb);

  if( getPntrToMultiColvar()->isPeriodic() ){
     std::string min, max; getPntrToMultiColvar()->retrieveDomain( min, max );
     double mlow, mhigh; Tools::convert( min,mlow ); Tools::convert( max, mhigh);
     hb.isPeriodic( mlow, mhigh );
  } else {
     hb.isNotPeriodic();
  }
  log.printf("  filtering colvar values and focussing only on those between %d and %d \n",lowb, highb);
  log.printf("  switching functions have a width of %s \n", hb.description().c_str() );

  checkRead();  
}

double FilterBetween::applyFilter( const double& val, double& df ) const {
  double f = hb.setBoundsAndCalculate( lowb, highb, val, df ); 
  return f;
}

}
}
