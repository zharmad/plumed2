/* +++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
   Copyright (c) 2013 The plumed team
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
#include "Analysis.h"
#include "core/ActionRegister.h"
#include "core/ActionWithValue.h"
#include <cmath>

namespace PLMD{
namespace analysis{

//+PLUMEDOC ANALYSIS AVERAGE 
/* 
Calculate the average and stdev of some values calculated in other functions 

\par Examples

The following input monitors two torsional angles during a simulation
and outputs a histogram as a function of them at the end of the simulation.
\verbatim
TORSION ATOMS=1,2,3,4 LABEL=r1
TORSION ATOMS=2,3,4,5 LABEL=r2
AVERAGE ...
  ARG=r1,r2 
  STRIDE=10
... AVERAGE 
\endverbatim


*/
//+ENDPLUMEDOC

class Average : public Analysis, public ActionWithValue {
private:
  unsigned nsample;
  std::vector<std::string> avgvec_name;
  std::vector<std::string> stdevvec_name;
  std::vector<double> avgvec;
  std::vector<double> stdevvec;
  bool hasStdev;
public:
  static void registerKeywords( Keywords& keys );
  Average(const ActionOptions&ao);
  void performAnalysis();
  void performTask(){plumed_error();}; // this is required since it inherits from and ActionWithVessels 
  void intermediateAnalysis();
  unsigned getNumberOfDerivatives(){return 0;}; // this is required since it inherits from ActionWithValues
  void updateAverage();
};

PLUMED_REGISTER_ACTION(Average,"AVERAGE")

void Average::registerKeywords( Keywords& keys ){
  Analysis::registerKeywords( keys );
  keys.remove("ATOMS"); // the base class require ATOMS: here we do not need 
  keys.reset_style("ARG","compulsory"); // the base does not require ARG but we do
  keys.addFlag("STDDEV",false," This includes the standard deviation ");
  ActionWithValue::useCustomisableComponents(keys);
}

Average::Average(const ActionOptions&ao):
PLUMED_ANALYSIS_INIT(ao),ActionWithValue(ao),hasStdev(false)
{
  nsample=0;
  parseFlag("STDDEV",hasStdev);  
  checkRead();
  // add as many values as arguments
  for(unsigned i=0;i<getNumberOfArguments();++i){
	// the periodicity
        std::string myname; myname=getPntrToArgument(i)->getName();	
	ActionWithValue::addComponent(myname+"_avg");avgvec_name.push_back(myname+"_avg");avgvec.push_back(0.); 
	if(hasStdev)ActionWithValue::addComponent(myname+"_stdev");stdevvec_name.push_back(myname+"_stdev");stdevvec.push_back(0.);
	bool isperiodic=getPntrToArgument(i)->isPeriodic();
	if(isperiodic){
		plumed_merror("Average does not support periodic calculations\n");	
                // TODO: for each value use two compoents rescaled onto a -pi:pi range
		// then do the avgs of sin/cos, at the end do the tan  and rescale back onto the original range 
	}else{
		componentIsNotPeriodic(myname+"_avg");
		if(hasStdev)componentIsNotPeriodic(myname+"_stdev");
	}
  }
  
}
// at this stage perform a final analysis: will not work if the calculation is killed 
void Average::performAnalysis(){
	updateAverage();	
}
// every STRIDE
void Average::intermediateAnalysis(){
  updateAverage();
  // accumulate  

}

void Average::updateAverage(){
  // calculate should have already fetched everything from the arguments: it is in the   
  //std::vector<ReferenceConfiguration*> data 
  std::vector<double> point( getNumberOfArguments() );
  double weight;weight=1.;
  getDataPoint( getNumberOfDataPoints()-1, point, weight );
  nsample=getNumberOfDataPoints();
  for(unsigned i=0;i<point.size();++i){
       avgvec[i]+=point[i];
       getPntrToComponent(avgvec_name[i])->set(avgvec[i]/float(nsample));
  }
  if(hasStdev){
	  for(unsigned i=0;i<point.size();++i){
	        stdevvec[i]+=point[i]*point[i];// accumulate the squares
	        getPntrToComponent(stdevvec_name[i])->set(sqrt(stdevvec[i]/float(nsample)-  getOutputQuantity(2*i)*getOutputQuantity(2*i)));
	  }
  }

}

}
}
