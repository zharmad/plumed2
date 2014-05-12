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
#include "Function.h"
#include "ActionRegister.h"

#include <cmath>

using namespace std;

namespace PLMD{
namespace function{

//+PLUMEDOC FUNCTION BAYESIANCS
/*
Calculate a polynomial combination of a set of other variables.

The functional form of this function is
\f[
C=\sum_{i=1}^{N_{arg}} c_i x_i^{p_i}
\f]

The coefficients c and powers p are provided as vectors.



\par Examples
The following input tells plumed to print the distance between atoms 3 and 5
its square (as computed from the x,y,z components) and the distance
again as computed from the square root of the square.
\verbatim
DISTANCE LABEL=dist      ATOMS=3,5 COMPONENTS
COMBINE  LABEL=distance2 ARG=dist.x,dist.y,dist.z POWERS=2,2,2 PERIODIC=NO
COMBINE  LABEL=distance  ARG=distance2 POWERS=0.5 PERIODIC=NO
PRINT ARG=distance,distance2
\endverbatim
(See also \ref PRINT and \ref DISTANCE).


*/
//+ENDPLUMEDOC


class BayesianCS :
  public Function
{
  const double pi = 3.14159265359;
  const double sqrt2 = 1.41421356237;
  double sigma_;
  double sigma0_;
  double sigma_min_;
  double sigma_max_;
  double Dsigma_;
  double temp_;
  int MCsteps_;
  int MCstride_;
  int MCaccept_;
  Value* valueSigma;
  Value* valueAccept;
  
  void do_MC();
  
public:
  BayesianCS(const ActionOptions&);
  void calculate();
  static void registerKeywords(Keywords& keys);
};


PLUMED_REGISTER_ACTION(BayesianCS,"BAYESIANCS")

void BayesianCS::registerKeywords(Keywords& keys){
  Function::registerKeywords(keys);
  keys.use("ARG");
  keys.add("compulsory","SIGMA0",   "100.0", "initial value of the uncertainty parameter");
  keys.add("compulsory","SIGMA_MIN","0.0001","minimum value of the uncertainty parameter");
  keys.add("compulsory","SIGMA_MAX","100.0", "maximum value of the uncertainty parameter");
  keys.add("compulsory","DSIGMA",   "0.0001","maximum MC move of the uncertainty parameter");
  keys.add("compulsory","TEMP",     "300.0", "temperature of the system");
  keys.add("compulsory","MC_STEPS", "1",     "number of MC steps");
  keys.add("compulsory","MC_STRIDE","1",     "MC stride");
  componentsAreNotOptional(keys); 
}

BayesianCS::BayesianCS(const ActionOptions&ao):
Action(ao),
Function(ao),
sigma0_(100.0), sigma_min_(0.0001), sigma_max_(100.0), Dsigma_(0.001), temp_(300.0),
MC_steps_(1), MCstride_(1), MCaccept_(0)
{
  parse("SIGMA0",   sigma0_);
  parse("SIGMA_MIN",sigma_min_);
  parse("SIGMA_MAX",sigma_max_);
  parse("DSIGMA",   Dsigma_);
  parse("TEMP",     temp_);
  parse("MC_STEPS", MCsteps_);
  parse("MC_STRIDE",MCstride_);
  checkRead();
 
  log.printf("  initial value of uncertainty %f\n",sigma0_);
  log.printf("  minimum value of uncertainty %f\n",sigma_min_);
  log.printf("  maximum value of uncertainty %f\n",sigma_max_);
  log.printf("  maximum MC move of the uncertainty parameter %f\n",Dsigma_);
  log.printf("  temperature of the system %f\n",temp_);
  log.printf("  number of MC steps %d\n",MCsteps_);
  log.printf("  do MC every %d steps\n", MCstride_);    

  addComponent("sigma");  componentIsNotPeriodic("sigma");
  addComponent("accept"); componentIsNotPeriodic("accept");
  valueSigma=getPntrToComponent("sigma");
  valueAccept=getPntrToComponent("accept");
  
  // initialize sigma_
  sigma_ = sigma0_;
  // multiply by Boltzmann
  temp_ *= plumed.getAtoms().getKBoltzmann();
  // initialize random seed
  srand (time(NULL));
}

double BayesianCS::get_energy(double sigma){
 double ene = 0.0;
  for(unsigned i=0;i<getNumberOfArguments();++i){
    ene += log( getArgument(i) + 2.0 * sigma * sigma );
  }
  ene += log(sigma) + static_cast<double>(getNumberOfArguments())*log(pi/sqrt2/sigma);
  return temp_ * ene;
}

void BayesianCS::do_MC(){
 
 // store old energy
 double old_energy = get_energy(sigma_);
 
 for(unsigned i=0;i<MCsteps_;++i){
  // propose move
  double r = (double)rand() / RAND_MAX;
  double ds = -Dsigma_ + r * 2.0 * Dsigma_;
  double new_sigma = sigma_ + ds;
  // check boundaries
  if(new_sigma > sigma_max_){new_sigma = 2.0 * sigma_max_ - new_sigma;}
  if(new_sigma < sigma_min_){new_sigma = 2.0 * sigma_min_ - new_sigma;}
  // calculate new energy
  new_energy = get_energy(new_sigma);
  // accept or reject
  double s = (double)rand() / RAND_MAX;
  double delta = exp(-(new_energy-old_energy)/temp_);
  if(s<delta){
   old_energy = new_energy;
   sigma_ = new_sigma;
   MCaccept_++;
  }
 }
}

void BayesianCS::calculate(){
  
  // do MC stuff, if the right time step
  long int step = getStep();
  if(step%MCstride_==0){do_MC();}

  // cycle on dataset
  double ene = 0.0;
  for(unsigned i=0;i<getNumberOfArguments();++i){
    double v = getArgument(i) + 2.0 * sigma_ * sigma_;
    // increment energy
    ene += log(v); 
    // set derivatives
    setDerivative(i, temp_/v);
  };
  ene += log(sigma_)+static_cast<double>(getNumberOfArguments())*log(pi/sqrt2/sigma_);
  // set value
  setValue(temp_*ene);

  // set value of uncertainty
  valueSigma->set(sigma_);
  // calculate acceptance
  int MCtrials = step / MCstride_;
  double accept = MCaccept_ / MCsteps_ / MCtrials; 
  // set value of acceptance 
  valueAccept->set(accept);
}


}
}


