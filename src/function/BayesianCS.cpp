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
Calculate a Bayesian Score to use with Chemical Shifts CV.

The functional form of this function is
\f[
C=k_BT \sum_{i=1}^{N_{arg}} \log{ \left[ \frac{\pi}{\sqrt{2} \sigma} (x_i+2\sigma^2) \right]} + k_BT\log{\sigma}
\f]

where sigma is an uncertainty parameter,
sampled by a MC algorithm in the bounded interval specified by SIGMA_MIN and SIGMA_MAX.
The initial value of is set by SIGMA0. The MC move is a random displacement
of maximum value specified by DSIGMA.



\par Examples
The following input tells plumed to use all the HA chemical shifts with the Bayesian Score and
to print the values of the uncertainty parameter, of the MC acceptance, and of the Bayesian score.
\verbatim
WHOLEMOLECULES ENTITY0=1-174
cs:  CS2BACKBONE ATOMS=1-174 DATA=data/ FF=a03_gromacs.mdb FLAT=0.0 NRES=13 ENSEMBLE
csb: BAYESIANCS ARG=cs.ha#* SIGMA0=10.0 SIGMA_MIN=0.00001 SIGMA_MAX=10.0 DSIGMA=0.0001 KBT=2.494 MC_STEPS=10 MC_STRIDE=10 MC_SEED=1234
cse: RESTRAINT  ARG=csb SLOPE=1.0 KAPPA=0 AT=0.
PRINT ARG=csb.sigma,csb.accept,cse.bias
\endverbatim
(See also \ref CS2BACKBONE and \ref RESTRAINT).


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
  unsigned int MCseed_;
  unsigned int MCaccept_;
  Value* valueSigma;
  Value* valueAccept;
  
  void do_MC();
  double get_energy(double sigma);
  
public:
  BayesianCS(const ActionOptions&);
  void calculate();
  static void registerKeywords(Keywords& keys);
};


PLUMED_REGISTER_ACTION(BayesianCS,"BAYESIANCS")

void BayesianCS::registerKeywords(Keywords& keys){
  Function::registerKeywords(keys);
  keys.use("ARG");
  keys.add("compulsory","SIGMA0",   "10.0",   "initial value of the uncertainty parameter");
  keys.add("compulsory","SIGMA_MIN","0.00001","minimum value of the uncertainty parameter");
  keys.add("compulsory","SIGMA_MAX","10.0",   "maximum value of the uncertainty parameter");
  keys.add("compulsory","DSIGMA",   "0.0001", "maximum MC move of the uncertainty parameter");
  keys.add("compulsory","KBT",      "2.494",  "temperature of the system");
  keys.add("compulsory","MC_STEPS", "1",      "number of MC steps");
  keys.add("compulsory","MC_STRIDE","1",      "MC stride");
  keys.add("compulsory","MC_SEED",  "1234",   "MC seed");
  componentsAreNotOptional(keys); 
}

BayesianCS::BayesianCS(const ActionOptions&ao):
Action(ao),
Function(ao),
sigma0_(10.0), sigma_min_(0.00001), sigma_max_(10.0), Dsigma_(0.0001), temp_(2.494),
MCsteps_(1), MCstride_(1), MCseed_(1234), MCaccept_(0)
{
  parse("SIGMA0",   sigma0_);
  parse("SIGMA_MIN",sigma_min_);
  parse("SIGMA_MAX",sigma_max_);
  parse("DSIGMA",   Dsigma_);
  parse("KBT",      temp_);
  parse("MC_STEPS", MCsteps_);
  parse("MC_STRIDE",MCstride_);
  parse("MC_SEED",  MCseed_);
  checkRead();
 
  log.printf("  initial value of uncertainty %f\n",sigma0_);
  log.printf("  minimum value of uncertainty %f\n",sigma_min_);
  log.printf("  maximum value of uncertainty %f\n",sigma_max_);
  log.printf("  maximum MC move of the uncertainty parameter %f\n",Dsigma_);
  log.printf("  temperature of the system %f\n",temp_);
  log.printf("  number of MC steps %d\n",MCsteps_);
  log.printf("  do MC every %d steps\n", MCstride_);
  log.printf("  MC seed %d\n", MCseed_);    

  addComponent("sigma");  componentIsNotPeriodic("sigma");
  addComponent("accept"); componentIsNotPeriodic("accept");
  valueSigma=getPntrToComponent("sigma");
  valueAccept=getPntrToComponent("accept");
  
  // initialize sigma_
  sigma_ = sigma0_;
  // initialize random seed
  srand (MCseed_);
}

double BayesianCS::get_energy(double sigma){
 double ene = 0.0;
  for(unsigned i=0;i<getNumberOfArguments();++i){
    ene += std::log( getArgument(i) + 2.0 * sigma * sigma );
  }
  ene += std::log(sigma) + static_cast<double>(getNumberOfArguments())*std::log(pi/sqrt2/sigma);
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
  double new_energy = get_energy(new_sigma);
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
  
  // do MC stuff at the right time step
  long int step = getStep();
  if(step%MCstride_==0){do_MC();}

  // cycle on dataset
  double ene = 0.0;
  for(unsigned i=0;i<getNumberOfArguments();++i){
    double v = getArgument(i) + 2.0 * sigma_ * sigma_;
    // increment energy
    ene += std::log(v); 
    // set derivatives
    setDerivative(i, temp_/v);
  };
  ene += std::log(sigma_)+static_cast<double>(getNumberOfArguments())*std::log(pi/sqrt2/sigma_);
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


