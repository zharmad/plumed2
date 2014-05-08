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
#include "Bias.h"
#include "ActionRegister.h"


using namespace std;


namespace PLMD{
namespace bias{

//+PLUMEDOC BIAS BAYESIANCS 
/*
Adds harmonic and/or linear restraints on one or more variables.  

Either or both
of SLOPE and KAPPA must be present to specify the linear and harmonic force constants
respectively.  The resulting potential is given by: 
\f[
  \sum_i \frac{k_i}{2} (x_i-a_i)^2 + m_i*(x_i-a_i)
\f].

The number of components for any vector of force constants must be equal to the number
of arguments to the action.

\par Examples
The following input tells plumed to restrain the distance between atoms 3 and 5
and the distance between atoms 2 and 4, at different equilibrium
values, and to print the energy of the restraint
\verbatim
DISTANCE ATOMS=3,5 LABEL=d1
DISTANCE ATOMS=2,4 LABEL=d2
BAYESIANCS ARG=d1,d2 AT=1.0,1.5 KAPPA=150.0,150.0 LABEL=restraint
PRINT ARG=restraint.bias
\endverbatim
(See also \ref DISTANCE and \ref PRINT).

*/
//+ENDPLUMEDOC

class BayesianCS : public Bias{
  double sigma_;
  double sigma0_;
  double sigma_min_;
  double sigma_max_;
  double Dsigma_;
  double temp_;
  int MCsteps_;
  int MCstride_;
  int MCaccept_;
  Value* valueBias;
  Value* valueForce2;
  Value* valueSigma;
  Value* valueAccept;

  void do_MC();

public:
  BayesianCS(const ActionOptions&);
  void calculate();
  static void registerKeywords(Keywords& keys);
};

PLUMED_REGISTER_ACTION(Restraint,"BayesianCS")

void BayesianCS::registerKeywords(Keywords& keys){
   Bias::registerKeywords(keys);
   keys.use("ARG");
   keys.add("compulsory","SIGMA0",   "100.0", "initial value of the uncertainty parameter");
   keys.add("compulsory","SIGMA_MIN","0.0001","minimum value of the uncertainty parameter");
   keys.add("compulsory","SIGMA_MAX","100.0", "maximum value of the uncertainty parameter");
   keys.add("compulsory","DSIGMA",   "0.0001","maximum MC move of the uncertainty parameter");
   keys.add("compulsory","TEMP",     "300.0", "temperature of the system");
   keys.add("compulsory","MC_STEPS", "1",     "number of MC steps");
   keys.add("compulsory","MC_STRIDE","1",     "MC stride");
   componentsAreNotOptional(keys);
   keys.addOutputComponent("bias","default","the instantaneous value of the bias potential");
   keys.addOutputComponent("force2","default","the instantaneous value of the squared force due to this bias potential");
   keys.addOutputComponent("sigma","default","the instantaneous value of the uncertainty parameter");
   keys.addOutputComponent("accept","default","MC acceptance");
}

BayesianCS::BayesianCS(const ActionOptions&ao):
PLUMED_BIAS_INIT(ao),
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

  addComponent("bias");   componentIsNotPeriodic("bias");
  addComponent("force2"); componentIsNotPeriodic("force2");
  addComponent("sigma");  componentIsNotPeriodic("sigma");
  addComponent("accept"); componentIsNotPeriodic("accept");
  valueBias=getPntrToComponent("bias");
  valueForce2=getPntrToComponent("force2");
  valueSigma=getPntrToComponent("sigma");
  valueAccept=getPntrToComponent("accept");
  
  // initialize sigma_
  sigma_ = sigma0_;
  // multiply by Boltzmann
  temp_ *= plumed.getAtoms().getKBoltzmann();
  // initialize random seed
  srand (time(NULL));
}

void BayesianCS::get_energy(double sigma){

}

void BayesianCS::get_force(double sigma){

}

void BayesianCS::do_MC(){
 
 double old_energy = get_energy(sigma_);
 
 for(unsigned i=0;i<MCsteps_;++i){
  // propose move
  double r = (double)rand() / RAND_MAX;
  double ds = -Dsigma_ + r * 2.0 * Dsigma_;
  double new_sigma = sigma_ + ds;
  // check boundaries
  if(new_sigma > sigma_max_){new_sigma = 2.0 * sigma_max_ - new_sigma;}
  if(new_sigma < sigma_min_){new_sigma = 2.0 * sigma_min_ - new_sigma;}
  // calculate energy
  new_energy = get_energy(new_sigma);
  // accept or reject 
  if(accepted){
   old_energy = new_energy;
   sigma_ = new_sigma;
   MCaccept_++;
  }
 }
}

void BayesianCS::calculate(){
  double ene=0.0;
  double totf2=0.0;
  
  // do MC stuff, if the right time step
  long int step = getStep();
  if(step%MCstride_==0){do_MC();}

  // calculate restraint
  for(unsigned i=0;i<getNumberOfArguments();++i){
    // cycle on components of i-th argument

    const double cv=difference(i,at[i],getArgument(i));
    const double k=kappa[i];
    const double m=slope[i];
    const double f=-(k*cv+m);
    ene+=0.5*k*cv*cv+m*cv;
    setOutputForce(i,f);
    totf2+=f*f;

  };

  valueBias->set(ene);
  valueForce2->set(totf2);
  valueSigma->set(sigma_);
  // calculate acceptance
  int MCtrials = step / MCstride_;
  double accept = MCaccept_ / MCsteps_ / MCtrials;  
  valueAccept->set(accept);
}

}


}
