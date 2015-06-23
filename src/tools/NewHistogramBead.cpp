/* +++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
   Copyright (c) 2012-2014 The plumed team
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
#include "NewHistogramBead.h"
#include <vector>
#include <limits>
#include "Tools.h"
#include "Keywords.h"

namespace PLMD{

//+PLUMEDOC INTERNAL histogrambead 
/*
A function that can be used to calculate whether quantities are between fixed upper and lower bounds.

If we have multiple instances of a variable we can estimate the probability distribution (pdf)
for that variable using a process called kernel density estimation:

\f[
P(s) = \sum_i K\left( \frac{s - s_i}{w} \right)
\f] 

In this equation \f$K\f$ is a symmetric funciton that must integrate to one that is often
called a kernel function and \f$w\f$ is a smearing parameter.  From a pdf calculated using 
kernel density estimation we can calculate the number/fraction of values between an upper and lower
bound using:

\f[
w(s) = \int_a^b \sum_i K\left( \frac{s - s_i}{w} \right)
\f]     

All the input to calculate a quantity like \f$w(s)\f$ is generally provided through a single 
keyword that will have the following form:

KEYWORD={TYPE UPPER=\f$a\f$ LOWER=\f$b\f$ SMEAR=\f$\frac{w}{b-a}\f$}

This will calculate the number of values between \f$a\f$ and \f$b\f$.  To calculate
the fraction of values you add the word NORM to the input specification.  If the 
function keyword SMEAR is not present \f$w\f$ is set equal to \f$0.5(b-a)\f$. Finally,
type should specify one of the kernel types that is present in plumed. These are listed
in the table below:

<table align=center frame=void width=95%% cellpadding=5%%>
<tr> 
<td> TYPE </td> <td> FUNCTION </td> 
</tr> <tr> 
<td> GAUSSIAN </td> <td> \f$\frac{1}{\sqrt{2\pi}w} \exp\left( -\frac{(s-s_i)^2}{2w^2} \right)\f$ </td>
</tr> <tr>
<td> TRIANGULAR </td> <td> \f$ \frac{1}{2w} \left( 1. - \left| \frac{s-s_i}{w} \right| \right) \quad \frac{s-s_i}{w}<1 \f$ </td>
</tr>
</table>

Some keywords can also be used to calculate a descretized version of the histogram.  That
is to say the number of values between \f$a\f$ and \f$b\f$, the number of values between
\f$b\f$ and \f$c\f$ and so on.  A keyword that specifies this sort of calculation would look
something like

KEYWORD={TYPE UPPER=\f$a\f$ LOWER=\f$b\f$ NBINS=\f$n\f$ SMEAR=\f$\frac{w}{n(b-a)}\f$}
 
This specification would calculate the following vector of quantities: 

\f[
w_j(s) = \int_{a + \frac{j-1}{n}(b-a)}^{a + \frac{j}{n}(b-a)} \sum_i K\left( \frac{s - s_i}{w} \right) 
\f]

*/
//+ENDPLUMEDOC

void NewHistogramBead::registerKeywords( Keywords& keys ){
  SwitchingFunction::registerKeywords( keys ); keys.remove("D_0");
}

NewHistogramBead::NewHistogramBead():
init(false),
periodicity(unset),
min(0.0),
max(0.0),
max_minus_min(0.0),
inv_max_minus_min(0.0)
{
}

std::string NewHistogramBead::description() const {
  return sf.description();
}

void NewHistogramBead::set( const std::string& params, std::string& errormsg ){
  std::vector<std::string> data=Tools::getWords(params);
  if(data.size()<1) errormsg="No input has been specified";

  if( data[0]=="GAUSSIAN" ) plumed_merror("gaussian kernels are now specified using the CERF keyword"); 

  std::string sw_input;
  for(unsigned i=0;i<data.size();++i){
      if( data[i].find("D_0")!=std::string::npos){
          errormsg="D_0 keyword should not be present in histogram bead input";
          return;
      }
      sw_input += " " + data[i];
  }
  sf.set(sw_input,errormsg); init=true;  
}

double NewHistogramBead::setBoundsAndCalculate( const double& lowb, const double& highb, const double& x, double& df ) const {
  plumed_dbg_assert(init && periodicity!=unset ); 
  double lowB, upperB, f;
  double df1, f1, diff1 = difference( x, lowb );
  if( diff1>0 ) f1 = 1. - sf.calculate( diff1, df1 ); 
  else f1 = -(1. - sf.calculate( fabs(diff1), df1 ) );
  if( fabs(diff1)>0 ) df1*=fabs(diff1);

  double df2, f2, diff2 = difference( x, highb );
  if( diff2>0 ) f2 = 1. - sf.calculate( diff2, df2 ); 
  else f2 = -(1. - sf.calculate( fabs(diff2), df2 ));
  if( fabs(diff2)>0 ) df2*=fabs(diff2);

  df = 0.5*(df2 - df1);
  return 0.5*(f2-f1);
}

double NewHistogramBead::boundDerivative( const double& highb, const double& x ) const { 
  plumed_dbg_assert(init && periodicity!=unset );
  double df2, f2, diff2 = difference( x, highb );
  if( diff2>0 ) f2 = 1. - sf.calculate( diff2, df2 );
  else f2 = -(1. - sf.calculate( fabs(diff2), df2 ));
  return df2*fabs(diff2); 
}

}
