/* +++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
   Copyright (c) 2012 The plumed team
   (see the PEOPLE file at the root of the distribution for a list of names)

   See http://www.plumed-code.org for more information.

   This file is part of plumed, version 2.0.

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
#include "AnalysisWithLandmarks.h"
#include "ClassicalScaling.h"
#include "SMACOF.h"
#include "reference/PointWiseMapping.h"
#include "core/ActionRegister.h"
#include "tools/SwitchingFunction.h"

namespace PLMD {
namespace analysis {

class SketchMap : public AnalysisWithLandmarks {
private:
  unsigned nlow;
  bool nosmacof;
  std::string ofilename;
  std::string efilename;
  PointWiseMapping* myembedding;
  SwitchingFunction lowdf, highdf;
  double recalculateWeights( const Matrix<double>& Distances, const Matrix<double>& F, Matrix<double>& weights );
public:
  static void registerKeywords( Keywords& keys );
  SketchMap( const ActionOptions& ao );
  ~SketchMap();
  void analyzeLandmarks();
};

PLUMED_REGISTER_ACTION(SketchMap,"SketchMap")

void SketchMap::registerKeywords( Keywords& keys ){
  AnalysisWithLandmarks::registerKeywords( keys );
  keys.add("compulsory","NLOW_DIM","number of low-dimensional coordinates required");
  keys.add("compulsory","HIGH_DIM_FUNCTION","the parameters of the switching function in the high dimensional space");
  keys.add("compulsory","LOW_DIM_FUNCTION","the parameters of the switching function in the low dimensional space");
  keys.add("compulsory","OUTPUT_FILE","file on which to output the final embedding coordinates");
  keys.add("compulsory","EMBEDDING_OFILE","dont output","file on which to output the embedding in plumed input format");
}

SketchMap::SketchMap( const ActionOptions& ao ):
Action(ao),
AnalysisWithLandmarks(ao)
{
  myembedding = new PointWiseMapping( getMetricName(), false );
  setDataToAnalyze( dynamic_cast<MultiReferenceBase*>(myembedding) );

  parse("NLOW_DIM",nlow); 
  if( nlow<1 ) error("dimensionality of low dimensional space must be at least one");

  // Read in the switching functions
  std::string linput,hinput, errors;
  parse("HIGH_DIM_FUNCTION",hinput);
  highdf.set(hinput,errors);
  if(errors.length()>0) error(errors);
  parse("LOW_DIM_FUNCTION",linput);
  lowdf.set(hinput,errors);
  if(errors.length()>0) error(errors);

  // Setup the property names
  std::vector<std::string> propnames( nlow ); std::string num;
  for(unsigned i=0;i<propnames.size();++i){
     Tools::convert(i+1,num); std::string lab=getLabel();
     if(lab.find("@")!=std::string::npos) propnames[i]=getName() + "." + num;
     else propnames[i]=getLabel() + "." + num;
  }
  myembedding->setPropertyNames( propnames, false );

  parseOutputFile("EMBEDDING_OFILE",efilename);
  parseOutputFile("OUTPUT_FILE",ofilename);
}

SketchMap::~SketchMap(){
  delete myembedding;
}

void SketchMap::analyzeLandmarks(){


   // The new algorithm goes in here.  To calculate the switching functions use the following code:
   //
   // For the switching function on the distances in the high dimensional space, F(D_ij), use:
   //
   // Fr = highdf.calculate( r, dr )
   //
   // where fr is the value of the value of the function, r is the value of the distance and r*dr is the
   // first derivative of the function.
   //
   // For the switching funciton on the distances in the low dimensional space, f(d_ij), use:
   //
   // fr = lowdf.calculate( r, dr )

  //Inputting the weights
    unsigned M = myembedding->getNumberOfReferenceFrames();
    Matrix<double> NiNj(M, M); NiNj=0.0;
    for(unsigned i=0; i<M; ++i){
        for(unsigned j=0; j<M; ++j){
            double n_i,n_j=0;
            n_i=myembedding->getWeight(i);
            n_j=myembedding->getWeight(j);
            if(i==j) n_i = n_j = 0;
            NiNj(i,j) = n_i * n_j;
                        
        }
    }
    
    // Calculate matrix of dissimilarities (High dimensional space) 
    myembedding->calculateAllDistances( getPbc(), getArguments(), comm, myembedding->modifyDmat(), true );
    Matrix<double> Distances( mymap->modifyDmat() ); //making a copy
    Matrix<double> F(M,M); double dr;
    for(unsigned i=1; i<M; ++i){
        for(unsigned j=0; j<i; ++j) {
            F(i,j) = F(j,i) = highdf.calculate( Distances(i,j), dr ); // high dim space
        }
        
    }
    // Calculates the first guess of projections in LD space
    ClassicalScaling::run( myembedding );
    
    
    // Calculate the value of sigma and the weights
    Matrix<double> weights(M,M); double filt = recalculateWeights( Distances, F, weights );

    unsigned MAXSTEPS=100; double newsig;
    for(unsigned i=0;i<MAXSTEPS;++i){
        // Run the smacof algorithm
        SMACOF::run( weights, myembedding );
        // Recalculate weights matrix and sigma
        newsig = recalculateWeights( Distances, F, weights );       
        // Test whether or not the algorithm has converged
        if( fabs( newsig - filt )<tol ) break;
        // Make initial sigma into new sigma so that the value of new sigma is used every time so that the error can be reduced
        filt=newsig;
    } 
          
//   //Computing whether the algorithm has converged (has the mass of the potato change
//     //when we put it back in the oven!
//     if( fabs( newsig - myfirstsig )<tol ) break;    //An f means real and by convention is put in before abs. abs means absolute value.
//     // Make initial sigma into new sigma so that the value of new sigma is used every time so that the error can be reduced
//     myfirstsig=newsig;
//     // Set InitialZ equal to newZ (step 7)
//     InitialZ = newZ;
//         
//         
// 
//     //need to call in dissimilarities and F(Dij)
//     Fr = highdf.calculate( r, dr ) // high dim space
//     
//     
//     // Computing the switching functions
//     // r is the distance value, dr= the derivative and df(r)/dr x 1/r 
//     double r=0; //r is the distance value
//     fr = lowdf.calculate( r, dr ) // low dim space
//     
// 
//     
//     
//     
//     
//     Fr = highdf.calculate( r, dr ) // high dim space
//     fr = lowdf.calculate( r, dr ) // low dim space


  // Output the embedding as long lists of data
  OFile gfile; gfile.link(*this); 
  gfile.setBackupString("analysis");
  gfile.fmtField(getOutputFormat()+" ");
  gfile.open( ofilename.c_str() );
  
  // Print embedding coordinates
  for(unsigned i=0;i<myembedding->getNumberOfReferenceFrames();++i){
      for(unsigned j=0;j<nlow;++j){
          std::string num; Tools::convert(j+1,num);
          gfile.printField( getLabel() + "." + num , myembedding->getProjectionCoordinate(i,j) );
      }
      gfile.printField();
  }  
  gfile.close();

  // Output the embedding in plumed format
  if( efilename!="dont output"){
     // std::string ifname=saveResultsFromPreviousAnalyses( efilename );
     OFile afile; afile.link(*this); afile.setBackupString("analysis");
     afile.open( efilename.c_str() );
     myembedding->print( "classical mds", getTime(), afile, getOutputFormat() );
     afile.close();
  }
}

double SketchMap::recalculateWeights( const Matrix<double>& Distances, const Matrix<double>& F, Matrix<double>& weights ){
  double filt=0; Matrix<double> weights(M,M);
  for(unsigned i=1; i<M; ++i){
      for(unsigned j=0; j<i; ++j){
          double tempd=0;
          for(unsigned k=0;k<nlow;++k){
             double tmp = myembedding->getProjectionCoordinate( i, k ) - myembedding->getProjectionCoordinate( j, k );
             tempd += tmp*tmp;
          }
          double dij=sqrt(tempd);
          double fij=lowdf.calculate( dij, dr );
          double filter=F(i,j)-fij;
          weights(i,j)=weights(j,i) = ( filter*dij*dr ) / ( Distances(i,j) - dij );
          filt+=NiNj(i,j)*filter*filter;
      }
  }
  return filt;
}

}
}
