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
#include "SMACOF.h"
#include "reference/PointWiseMapping.h"

namespace PLMD {
namespace analysis {

void SMACOF::run( const Matrix<double>& Weights, PointWiseMapping* mymap ){
    Matrix<double> Distances( mymap->modifyDmat() ); unsigned M = Distances.nrows();
    Matrix<double> InitialZ( M, mymap->getNumberOfProperties() );     //want Mx2 Z matrix but you set columns to M dimension first and then changes to 2 underneath
    for(unsigned i=0; i<M; ++i){
       for(unsigned j=0; j<mymap->getNumberOfProperties(); ++j){
            InitialZ(i,j) =  mymap->getProjectionCoordinate( i, j );    
       }
    }
   
    // Calculate V
    Matrix<double> V(M,M);
    for(unsigned i=0; i<M; ++i){
      for(unsigned j=0; j<M; ++j){
         if(i==j) continue;
         V(i,j)=-Weights(i,j);
      }
      for(unsigned j=0; j<M; ++j){
        if(i==j)continue;
        V(i,i)-=V(i,j);
      }   
    }    
    
    // And invert V
    Matrix<double> mypseudo(M, M);
    pseudoInvert(V, mypseudo);
    double myfirstsig = calculateSigma( Weights, Distances, InitialZ );    //this is a function that outputs the initial sigma values
    // initial sigma is made up of the original distances minus the distances between the projections all squared.
    unsigned MAXSTEPS=1000; double tol=1.E-2; Matrix<double> temp( M, M ), BZ( M, M ), newZ( M, mymap->getNumberOfProperties() );
    for(unsigned i=0;i<MAXSTEPS;++i){
        if(i==MAXSTEPS-1) plumed_merror("ran out of steps in SMACOF algorithm");
//        printf("Doing step %d of SMACOF %d %d \n",i,M,mymap->getNumberOfProperties() );        

        // Calculate B(Z)  which is a M x M square matrix
        //Off diagonals considered first 
        // Matrix<double> BZ(M,M); //create the space memory for the MxM matrix
        for(unsigned i=0; i<M; ++i){
           for(unsigned j=0; j<M; ++j){
            if(i==j) continue;  //skips over the diagonal elements
            //Calculate the distance between point i and point j.
            //Create the space memory for a scalar, as this is what the distance is.
            double dist=0;
            for(unsigned k=0; k<mymap->getNumberOfProperties(); ++k){  //mymap->getNumberOfProperties() is the range that the columns in InitialZ go up to therefore we use it here also
              double tmp=InitialZ(i,k) - InitialZ(j,k);
              dist+=tmp*tmp; //this is squaring the line above. Quicker way to do things than to actually write tmp^2
            } 
            //off diagonal elements are [-1x(dissimilarities/distances)]/M i.e. BZ isn't symmetric
            //If it was symmetric that line would simply read BZ(i,j)=BZ(j,i)=sqrt(dist)
            BZ(i,j)=-(Weights(i,j)*Distances(i,j) / sqrt(dist)); 
            //static_cast<double>(M) is the M that we want to dive by in the code 
           }
           //the diagonal elements are -off diagonal elements BZ(i,i)-=BZ(i,j)   (Equation 8.25)
           BZ(i,i)=0; //create the space memory for the diagonal elements which are scalars
           for(unsigned j=0; j<M; ++j){ 
              if(i==j) continue;
              BZ(i,i)-=BZ(i,j); 
           }
        }
        
        // Matrix matrix multiply B(Z) times Z and multiply by 1/M. This is the Guttman transform for w(ij)=1, which calculates X which we know is exactly equal to the newZ value
        mult( mypseudo, BZ, temp); mult(temp, InitialZ, newZ);   
        //Compute new sigma
        double newsig = calculateSigma( Weights, Distances, newZ );
        printf("CHECKING VALUE OF SIGMA %d %f %f %f \n",i,myfirstsig,newsig,myfirstsig-newsig);
        //Computing whether the algorithm has converged (has the mass of the potato change
        //when we put it back in the oven!)
        if( fabs( newsig - myfirstsig )<tol ) break;    //An f means real and by convention is put in before abs. abs means absolute value.
        // Make initial sigma into new sigma so that the value of new sigma is used every time so that the error can be reduced
        myfirstsig=newsig;       
        // Set InitialZ equal to newZ (step 7)
        InitialZ = newZ;
    } 

    // Pass final projections to map object
    for(unsigned i=0;i<M;++i){
        for(unsigned j=0;j<mymap->getNumberOfProperties();++j) mymap->setProjectionCoordinate( i, j, newZ(i,j) ); 
        //is the line above right? i havent changed it but i dont think it is right.
    }

}

double SMACOF::calculateSigma( const Matrix<double>& Weights, const Matrix<double>& Distances, const Matrix<double>& InitialZ ){    //& sign always put in a function as routine but dont need to know why
    unsigned M = Distances.nrows();
    double sigma=0;  //Sigma is a scalar hence why it needs to be written in this format
    for(unsigned i=1; i<M; ++i){
       for(unsigned j=0; j<i; ++j){
       double dlow=0;  //dlow are the distances in the lower triangular part of the matrix. Don't need to include all distances because symmentric
          for(unsigned k=0; k<InitialZ.ncols();++k){
          //this is doing the summation again with ... being the highest value in the summation.
          double tmp2=InitialZ(i,k) - InitialZ(j,k); //this is the maths within the summation sign
          dlow+=tmp2*tmp2;
          }
       double tmp3= Distances(i,j) - sqrt(dlow);
       sigma+= Weights(i,j)*tmp3*tmp3;         
       }
    }
    return sigma;     //Returns the value you have specified in the function- a function always needs to return something
}      

}
}
