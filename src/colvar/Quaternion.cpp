/* +++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
   Copyright (c) 2011-2015 The plumed team
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
#include "Colvar.h"
#include "core/PlumedMain.h"
#include "ActionRegister.h"
#include "tools/PDB.h"
#include "reference/RMSDBase.h"
#include "reference/MetricRegister.h"
#include "core/Atoms.h"


#include <string>
#include <cmath>

using namespace std;

namespace PLMD{
namespace colvar{

//+PLUMEDOC COLVAR QUATERNION
/*
 * Calculates quaternion rotation to a reference.
 * Version: 0.1 - Poker Chen 12.11.2015
 * This functionality is intended as a near-clone to RMSD,
 * in which one can restrain the orientation of a molecule in an arbitrary
 * rotation relative to a reference frame, or relative orientation of
 * two molecules, given respective reference coordinates that
 * would place them in the same frame.
 * 
 * There are several KEYWORDS that are necessary. One must have wither one or two
 * reference coordinate files, with atoms to be fitted having an occupancy of > 0.0.
 * If only one is given in REFERENCE, then the absolute rotation with the system frame
 * is calculated.
 * If both REFERENCE and REFERENCE_B are given, then the relative orientation
 * between the two groups are given is calculated, using the convention
 *  qB = qr.qA, or, qr=qB.qA^-1 = qB.(-qA)
 * 
\par Examples

<!---You should put an example of how to use your CV here--->

*/
//+ENDPLUMEDOC

//Declare class and content variables
class Quaternion : public Colvar {
    private:
    //bool pbc; No PBC.
    //MultiValue posder1; //a flexible vector to store positions and derivatives.
    //MultiValue posder2; //a flexible vector to store positions and derivatives.
    //ReferenceValuePack mypack1; //Object to carry derivatives from calculations to system.
    //ReferenceValuePack mypack2;
    PLMD::RMSDBase* rmsd1;
    PLMD::RMSDBase* rmsd2;
    std::vector<Vector> refpos1;
    std::vector<Vector> refpos2;
    std::vector<double> refw1;
    std::vector<double> refw2;
    std::vector<double> refdisp1;
    std::vector<double> refdisp2;
    Vector refc1;
    bool bRefc1_is_calculated;
    bool bRef1_is_centered;
    Vector refc2;
    bool bRefc2_is_calculated;
    bool bRef2_is_centered;
    
    // Relevant variables for outputs and derivatives.
    double dist;
    double lambda;
    Vector4d quat;
    Vector4d normquat;
    std::vector<double> nqdummy;
    //std::vector<double> eigenvals;
    //Matrix<double> eigenvecs;
    double rr00; //  sum of positions squared (needed for dist calc)
    double rr11; //  sum of reference squared (needed for dist calc)
    Tensor rr01;
    Tensor rotation; // rotation derived from the eigenvector having the smallest eigenvalue
    Tensor drotation_drr01[3][3]; // derivative of the rotation only available when align!=displace
    Tensor ddist_drr01;
    Tensor ddist_drotation;
    std::vector<Vector> diff; // difference of components

public:
  explicit Quaternion(const ActionOptions&);
  ~Quaternion();
  void clear();
  // set reference coordinates, remove the com by using uniform weights
  void setReference1(const std::vector<Vector> & inppos, const std::vector<double> & inpw, const std::vector<double> & inpdisp );
  void setReference2(const std::vector<Vector> & inppos, const std::vector<double> & inpw, const std::vector<double> & inpdisp );
  void setReference1(const PDB inppdb );
  void setReference2(const PDB inppdb ); 
// active methods:
  //Because I don't have the time to fuck around with parseVector from Bias.
  void stringToQuat(const string str, std::vector<double> q, const char *fs);
  void optimalAlignment1(const std::vector<Vector> & currpos);
  virtual void calculate();
  static void registerKeywords(Keywords& keys);
};

//Functions to parse options and keywords in plumed.dat
PLUMED_REGISTER_ACTION(Quaternion,"QUATERNION")

void Quaternion::registerKeywords(Keywords& keys){
  Colvar::registerKeywords(keys);
  keys.add("compulsory","REFERENCE","A file in pdb format containing the reference structure and the atoms involved in the CV.");  
  keys.add("optional","REFERENCE_B","A second reference file for the second group to calculate relative orientation.");
  keys.add("compulsory","NORM_DIRECTION","w","q-space is double defined such that q=-q, so it is conventional to define an alignment direction."
                      "This defaults to (1,0,0,0) in general literature, but should be assigned to a similar direction to a restraint."
                      "Options: x=(0,1,0,0), y, z, or an arbitrary quaternion.");
  //This must always be optimal rotations.
  //keys.add("compulsory","TYPE","SIMPLE","the manner in which RMSD alignment is performed.  Should be OPTIMAL or SIMPLE.");
  //keys.addFlag("SQUARED",false," This should be setted if you want MSD instead of RMSD ");
  keys.remove("NOPBC");
  keys.addFlag("COMPONENTS",true,"(Compulsary) Calculate the quaternion as 4-individual colvars, stored as label.w, label.x, label.y and label.z");    
  keys.addOutputComponent("w","COMPONENTS","the w-component of the rotational quaternion");  
  keys.addOutputComponent("x","COMPONENTS","the x-component of the rotational quaternion");
  keys.addOutputComponent("y","COMPONENTS","the y-component of the rotational quaternion");
  keys.addOutputComponent("z","COMPONENTS","the z-component of the rotational quaternion");  
  //keys.addFlag("TEMPLATE_DEFAULT_OFF_FLAG",false,"flags that are by default not performed should be specified like this");
  //keys.addFlag("TEMPLATE_DEFAULT_ON_FLAG",true,"flags that are by default performed should be specified like this");
  //keys.add("compulsory","TEMPLATE_COMPULSORY","all compulsory keywords should be added like this with a description here");
  //keys.add("optional","TEMPLATE_OPTIONAL","all optional keywords that have input should be added like a description here");
  //keys.add("atoms","TEMPLATE_INPUT","the keyword with which you specify what atoms to use should be added like this");
}

void Quaternion::stringToQuat(const string str, std::vector<double> q, const char* fs)
{
    //Use find to location positions of field separators
    //Cat each substring and give back q.
    //str.substr(0,3);
    //q[0];
}

//Constructor
Quaternion::Quaternion(const ActionOptions&ao):
//PLUMED_COLVAR_INIT(ao)
//Initialise posre1 to hold 1 position and 0 derivatives.
//Initialise mypack to be empty with no arguments and no atoms.
//PLUMED_COLVAR_INIT(ao),posder1(1,0), mypack1(0,0,posder1), posder2(1,0)
PLUMED_COLVAR_INIT(ao),
nqdummy(4,0.0)
{
    fprintf (stderr, "= = Debug: Constructing the QUATERNION colvar...\n");
    bRefc1_is_calculated = false;
    bRef1_is_centered = false;
    bRefc2_is_calculated = false;
    bRef2_is_centered = false;
    
    string reffile1;
    string reffile2;
    parse("REFERENCE",reffile1);
    parse("REFERENCE_B", reffile2);
    string type;
    type.assign("OPTIMAL");
    string normdir;
    parse("NORM_DIRECTION", normdir);
    //parse("TYPE",type);

    //Setup the normalisation direction of q. Initialise a vector
    //and then give its address to the COLVAR.
    Vector4d nq(1.0,0.0,0.0,0.0);
    fprintf (stderr, "= = = Debug: normdir argument: %s, comparisons %i %i %i %i\n",
            normdir.c_str(),
            normdir.compare("w"),normdir.compare("x"),normdir.compare("y"),normdir.compare("z"));
    if ( normdir.compare("w")==0 ) {Vector4d dq(1.0, 0.0, 0.0, 0.0); nq=dq;}
    else if ( normdir.compare("x")==0 ) {Vector4d dq(0.0, 1.0, 0.0, 0.0); nq=dq;}
    else if ( normdir.compare("y")==0 ) {Vector4d dq(0.0, 0.0, 1.0, 0.0); nq=dq;}
    else if ( normdir.compare("z")==0 ) {Vector4d dq(0.0, 0.0, 0.0, 1.0); nq=dq;}
    else {
        fprintf (stderr, "= = = Debug: Detected custom q input...\n");
        std::vector<double> qvals(4);
        stringToQuat(normdir,qvals,",");
        //parseVector("NORM_DIRECTION",qvals);
        if ( qvals.size() != 4 ) {
            fprintf (stderr, "= = = Debug: WRONG qvals: %g %g %g %g\n",
                     qvals[0],qvals[1],qvals[2],qvals[3]);
            //error("NORM_DIRECTION does not contain a q ");
        }
        Vector4d dq(qvals[0],qvals[1],qvals[2],qvals[3]);
        dq /= dq.modulo(); //Normalise for safety & flexibility.
    }
    normquat = nq;
    fprintf (stderr, "= = = Debug: normalisation-q: %g %g %g %g\n", nq[0], nq[1], nq[2], nq[3]);    
    Vector4d qq; qq.zero();
    quat = qq;
    checkRead();
    
    fprintf (stderr, "= = = Debug: Will now add components...\n");
    addComponentWithDerivatives("w"); componentIsNotPeriodic("w");
    addComponentWithDerivatives("x"); componentIsNotPeriodic("x");
    addComponentWithDerivatives("y"); componentIsNotPeriodic("y");
    addComponentWithDerivatives("z"); componentIsNotPeriodic("z");
    log<<"NOTE: q is stored as four components (w x y z).\n";
  
    PDB pdb1; //PDB storage of reference coordinates 1
    PDB pdb2; //PDB storage of reference coordinates 2
  
    fprintf (stderr, "= = = Debug: Will now read REFERENCE pdb file...\n");  
    // read everything in ang and transform to nm if we are not in natural units
    if( !pdb1.read(reffile1,plumed.getAtoms().usingNaturalUnits(),0.1/atoms.getUnits().getLength()) )
        error("missing input file " + reffile1 );
  
    //Store reference positions with zeroed barycenters.
    setReference1( pdb1 );
    fprintf(stderr, "= = = = Debug: read and parse finished.\n"); 
    
    if ( !reffile2.empty() ) {
        fprintf (stderr, "= = = Debug: REFERENCE_B found, parsing second pdb file...\n");          
        if( !pdb1.read(reffile2,plumed.getAtoms().usingNaturalUnits(),0.1/atoms.getUnits().getLength()) )
        error("missing input file " + reffile2 );
  
        //Store reference positions with zeroed barycenters.
        setReference2( pdb2 );
        fprintf (stderr, "= = = = Debug: read and parse finished ...\n");  
    }
  
    fprintf (stderr, "= = = Debug: Requesting atoms for tracking...\n");
    std::vector<AtomNumber> atoms;
    rmsd1 = metricRegister().create<RMSDBase>(type,pdb1);
    rmsd1->getAtomRequests( atoms );
    requestAtoms( atoms );
    fprintf (stderr, "= = = = Debug: Request finished.\n");    
    // Setup the derivative pack
    //posder1.resize( 1, 3*atoms.size()+9 ); mypack1.resize( 0, atoms.size() );
    //for(unsigned i=0;i<atoms.size();++i) mypack1.setAtomIndex( i, i );
  
    log.printf("  reference from file %s\n",reffile1.c_str());
    log.printf("  which contains %d atoms\n",getNumberOfAtoms());
   
  /*if ( !reference2.empty() ) {
    if ( !pdb2.read(reference2,plumed.getAtoms().usingNaturalUnits(),0.1/atoms.getUnits().getLength()) )
      error("missing input file " + reference2 );      
    
    rmsd2 = metricRegister().create<RMSDBase>(type,pdb2);
    
    rmsd2->getAtomRequests( atoms );
    //  rmsd->setNumberOfAtoms( atoms.size() );
    requestAtoms( atoms );
    // Setup the derivative pack
    posder2.resize( 1, 3*atoms.size()+9 ); mypack2.resize( 0, atoms.size() );
    for(unsigned i=0;i<atoms.size();++i) mypack2.setAtomIndex( i, i );
  
    log.printf("  reference from file %s\n",reference2.c_str());
    log.printf("  which contains %d atoms\n",getNumberOfAtoms());
  }*/
    fprintf (stderr, "= = = Debug: finished constructing QUATERNION.\n");
}

Quaternion::~Quaternion(){
  delete rmsd1;
  if (rmsd2) delete rmsd2;
}

void Quaternion::clear(){
    refpos1.clear();
    refpos2.clear();
}

void Quaternion::setReference1(const std::vector<Vector> & inppos, const std::vector<double> & inpw, const std::vector<double> & inpdisp ){
  unsigned nvals=inppos.size();
  refpos1=inppos;
  plumed_massert(refw1.empty(),"you should first clear() an RMSD object, then set a new reference");
  plumed_massert(refdisp1.empty(),"you should first clear() an RMSD object, then set a new reference");
  refw1.resize(nvals,0.0);
  refdisp1.resize(nvals,0.0);
  double wtot=0, disptot=0;
  for(unsigned i=0;i<nvals;i++) {refc1+=inpw[i]*refpos1[i]; wtot+=inpw[i]; disptot+=inpdisp[i];}
  fprintf(stderr, "= = = = = Debug setReference1(): wtot is %g and disptot is %g\n", wtot, disptot);
  refc1/=wtot;
  for(unsigned i=0;i<nvals;i++) {refpos1[i]-=refc1; refw1[i]=inpw[i]/wtot ; refdisp1[i]=inpdisp[i]/disptot; }
  bRefc1_is_calculated=true;
  bRef1_is_centered=true;
}
void Quaternion::setReference1( const PDB inppdb ){
    setReference1( inppdb.getPositions(), inppdb.getOccupancy(), inppdb.getBeta() );
}

void Quaternion::setReference2(const std::vector<Vector> & inppos, const std::vector<double> & inpw, const std::vector<double> & inpdisp ){
  unsigned nvals=inppos.size();
  refpos2=inppos;
  plumed_massert(refw2.empty(),"you should first clear() an RMSD object, then set a new reference");
  plumed_massert(refdisp2.empty(),"you should first clear() an RMSD object, then set a new reference");
  refw2.resize(nvals,0.0);
  refdisp2.resize(nvals,0.0);
  double wtot=0, disptot=0;
  for(unsigned i=0;i<nvals;i++) {refc2+=inpw[i]*refpos2[i]; wtot+=inpw[i]; disptot+=inpdisp[i];}
  refc2/=wtot;
  for(unsigned i=0;i<nvals;i++) {refpos2[i]-=refc2; refw2[i]=inpw[i]/wtot ; refdisp2[i]=inpdisp[i]/disptot; }
  bRefc2_is_calculated=true;
  bRef2_is_centered=true;
}
void Quaternion::setReference2(const PDB inppdb ){
    setReference1( inppdb.getPositions(), inppdb.getOccupancy(), inppdb.getBeta() );
}

// Quaternion version of Optimal rotations. Copied from RMSD section to simplify
// class issues that I haven't mastered yet.
// Should return the quaternion
// This first case deals with only one absolute quaternion calculation.
void Quaternion::optimalAlignment1(const std::vector<Vector> & currpos)
{
    /* General procedure to calculate the quaternion value is to provide first a best fit.
    * Reference notation in NAMD, PLUMED, and Dill.
    * NOTES: Copy over RMSDCoreData::doCoreCalc first, and then do relevant pruning.
    * 1) Make the barycenters zero.
    * 2) gather the second moment/correlation matrix C = sum_1^N X1 X2T - is curly R in Ken Dill (2004), eq. 5
    * 3) Create the equivalent quaternion matrix/overlap matrix S  - is curly F, eq. 10
    * 4) Diagonalise this 4x4 matrix, the maximum eigenvector is almost always
    *    the best rotation, the first eigenvector by convention.
    * 5) Can use for both RMSD and rotation matrices.
    */
    //fprintf (stderr, "= = Debug: Starting optimalAlignment1...\n");
    
    const unsigned nvals = refpos1.size();
    //Make a copy of input vectors so as to center and modify them.
    std::vector<Vector> locpos;
    locpos = currpos;
    
    //double rr00=0.;
    //double rr11=0.;
    // This is positions*reference
    //Tensor rr01;

    // (1)
    //if (!bRef1_is_centered) { refc1 }
    Vector posc ; posc.zero();
    for(unsigned i=0;i<nvals;i++) {posc+=refw1[i]*locpos[i];}
    for(unsigned i=0;i<nvals;i++) {locpos[i]-=posc;}
    //fprintf (stderr, "= = = Debug: Current frame has been centered: %g %g %g\n",
    //                 posc[0], posc[1], posc[2]);
    // (2)
    // second expensive loop: compute second moments wrt centers
    // Already subtracted from the data.
    //Vector cp; cp.zero(); //if(!cpositions_is_removed)cp=cpositions;
    //Vector cr; cr.zero(); //if(!creference_is_removed)cr=creference;
    rr00=0.0; rr11=0.0; rr01.zero();
    for(unsigned iat=0;iat<nvals;iat++){
      double w=refw1[iat];
      rr00+=dotProduct(locpos[iat],locpos[iat])*w;
      rr11+=dotProduct(refpos1[iat],refpos1[iat])*w;
      rr01+=Tensor(locpos[iat],refpos1[iat])*w;
    }
    // (3) 
    // the quaternion matrix: this is internal
    Matrix<double> S=Matrix<double>(4,4);
    S[0][0]=2.0*(-rr01[0][0]-rr01[1][1]-rr01[2][2]);
    S[1][1]=2.0*(-rr01[0][0]+rr01[1][1]+rr01[2][2]);
    S[2][2]=2.0*(+rr01[0][0]-rr01[1][1]+rr01[2][2]);
    S[3][3]=2.0*(+rr01[0][0]+rr01[1][1]-rr01[2][2]);
    S[0][1]=2.0*(-rr01[1][2]+rr01[2][1]);
    S[0][2]=2.0*(+rr01[0][2]-rr01[2][0]);
    S[0][3]=2.0*(-rr01[0][1]+rr01[1][0]);
    S[1][2]=2.0*(-rr01[0][1]-rr01[1][0]);
    S[1][3]=2.0*(-rr01[0][2]-rr01[2][0]);
    S[2][3]=2.0*(-rr01[1][2]-rr01[2][1]);
    S[1][0] = S[0][1];
    S[2][0] = S[0][2];
    S[2][1] = S[1][2];
    S[3][0] = S[0][3];
    S[3][1] = S[1][3];
    S[3][2] = S[2][3];
    //fprintf (stderr, "= = = Debug: S overlap matrix calculated:\n");
    //for(unsigned i=0;i<4;i++)
    //    fprintf (stderr, "           [ %g %g %g %g ]\n", S[i][0], S[i][1], S[i][2], S[i][3]);
    
    // (4) Diagonalisation
    std::vector<double> eigenvals;
    Matrix<double> eigenvecs;

    //fprintf (stderr, "= = = Debug: Starting diagonalisation of matrix S...\n");
    int diagerror=diagMat(S, eigenvals, eigenvecs);
    if (diagerror!=0){
        std::string sdiagerror;
        Tools::convert(diagerror,sdiagerror);
        std::string msg="DIAGONALIZATION FAILED WITH ERROR CODE "+sdiagerror;
        plumed_merror(msg);
    } else {
        //fprintf(stderr, "= = = = Debug: diagMat() successful.\n");
    }
    std::vector<double> l = eigenvals; //Make copy for ease of reading.
    Matrix<double> q = eigenvecs; //Make copy for ease of reading.
    // Retrieve Quaternion and lambdas here, with a first check for normalisation
    // conditions
    double dot;
    //Normalise each eigenvector in the direction closer to norm
    for (unsigned i=0;i<4;i++) {
        dot=0.0;
        for (unsigned j=0;j<4;j++) {
            dot+=normquat[j]*q[i][j];
        }
        //printf(stderr,"Dot: %g q: %g %g %g %g\n",d dot, q[i][0], q[i][1], q[i][2], q[i][3]);
        if (dot < 0.0)
            for (unsigned j=0;j<4;j++)
                q[i][j]=-q[i][j];
    }
    for (unsigned i=0;i<4;i++) quat[i]=q[0][i];
    
    //setup ptr to components as a vector.
    Value* valuew=getPntrToComponent("w");
    Value* valuex=getPntrToComponent("x");
    Value* valuey=getPntrToComponent("y");
    Value* valuez=getPntrToComponent("z");
    std::vector<Value*> qptr;
    qptr.push_back(valuew); qptr.push_back(valuex);
    qptr.push_back(valuey); qptr.push_back(valuez);
    
    //Assign q-values to the COLVAR component.
    for (unsigned i=0;i<4;i++) qptr[i]->set(quat[i]);

    //fprintf(stderr, "= = = = Debug: Printing eigenvalues and eigenvectors:\n");
    //fprintf(stderr, "       lambda: [ %g %g %g %g ]\n", l[0], l[1], l[2], l[3]);
    //for(unsigned i=0;i<4;i++)
    //    fprintf(stderr, "       vec-%i: [ %g %g %g %g ]\n",
    //                i, eigenvecs[i][0],eigenvecs[i][1],eigenvecs[i][2], eigenvecs[i][3]);
    
    //Basic form, which we will reverse in the main for loop:
    //for (unsigned i=0;i<4;i++)
    //    for(unsigned j=0;j<getNumberOfAtoms();j++)
    //        setAtomsDerivatives(qptr[i], j, dqi/dxj );


    // (6) Finally calculate the atomic derivatives with respect to q
    // ... Since Sij = sum( x_k*xref_k ) if i=j, dx_j = xref_j ....
    
    // Matrix<double> dS=Matrix<double>(4,4,3);
    //Matrix<Vector> dSdx = Matrix<Vector>(4,4,3);
    std::vector<std::vector<Vector > > dSdx;
    dSdx.resize(4, std::vector<Vector>( 4 ));
    
    //Calculate derivatives for atom j.
    //This is ported in from NAMD code.
    //fprintf(stderr, "= = = Debug: Now calculating dSij/dx & dq/dx over %i atoms.\n", nvals);
    for (unsigned j=0;j<nvals;j++) {
        double const rx = refpos1[j][0];
        double const ry = refpos1[j][1];
        double const rz = refpos1[j][2];

        // dSijdx : derivative of the S matrix, w.r.t. atom x_j 
        dSdx[0][0] = Vector(  rx,  ry,  rz);
        dSdx[1][0] = Vector( 0.0,  rz, -ry);
        dSdx[0][1] = dSdx[1][0];
        dSdx[2][0] = Vector( -rz, 0.0,  rx);
        dSdx[0][2] = dSdx[2][0];
        dSdx[3][0] = Vector(  ry, -rx,  0.0);
        dSdx[0][3] = dSdx[3][0];
        dSdx[1][1] = Vector(  rx, -ry, -rz);
        dSdx[2][1] = Vector(  ry,  rx,  0.0);
        dSdx[1][2] = dSdx[2][1];
        dSdx[3][1] = Vector(  rz,  0.0,  rx);
        dSdx[1][3] = dSdx[3][1];
        dSdx[2][2] = Vector(- rx,  ry, -rz);
        dSdx[3][2] = Vector( 0.0,  rz,  ry);
        dSdx[2][3] = dSdx[3][2];
        dSdx[3][3] = Vector( -rx, -ry,  rz);

        // dqi/dxj = Sum_i Sum_j q1i dSijdx q0j /(Norm) * qi...
        Vector dqidxj;
        for (unsigned i=0;i<4;i++) {
            //Calculate and append the derivatives due to each q-component separately
            dqidxj.zero();
            for (unsigned a=0;a<4;a++) {
                for (unsigned b=0;b<4;b++) {
                    double fact=q[1][a]*q[0][b]/(l[0]-l[1])*q[1][i];
                    fact+=q[2][a]*q[0][b]/(l[0]-l[2])*q[2][i];
                    fact+=q[3][a]*q[0][b]/(l[0]-l[3])*q[3][i];
                    dqidxj+= -1.0*fact*dSdx[a][b];
                }
            }
            //
            setAtomsDerivatives(qptr[i], j, dqidxj);
        }
    }
    
    //Not that all derivatives have been calculated set the system box derivatives.
    for (unsigned i=0;i<4;i++)
        setBoxDerivativesNoPbc(qptr[i]);
    
    //fprintf (stderr, "= = = Debug: Finished optimalAlignment1. q = (%g %g %g %g)\n",
    //                 qq[0], qq[1], qq[2], qq[3]);
    return;
}

// calculator
void Quaternion::calculate(){
    //fprintf (stderr, "= = Debug: QUATERNION::calculate has been called.\n");    
    std::vector<Vector> q1(4);
    
    // Obtain current position, reference position, and 
    optimalAlignment1( getPositions() );
    //The values and derivations of the components will have been set
    //within the above function.

  //At the end, transfer atom derivatice from sub-function to the colvar.
  //for(unsigned i=0;i<getNumberOfAtoms();i++) setAtomsDerivatives( i, mypack1.getAtomDerivative(i) );
  //Transfer Virials.
  //Tensor virial; plumed_dbg_assert( !mypack.virialWasSet() );
  //setBoxDerivativesNoPbc();
}

}
}



