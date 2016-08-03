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
    std::vector<Vector> refpos1;
    std::vector<double> refw1;
    std::vector<double> refdisp1;
    Vector refc1;
    bool bRefc1_is_calculated;
    bool bRef1_is_centered;
    Matrix<double> S01;
    Vector4d quat1;
    double lambda01;
    
    PLMD::RMSDBase* rmsd2;
    std::vector<Vector> refpos2;
    std::vector<double> refw2;
    std::vector<double> refdisp2;
    Vector refc2;
    bool bRefc2_is_calculated;
    bool bRef2_is_centered;
    Matrix<double> S12;
    Vector4d quat2;
    double lambda12;
    
    // Relevant variables for outputs and derivatives.
    double dist;
    double lambda;
    Vector4d quat;
    Vector4d normquat; //Direction of normalisation.
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
  
  std::vector<Vector> rotateCoordinates(const Vector4d q, const std::vector<Vector> pos);
  
  void alignDomain(const unsigned nvals,
        const std::vector<Vector> ref, const std::vector<double> w,
        const std::vector<Vector> loc, const unsigned offset,
        Matrix<double> S, std::vector<double> eigenvals, Matrix<double> eigenvecs);
  
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
  
    // The colvar module takes in only 1 set of atoms.
    // Therefore, will need additional code to handle relative orientation from two diomains.
    fprintf (stderr, "= = = Debug: Requesting atoms for tracking...\n");
    if (!bRefc2_is_calculated) {
        // Simple case of 1 reference, absolute quaternion.
        std::vector<AtomNumber> atoms;
        
        rmsd1 = metricRegister().create<RMSDBase>(type,pdb1);
        rmsd1->getAtomRequests( atoms );
        fprintf(stderr,"= = = = Debug: First 3 entries of pdb1 atoms: %i %i %i\n",
                atoms[0], atoms[1], atoms[2]);
        
        requestAtoms( atoms );
        log.printf("  reference from file %s\n",reffile1.c_str());
        log.printf("  which contains %d atoms\n",getNumberOfAtoms());
    } else {
    //if (bRefc2_is_calculated) {        
        // 2 Reference files, relative quaternions.
        std::vector<AtomNumber> atoms;
        std::vector<AtomNumber> temp1;
        std::vector<AtomNumber> temp2;

        rmsd1 = metricRegister().create<RMSDBase>(type,pdb1);
        rmsd1->getAtomRequests( temp1 );
        rmsd2 = metricRegister().create<RMSDBase>(type,pdb2);
        rmsd2->getAtomRequests( temp2 );
        
        fprintf(stderr,"= = = = Debug: First 3 entries of pdb1 atoms: %i %i %i\n",
                temp1[0], temp1[1], temp1[2]);
        fprintf(stderr,"= = = = Debug: First 3 entries of pdb2 atoms: %i %i %i\n",
                temp2[0], temp2[1], temp2[2]);               
        
        atoms.reserve( temp1.size()+temp2.size() );
        atoms.insert( atoms.end(), temp1.begin(), temp1.end());
        atoms.insert( atoms.end(), temp2.begin(), temp2.end());

        requestAtoms( atoms );
        log.printf("  reference from files %s and %s\n",
                reffile1.c_str(), reffile2.c_str());
        log.printf("  which contains %d atoms\n",getNumberOfAtoms());    
    }
    fprintf (stderr, "= = = = Debug: Request finished.\n");    
    // Setup the derivative pack
    //posder1.resize( 1, 3*atoms.size()+9 ); mypack1.resize( 0, atoms.size() );
    //for(unsigned i=0;i<atoms.size();++i) mypack1.setAtomIndex( i, i );
  
   
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

//Returns a copy of coordinates pos, rotated according to q.
std::vector<Vector> Quaternion::rotateCoordinates(const Vector4d q, const std::vector<Vector> pos)
{
    std::vector<Vector> rot;
    unsigned ntot = pos.size();
    rot.resize(ntot,Vector(0,0,0));
    for (unsigned i=0;i<ntot;i++)
        rot[i]=quaternionRotate(q, pos[i]);
    return rot;
}

// Generic quaternion calculator. Able to function on submatrixes 
// Using formalism q being the rotation from local to reference frame.
void Quaternion::alignDomain(const unsigned nvals,
        const std::vector<Vector> ref, const std::vector<double> w,
        const std::vector<Vector> loc, const unsigned offset,
        Matrix<double> S, std::vector<double> l, Matrix<double> q)
{
    Vector posc ; posc.zero();
    unsigned ntot=nvals+offset;
    if (ntot>loc.size())
    {
        std::string msg="ALIGNDOMAIN FAILED: OFFSET + NVALS EXCEEDS CURRENT POSITION SIZE";
        plumed_merror(msg);
    }

    for(unsigned i=offset;i<ntot+;i++) {posc+=w[i]*loc[i];}
    for(unsigned i=offset;i<ntot;i++) {loc[i]-=posc;}
    //fprintf (stderr, "= = = Debug: Current frame has been centered: %g %g %g\n",
    //                 posc[0], posc[1], posc[2]);
    // (2)
    // second expensive loop: compute second moments wrt centers
    // Already subtracted from the data.
    //Vector cp; cp.zero(); //if(!cpositions_is_removed)cp=cpositions;
    //Vector cr; cr.zero(); //if(!creference_is_removed)cr=creference;
    rr00=0.0; rr11=0.0; rr01.zero();
    unsigned iloc;
    for(unsigned iat=0;iat<nvals;iat++){
        iloc=iat+offset;
        rr00+=dotProduct(loc[iloc],loc[iloc])*w[iat];
        rr11+=dotProduct(ref[iat],ref[iat])*w[iat];
        rr01+=Tensor(loc[iloc],ref[iat])*w[iat];
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
    //fprintf (stderr, "= = = Debug: Starting diagonalisation of matrix S...\n");
    
    // l is the eigenvalues, q are the eigenvectors/quaternions.
    int diagerror=diagMat(S, l, q);
    if (diagerror!=0){
        std::string sdiagerror;
        Tools::convert(diagerror,sdiagerror);
        std::string msg="DIAGONALIZATION FAILED WITH ERROR CODE "+sdiagerror;
        plumed_merror(msg);
    } else {
        //fprintf(stderr, "= = = = Debug: diagMat() successful.\n");
    }
    
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
}

// Quaternion version of Optimal rotations. Copied from RMSD section to simplify
// class issues that I haven't mastered yet.
// Should return the quaternion
// The first case deals with only one absolute quaternion calculation.
// The second case deals with a relative quaternion calculation, between two groups.
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
    * 6) Finally calculate the atomic derivatives with respect to q
    *      Sij = sum( x_k*xref_k ), so if i=j, dSdx_j = xref_j ....
    */
    //fprintf (stderr, "= = Debug: Starting optimalAlignment1...\n");

    /* When a second reference is given:
     * 0) Assume that two reference frames are in the same orientation (1,0,0,0).
     * 1) rotate the second reference PDB according to q1,
     * 2) calculate the quaternion value to this new frame. This is q1->2.
     * 3) Then apply forces.
     */
    
    /* Pseudocode for this function */
    /*
    alignDomain1();
    if (!bRefc2_is_calculated) {
        setSQLto1();
    } else {
        rotateReference2();
        alignDomain2();
        set SQLto12();
    }
    reportQuaterniontoOutput();
    calculateDerivatives(); That is...
    for (unsigned i=0;i<4;i++)
        for(unsigned j=0;j<getNumberOfAtoms();j++)
            setAtomsDerivatives(qptr[i], j, dqi/dxj );
    */

    Matrix<double> S;   
    Matrix<double> qmat;
    std::vector<double> l;    
    
    std::vector<Vector> rotref1;
    std::vector<Vector> rotref2;
    
    
    //Align the first domain,
    //cheat with fact that currpos is concatenated from domain 1 and domain 2.
    alignDomain(refpos1.size(), refpos1, refw1, currpos, 0, S01, l, qmat);
    lambda01 = l[0];
    for (unsigned i=0;i<4;i++) quat1[i]=qmat[0][i];
    
    if (!bRefc2_is_calculated) {
        S=S01;
        quat=quat1;
    } else {
        //Rotate reference to match q1, then calculate q_12.
        rotref2 = rotateCoordinates(quat1, refpos2);
        
        //cheat with fact that currpos is concatenated of domain 1 and domain 2.        
        alignDomain(refpos2.size(), rotref2, refw2, currpos, refpos1.size(), S12, l, qmat);
        
        lambda12 = l[0];
        // From q_12 we can obtain directly the rotation from reference q02,
        // and q_21.
        S=S12;
        for (unsigned i=0;i<4;i++) quat[i]= qmat[0][i];
        quat2 = quat*quaternionInvert(quat1);
    }
    //fprintf(stderr, "= = = = Debug: Printing eigenvalues and eigenvectors:\n");
    //fprintf(stderr, "       lambda: [ %g %g %g %g ]\n", lambda1[0], lambda1[1], lambda1[2], lambda1[3]);
    //for(unsigned i=0;i<4;i++)
    //    fprintf(stderr, "       vec-%i: [ %g %g %g %g ]\n",
    //                i, qmat1[i][0], qmat1[i][1], qmat1[i][2], qmat1[i][3]);
    
    //setup ptr to components as a vector, in order to manipulate outputs as a vector.
    Value* valuew=getPntrToComponent("w");
    Value* valuex=getPntrToComponent("x");
    Value* valuey=getPntrToComponent("y");
    Value* valuez=getPntrToComponent("z");
    std::vector<Value*> qptr;
    qptr.push_back(valuew); qptr.push_back(valuex);
    qptr.push_back(valuey); qptr.push_back(valuez);
    //Assign q-values to the COLVAR component.
    for (unsigned i=0;i<4;i++) qptr[i]->set(quat[i]);    
    
    // Matrix<double> dS=Matrix<double>(4,4,3);
    //Matrix<Vector> dSdx = Matrix<Vector>(4,4,3);
    std::vector<std::vector<Vector > > dSdx;
    dSdx.resize(4, std::vector<Vector>( 4 ));
    
    //Calculate derivatives for atom j
    //This is ported in from NAMD code.
    //fprintf(stderr, "= = = Debug: Now calculating dSij/dx & dq/dx over %i atoms.\n", nvals);

    if (!bRefc2_is_calculated) {
        //Calculate for group 1 versus reference.
        for (unsigned j=0;j<refpos1.size();j++) {
            if (refw1[j]==0) continue; //Only apply forces to weighted atoms in the RMSD calculation.
            
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
                        double fact=qmat[1][a]*qmat[0][b]/(l[0]-l[1])*qmat[1][i];
                              fact+=qmat[2][a]*qmat[0][b]/(l[0]-l[2])*qmat[2][i];
                              fact+=qmat[3][a]*qmat[0][b]/(l[0]-l[3])*qmat[3][i];
                        dqidxj+= -1.0*fact*dSdx[a][b];
                    }
                }
                //
                setAtomsDerivatives(qptr[i], j, dqidxj);
            }
        }
    } else {
        //We have the matrix for q_12. Need to do the reverse, and using the rotated reference atoms as the base.
        //Cheat with the concatenation again.
        //Do group 1 first, then group 2.
        for (unsigned j=0;j<rotref1.size();j++) {
            if (refw1[j]==0) continue; //Only apply forces to weighted atoms in the RMSD calculation.
            
            double const rx = rotref1[j][0];
            double const ry = rotref1[j][1];
            double const rz = rotref1[j][2];

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
                        double fact=qmat[1][a]*qmat[0][b]/(l[0]-l[1])*qmat[1][i];
                              fact+=qmat[2][a]*qmat[0][b]/(l[0]-l[2])*qmat[2][i];
                              fact+=qmat[3][a]*qmat[0][b]/(l[0]-l[3])*qmat[3][i];
                        dqidxj+= 1.0*fact*dSdx[a][b];
                    }
                }
                setAtomsDerivatives(qptr[i], j, dqidxj);
            }
        }
        
        unsigned offset=rotref1.size();
        for (unsigned j=0;j<rotref2.size();j++) {
            if (refw1[j]==0) continue; //Only apply forces to weighted atoms.
            
            double const rx = rotref2[j][0];
            double const ry = rotref2[j][1];
            double const rz = rotref2[j][2];

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
                        double fact=qmat[1][a]*qmat[0][b]/(l[0]-l[1])*qmat[1][i];
                              fact+=qmat[2][a]*qmat[0][b]/(l[0]-l[2])*qmat[2][i];
                              fact+=qmat[3][a]*qmat[0][b]/(l[0]-l[3])*qmat[3][i];
                        dqidxj+= -1.0*fact*dSdx[a][b];
                    }
                }
                setAtomsDerivatives(qptr[i], j+offset, dqidxj);
            }
        }
    }
    //Now that all derivatives have been calculated set the system box derivatives.
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



