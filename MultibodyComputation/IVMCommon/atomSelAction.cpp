
#include "atomSelAction.h"

#include <simulation.h>
#include <matrixTools.h>
#include <fixedSymMatrix.h>
#include "Mat33.h"


void 
AtomSelAction::SetVelAction::run(Simulation* sim,
				 int         index)
{ 
 sim->setAtomVel(index,vel);
} /* SetVelAction */

void
AtomSelAction::Translate::run(Simulation* sim,
			      int         index)
{
 Vec3 newPos = sim->atomPos(index) + trans;
 sim->setAtomPos(index,newPos);
} /* Translate */

void
AtomSelAction::Rotate::run(Simulation* sim,
			   int         index)
{
 Vec3 newPos = center + rot * (sim->atomPos(index) - center);
 sim->setAtomPos(index,newPos);
} /* Rotate */

void 
AtomSelAction::Fit::init(const AtomSel& sel)
  // algorithm from Kabsch, Acta Cryst. A34, 827 (1978).
{
 if ( !fitBy.size() )
   return;

 CDSVector<Vec3> fitted;
 if ( altCoords ) {
   fitted = *altCoords;
   if ( altCoords->size() != sel.simulation()->numAtoms() )
     throw CDS::exception("Fit::init: altCoords has incorrect dimension");
 } else
   fitted = sel.simulation()->atomPosArr();

 Vec3 qcm(0.);
 Vec3 qcmTo(0.);
 float_type totMass=0.;
 for (int i=0 ; i<fitBy.size() ; i++) {
   totMass += fitBy[i].mass();
   qcm   += fitBy[i].mass() * fitted[ fitBy[i].index() ];
   qcmTo += fitBy[i].mass() * fitTo[ fitBy[i].index() ];
 }
 qcm /= totMass;
 qcmTo /= totMass;

 Mat33 R(0.);
 for (int i=0 ; i<fitBy.size() ; i++)
   for (int a=0 ; a<3 ; a++)
     for (int b=0 ; b<3 ; b++)
       R(a,b) += fitBy[i].mass() * (fitted[fitBy[i].index()]-qcm)[b] * 
		 (fitTo[ fitBy[i].index() ]-qcmTo)[a];

 Mat33 S(0.);
 for (int i=0 ; i<fitBy.size() ; i++)
   for (int a=0 ; a<3 ; a++)
     for (int b=0 ; b<3 ; b++)
       S(a,b) += fitBy[i].mass() * 
		 (fitted[fitBy[i].index()]-qcm)[a] * 
		 (fitted[fitBy[i].index()]-qcm)[b];
 
 typedef FixedSymMatrix<double,3> SymMat33;
 using MatrixTools::transpose;
 using MatrixTools::EigenResults;
 using MatrixTools::eigen;
 // FIX: SymMat3 RTR(transpose(R) * R): want generic matrix stuff
 Mat33 RTR = transpose(R) * R;
 SymMat33 tmp;
 for (int i=0 ; i<3 ; i++)
   for (int j=i ; j<3 ; j++)
     tmp(i,j) = RTR(i,j);

 EigenResults<SymMat33::MatrixType> results = eigen(tmp);

 //resulting eigenpairs are sorted from smallest to largest

 FixedVector<Vec3,3> a;
 a[0] = results.eigenPairs[2].vector.vector();
 a[1] = results.eigenPairs[1].vector.vector();
 a[2] = cross(a[0],a[1]);

 FixedVector<Vec3,3> b;
 b[0] = unitVec( R*a[0] );
 b[1] = unitVec( R*a[1] );
 b[2] = cross(b[0],b[1]);

 if ( dot(b[2],R*a[2])<0. ) b[2] *= -1;

 rot.set(0.);
 for (int i=0 ; i<3 ; i++)
   for (int j=0 ; j<3 ; j++)
     for (int k=0 ; k<3 ; k++)
       rot(i,j) += b[k][i] * a[k][j];

 trans = qcmTo - rot * qcm;
} /* Fit::init */

void
AtomSelAction::Fit::run(Simulation* sim,
			int         index)
{
 Vec3 pos = sim->atomPos(index);

 if ( altCoords )
   pos = (*altCoords)[index];

 pos = rot * pos + trans;

 if ( altCoords )
   (*altCoords)[index] = pos;
 else
   sim->setAtomPos(index, pos);
} /* Fit::run */

void
AtomSelAction::RMSD::init(const AtomSel& sel) 
{
 sum = 0; 
 size=sel.size();
 count = CDSMap<int,int>();
 byResidue_ = CDSMap<int,float_type>();
} /* RMSD::init */

void
AtomSelAction::RMSD::run(Simulation* sim,
			 int         index)
{
 Atom atom(sim,index);

 int resid = atom.residueNum();

 if ( !count.exists(resid) ) {
   count[resid] = 0;
   byResidue_[resid] = 0.;
 }

 count[resid]++;
 
 float_type val = abs2( atom.pos() - compare[ atom.index() ] );
 byResidue_[resid] += val;
 sum += val;

} /* RMSD::run */

void
AtomSelAction::RMSD::finish() 
{
 if ( size>0 && sum>0. ) sum = sqrt(sum/size); 

 for (int i=0 ; i<count.keys().size() ; i++) {
   int resid = count.keys()[i];
   byResidue_[resid] = sqrt( byResidue_[resid] / count[resid] );
 }
 
} /* RMSD::finish */

#include <cdsSStream.h>
