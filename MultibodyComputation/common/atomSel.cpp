

#include "atomSel.h"
#include "atomSelAction.h"

#include <simulation.h>

#include <cdsList.h>
#include <cdsString.h>
#include <cdsSStream.h>
#include <cdsMath.h>

typedef CDSList<int> IList;

AtomSel::AtomSel(const char*       selStr,
		       Simulation* sim) :
  sim_( sim==0 ? Simulation::currentSimulation() : sim ),
  ok_(0),
  sel(selStr), 
  boolList_(0,0)
{
 if ( sel.length() )
   aList_ = sim_->select( selStr );
} /* constructor */

AtomSel::AtomSel(const CDSList<int>& indices,
		       Simulation* sim) :
  sim_( sim==0 ? Simulation::currentSimulation() : sim ),
  ok_(0),
  sel("[from list of atom indices]"), 
  boolList_(0,0)
{
 aList_ = indices;
} /* constructor */

void 
AtomSel::setIndices(const CDSList<int>& indices)
{
 sel = "[from list of atom indices]";
 aList_ = indices;
} /* setIndices */

//void 
//AtomSel::apply(AtomSelAction::Base& cb)
//{
// cb.init(*this);
//
// for (int i=0 ; i<indices().size() ; i++)
//   cb.run(sim_,aList_[i]);
//
// cb.finish();
//} /* apply */

void 
AtomSel::apply(AtomSelAction::Base& cb) const
{
 cb.init(*this);

 for (int i=0 ; i<indices().size() ; i++)
   cb.run(sim_,aList_[i]);

 cb.finish();
 sim_->sync();
} /* apply */

bool
AtomSel::containsIndex(const int index) const
{
 if ( !boolList_.size() ) {
   boolList_.resize( sim_->numAtoms() );
   for (int i=0 ; i<sim_->numAtoms() ; i++)
     boolList_[i] = 0;
   for (int i=0 ; i<size() ; i++)
     boolList_[ aList_[i] ] = 1;
 }
 
 return boolList_[index];
} /* containsIndex */



ostream& 
operator<<(      ostream& s, 
	   const AtomSel& x)
{
 // s << x.simulation()->name();
 // s << " ( ";
 for (int i=0 ; i<x.indices().size() ; i++) {
   if (i>0) s<< " , ";
   s << x.simulation()->atomString(x.indices()[i]);
 }
 // s << " )";
 return s;
} /* operator<<(AtomSel) */



//static const String&
//getWord(istream &s)
//{
// static String ret;
// ret="";
// char c;
// s >> c;
// while ( s && c!=' ' && c!='\n' && c!='\t' && c!=')' && c!='(' ) {
//   if ( c=='"') {  //quoted input
//     ret += "\\\"";
//     s >> c;
//     while ( s && c!='"' ) {
//	 ret += c;
//	 s >> c;
//     }
//     ret += "\\\"";
//   } else
//     ret += c;
//   s >> c;
// }
// if ( s )
//   s.putback(c);
// return ret;
//} /* getWord */

//static const String&
//getSpace(istream &s)
//{
// static String ret;
// ret="";
// char c;
// s >> c;
// while ( s && (c==' ' || c=='\n' || c=='\t' || c==')' || c=='(') ) {
//   ret += c;
//   s >> c;
// }
// if ( s )
//   s.putback(c);
// return ret;
//} /* getSpace */

istream& 
operator>>(istream& s, AtomSel& x)
{
 //gobble up characters until find the beginning '('
 char c=s.get();
 while ( (c!='(') && !s.eof() )
   c=s.get();

 // read all the text until find the final ')'
 StringStream sstr;
 int np=0;
 s.unsetf(ios::skipws);
 while ( s ) {
   char c;
   s >> c;
   if ( c=='(' ) np++;   //allow for nested parentheses
   if ( c==')' )
     if ( np ) np--;
     else break;
   sstr << c;
 }
 sstr << ends;

 s.setf(ios::skipws);
 if (s.eof()) {
   cerr << "eof encountered in selection\n";
   return s;
 }

 x.sel = sstr.str();
 x.aList_ = x.simulation()->select( x.sel );

// selList.resize( max(selList.numMols(),x.mol()+1) );
//
// int i = -1;
// for (int j=0 ; j<selList.list[x.mol()].size() ; j++)
//   if ( selList.list[x.mol()][j].str == selstr ) {
//     i = j;
//     break;
//   }
// if ( i == -1 ) {
//   x.selIndex_ = selList.list[x.mol()].size();
//   selList.list[x.mol()].resize( x.selIndex_+1 );
//   selList.list[x.mol()][x.selIndex_].str = selstr;
// } else
//   x.selIndex_ = i;

 x.ok_ = 1; //guess everything's ok...
 return s;
} /* operator>>(AtomSel) */

