
#include <pdbTool.h>

#include <cdsSStream.h>
#include <cdsMap.h>
#include <simulation.h>

#include <cdsFstream.h>
#include <cdsIomanip.h>

#include <cstdio>
#include <cstddef>
#include <sys/stat.h>
//#include <unistd.h>



PDBTool::PDBTool(const char*    filename,
		 const AtomSel& sel     )
  : selection(sel), filename(filename), makeBackup(1), useXplorNames_(1) 
{
 if ( sel.size() == 0 ) {
   AtomSel tmp("all",Simulation::currentSimulation()); 
   selection.set( tmp );
 }
 
 aux1_.resize( sel.simulation()->numAtoms() ); aux1_.set(1.0);
 aux2_.resize( sel.simulation()->numAtoms() ); aux2_.set(1.0);
} /* constructor */

void
PDBTool::addRemark(const CDSString &newVal)
{
 myRemarksList_.append(newVal);
}

void
PDBTool::clearRemarks()
{
 myRemarksList_.resize(0);
}

CDSList< CDSString > 
PDBTool::remarks() const
{
 return myRemarksList_;
}

void
PDBTool::useXplorNames()
{
 useXplorNames_ = 1;
}

void
PDBTool::useIupacNames()
{
 useXplorNames_ = 0;
}

static CDSString
xplorAtomName(CDSString oldVal)
{
 //
 // Given an atom name, return a 4-char CDSString
 // in xplor name format, ie.,
 //
 // 3-char names have a leading space and numbers at end
 // 4-char names have numbers at end
 //

 CDSString newVal;
 CDSString aSpace(" ");
 CDSString numbers("0123456789");
  
 // strip out whitespace, limit length to 4 chars

 StringStream tempStream;
 tempStream << subString(oldVal, 0, 4);
 tempStream >> oldVal;
  
 if (numbers.contains(subString(oldVal, 0, 1)) && 
     subString(oldVal, 1, 2) == CDSString("H")) 
   oldVal = subString(oldVal, 1) + subString(oldVal, 0, 1);
 
 switch(oldVal.length()) {

   case 0: 
     newVal = aSpace + aSpace + aSpace + aSpace;
     break;
    
   case 1:
     newVal = aSpace + oldVal + aSpace + aSpace;
     break;

   case 2:
     newVal = aSpace + oldVal + aSpace;
     break;

   case 3:
     newVal = aSpace + oldVal;
     break;

   case 4:
     newVal = oldVal;
     break;
    
   default:
     newVal = subString(oldVal, 0, 4);
 }
  
 return newVal;
} /* xplorAtomName */

static CDSString
iupacAtomName(CDSString oldVal)
{
 //
 // Given an atom name, return a 4-char CDSString
 // in IUPAC name format, ie.,
 //
 // xplor names that that end in a number have that
 // number moved to the front
 //
 // names that don't then begin with a number get a leading space
 //

 CDSString newVal;
 CDSString aSpace(" ");
 CDSString numbers("0123456789");
 int trailingNumber;

 CDSString lastOldValChar = subString(oldVal, oldVal.length());

 if (numbers.contains(lastOldValChar))
   trailingNumber = 1;
 else
   trailingNumber = 0;
    
 if (trailingNumber) 
   newVal = subString(oldVal, oldVal.length()) + subString(oldVal, 0, oldVal.length() - 1);
 else 
   newVal = aSpace + oldVal;
  
 newVal = subString(newVal, 0, 4);
 return newVal;
} /* iupacAtomName */

static bool
fileExists(const CDSString& filename)
{
 struct stat statStruct;
 int ret = stat(filename,&statStruct);
 if ( ret==0 )
   return 1;
 else 
   return 0;
} /* fileExists */

static void
renameGreater(const CDSString& name,
		    int     cnt)
{
 OStringStream filename;
 filename << name << '_' << cnt << ends;

 if ( !fileExists( filename.str() ) )
   return;

 renameGreater(name,cnt+1);

 OStringStream gFilename;
 gFilename << name << '_' << (cnt+1) << ends;
 rename(filename.str(), gFilename.str() );
} /* renameGreater */

static void
renameBackupFiles(const CDSString& name)
{
 if ( !fileExists(name) ) 
   return;

 renameGreater(name,1);
		 
 OStringStream bFilename;
 bFilename << name << "_1" << ends;
 rename(name, bFilename.str());
} /* renameBackupFiles */

// 
// write PDB files
//

void
PDBTool::write()
{
 if ( makeBackup() )
   renameBackupFiles( filename() );

 ofstream outfile( filename() );

 if (!outfile) 
   throw CDS::exception(CDSString("unable to open output file ") 
			+ filename());

 for (int i = 0; i < myRemarksList_.size(); ++i) 
   outfile << "REMARK " << myRemarksList_[i] << endl;

 int nAtomLines = selection().size();

 for (int x = 0; x < nAtomLines; ++x) {

   const Atom curAtom = selection()[x];

   outfile << resetiosflags(ios::left) 
	   << setiosflags(ios::right);               //right justified
   outfile << "ATOM  ";
   outfile << setw(5) << x+1 << " ";

   CDSString curName(curAtom.atomName());
   if (useXplorNames_)
     curName = xplorAtomName(curName);
   else
     curName = iupacAtomName(curName);

   outfile << resetiosflags(ios::right) 
	   << setiosflags(ios::left);               //left justified
   if (curName.length() >= 4) 
     outfile << setw(4) << curName << " ";
   else
     outfile << " " << setw(3) << curName << " ";

   outfile << setw(4) << curAtom.residueName() << " ";

   outfile << resetiosflags(ios::left) 
	   << setiosflags(ios::right);               //right justified
   outfile << setw(4) << curAtom.residueNum() << "    ";

   outfile << setiosflags(ios::showpoint | ios::fixed)
	   << setprecision(3);

   outfile << setw(8) << curAtom.pos().x();
   outfile << setw(8) << curAtom.pos().y();
   outfile << setw(8) << curAtom.pos().z();

   outfile << setprecision(2)
	   << setiosflags(ios::showpoint | ios::fixed);

   outfile << setw(6) << aux1(curAtom)
	   << setw(6) << aux2(curAtom) << "      ";

   outfile << resetiosflags(ios::right) 
	   << setiosflags(ios::left);               //left justified
   outfile << setw(4) << curAtom.segmentName();

   outfile << endl;
 }

 outfile << "END " << endl;

 outfile.close();
} /* write */

static CDSString 
pdbHash(const CDSString &segName, 
	const int resNum, 
	const CDSString &atomName) 
{
 StringStream tempStream;
 tempStream << atomName;
 tempStream << "__";
 tempStream << resNum;
 tempStream << "__";
 tempStream << segName;
 tempStream << ends;

 return tempStream.str();
} /* pdbHash */

//
// read PDB files
//

void 
PDBTool::read() 
{

 // open the file

 ifstream infile( filename() );

 if (!infile) 
   throw CDS::exception(CDSString("unable to open input file ") + 
			filename());

 // clear the old remarks

 this->clearRemarks();

 // index the selection for fast lookups during line-by-line read

 CDSMap<CDSString,int> selIndex;

 for (int x = 0; x < selection().size(); x++) {
   const Atom curAtom = selection()[x];
   CDSString curKey = pdbHash(curAtom.segmentName(), 
			   curAtom.residueNum() , 
			   curAtom.atomName()   );
   selIndex[curKey] = curAtom.index();
 }


 Simulation* sim = selection.get().simulation();
 // for each line,
 while ( infile.good()) {

   CDSString curLine;
   readline(infile, curLine);

   // check the line type

   CDSString lineType = subString(curLine, 0, 6);    

   if (lineType == "REMARK") {

     // handle remarks

     CDSString *rem = new CDSString;
     *rem = subString(curLine, 6, -1);

     this->addRemark(*rem);

   } else if (lineType == "ATOM  ") {

     // handle atom lines

     int id, resNum;
     CDSString atomName, resName, chainName, insertion, insertion2, segName;
     float x, y, z, b, q;

     StringStream tempStream;

     tempStream.clear();
     tempStream << subString(curLine, 6, 11);
     tempStream >> id;

     tempStream.clear();
     tempStream << subString(curLine, 12, 16);
     tempStream >> atomName;

     tempStream.clear();
     tempStream << subString(curLine, 16, 17);
     tempStream >> insertion2;

     tempStream.clear();
     tempStream << subString(curLine, 17, 21);
     tempStream >> resName;

     tempStream.clear();
     tempStream << subString(curLine, 21, 22);
     tempStream >> chainName;

     tempStream.clear();
     tempStream << subString(curLine, 22, 26);
     tempStream >> resNum;

     tempStream.clear();
     tempStream << subString(curLine, 26, 27);
     tempStream >> insertion;

     tempStream.clear();
     tempStream << subString(curLine, 30, 38);
     tempStream >> x;
     
     tempStream.clear();
     tempStream << subString(curLine, 38, 46);
     tempStream >> y;

     tempStream.clear();
     tempStream << subString(curLine, 46, 54);
     tempStream >> z;

     tempStream.clear();
     tempStream << subString(curLine, 54, 60);
     tempStream >> q;

     tempStream.clear();
     tempStream << subString(curLine, 60, 66);
     tempStream >> b;

     tempStream.clear();
     tempStream << subString(curLine, 72, 76);
     tempStream >> segName;
 
     tempStream.clear();

     if (useXplorNames_) 
       tempStream << xplorAtomName(atomName);
     else
       tempStream << iupacAtomName(atomName);

     tempStream >> atomName;

     //
     // Have all the info from this line.
     // See if this atom is in the selection.
     // If it is, reset its position.
     //

     CDSString curKey = pdbHash(segName, resNum, atomName);
      
     if (selIndex.exists(curKey)) {
       int index = selIndex[curKey];
       sim->setAtomPos(index, CDSVec3(x,y,z));
       aux1_[index] = q;
       aux2_[index] = b;
     }
   }
 }

 sim->sync();
 infile.close();
} /* read */

float_type
PDBTool::aux1(const Atom& atom)
{
 return aux1_[atom.index()];
} /* aux1 */

float_type
PDBTool::aux2(const Atom& atom)
{
 return aux2_[atom.index()];
} /* aux2 */

void
PDBTool::setAux1(const Atom&       atom,
		 const float_type& val )
{
 aux1_[atom.index()] = val;
} /* setAux1 */

void
PDBTool::setAux2(const Atom&       atom,
		 const float_type& val )
{
 aux2_[atom.index()] = val;
} /* setAux2 */
