#ifndef FF_MOLECULE_H
#define FF_MOLECULE_H

#include "BasicTypes.h"          //For uint
#include "EnsemblePreprocessor.h"
#include "SubdividedArray.h"
#include "Geometry.h"            //members
#include "CBMC.h"

#include <string>
#include <vector>

#include <cassert>


namespace mol_setup
{ 
   class Atom;
   class MolKind; 
}
namespace ff_setup
{
   class Bond;
   class FFBase;
}
namespace cbmc { class TrialMol; }

class FFSetup;
class PRNG;
struct MolPick;
class System;
class Forcefield;
class Setup;

class MoleculeKind
{
 public:
   
    MoleculeKind();
   ~MoleculeKind();
   
   uint NumAtoms() const { return numAtoms; }
   uint NumBonds() const { return bondList.count; }
   uint NumAngles() const { return angles.Count(); }
   uint NumDihs() const { return dihedrals.Count(); }
   uint AtomKind(const uint a) const { return atomKind[a]; }
   double AtomCharge(const uint a) const { return atomCharge[a]; }
   
   //Initialize this kind
   //Exits program if param and psf files don't match
   void Init(std::string const& l_name,
             Setup const& setup,
             Forcefield const& forcefield,
	     System & sys);
   
   //Invoke CBMC, oldMol and newMol will be modified
   void Build(cbmc::TrialMol& oldMol, cbmc::TrialMol& newMol,
	      const uint molIndex)
   { builder->Build(oldMol, newMol, molIndex); }

   //Invoke CBMC, oldMol and newMol will be modified
   void BuildIDNew(cbmc::TrialMol& newMol, const uint molIndex)
   { builder->BuildIDNew(newMol, molIndex); }

   void BuildIDOld(cbmc::TrialMol& oldMol, const uint molIndex)
   { builder->BuildIDOld(oldMol, molIndex); }

   void BuildNew(cbmc::TrialMol& newMol, const uint molIndex)
   { builder->BuildNew2(newMol, molIndex); }

   void BuildOld(cbmc::TrialMol& oldMol, const uint molIndex)
   { builder->BuildOld2(oldMol, molIndex); }

   void BuildGrowNew(cbmc::TrialMol& newMol, const uint molIndex)
   { builder->BuildGrowNew(newMol, molIndex); }

   void BuildGrowOld(cbmc::TrialMol& oldMol, const uint molIndex)
   { builder->BuildGrowOld(oldMol, molIndex); }

   void Regrowth(cbmc::TrialMol& oldMol, cbmc::TrialMol& newMol,
	      const uint molIndex)
   { builder->Regrowth(oldMol, newMol, molIndex); }

   
   double PrintChargeInfo();
   SortedNonbond sortedNB, sortedNB_1_4, sortedNB_1_3, sortedEwaldNB;


   //these are used for total energy calculations, see Geometry.h/cpp
   Nonbond nonBonded;
   Nonbond_1_4 nonBonded_1_4;
   Nonbond_1_3 nonBonded_1_3;
   EwaldNonbond nonEwaldBonded;
   
   BondList bondList;
   GeomFeature angles;
   GeomFeature dihedrals;
   bool oneThree, oneFour;

   std::string name;
   std::vector<std::string> atomNames;
   double molMass;
   
   double * atomMass;
   XYZ cavDim;
   uint exchangeRatio;
   
#if ENSEMBLE == GCMC
   double chemPot;
#endif

 private:
   
   void InitAtoms(mol_setup::MolKind const& molData);
    
   //uses buildBonds to check if molecule is branched
   //bool CheckBranches();
   void InitCBMC(System& sys, Forcefield& ff,
		 Setup& set);
   
   cbmc::CBMC* builder;
   
   uint numAtoms;
   uint * atomKind;
   double * atomCharge;
};





#endif /*FF_MOLECULE_H*/
