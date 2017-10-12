#ifndef EWALD_H
#define EWALD_H

#include "BasicTypes.h"
#include "EnergyTypes.h"
#include "EwaldCached.h"
#include <vector>
#include <stdio.h>
#include <cstring>
#include <cassert>
#include <omp.h>
//
//    Calculating Electrostatic calculation without caching Fourier terms.
//    Energy Calculation functions for Ewald summation method
//    Calculating self, correction and reciprocate part of ewald    
//
//    Developed by Y. Li and Mohammad S. Barhaghi
// 
//

class StaticVals;
class System;
class Forcefield;
class Molecules;
class MoleculeLookup;
class MoleculeKind;
class Coordinates;
class COM;
class XYZArray;
class BoxDimensions;
class CalculateEnergy;


class Ewald : public EwaldCached 
{
  //friend class CalculateEnergy;
   public:

   Ewald(StaticVals const& stat, System & sys);

   virtual void Init();

   virtual void AllocMem();

   //initiliazie term used for ewald calculation
   virtual void RecipInit(uint box, BoxDimensions const& boxAxes);
   
   //setup reciprocate term for a box
   virtual void BoxReciprocalSetup(uint box, XYZArray const& molCoords);

   //calculate reciprocate energy term for a box
   virtual double BoxReciprocal(uint box) const;

   //calculate reciprocate term for displacement and rotation move
   virtual double MolReciprocal(XYZArray const& molCoords, const uint molIndex,
				const uint box, XYZ const*const newCOM = NULL);	

   //calculate reciprocate term for inserting one molecule in destination box
   virtual double SwapDestRecip(const cbmc::TrialMol &newMol, const uint box, 
			        const int molIndex);

   //calculate reciprocate term for removing one molecule in source box
   virtual double SwapSourceRecip(const cbmc::TrialMol &oldMol,
				  const uint box, const int molIndex);

   //calculate reciprocate term for inserting some molecules kind A in
   //destination box and removing molecule kind B from dest box
   virtual double SwapRecip(const std::vector<cbmc::TrialMol> &newMol,
			    const std::vector<cbmc::TrialMol> &oldMol);

   //calculate reciprocate term for changing coordinate of multiple molecule
   virtual double IntraSwapRecip(const std::vector<cbmc::TrialMol> &newMolA,
				 const std::vector<cbmc::TrialMol> &oldMolA,
				 const std::vector<cbmc::TrialMol> &newMolB,
				 const std::vector<cbmc::TrialMol> &oldMolB);

   //restore cosMol and sinMol
   virtual void RestoreMol(int molIndex);

   //update sinMol and cosMol
   virtual void exgMolCache();

};



#endif /*EWALD_H*/
