#include "NoEwald.h"
#include "EwaldCached.h"
#include "CalculateEnergy.h"        
#include "EnergyTypes.h"            //Energy structs
#include "EnsemblePreprocessor.h"   //Flags
#include "BasicTypes.h"             //uint
#include "System.h"                 //For init
#include "StaticVals.h"             //For init
#include "Forcefield.h"             //
#include "MoleculeLookup.h"
#include "MoleculeKind.h"
#include "Coordinates.h"
#include "BoxDimensions.h"
#include "TrialMol.h"
#include "GeomLib.h"
#include "NumLib.h"
#include <cassert>
#include <omp.h>

//
//   
//    Energy Calculation functions for Ewald summation method
//    Calculating self, correction and reciprocate part of ewald    
//
//    Developed by Y. Li and Mohammad S. Barhaghi
// 
//

using namespace geom;

NoEwald::NoEwald(StaticVals const& stat, System & sys) :
  EwaldCached(stat, sys) {}

void NoEwald::Init()
{

   electrostatic = forcefield.electrostatic;
   ewald = forcefield.ewald; 
   SetNull();
      
}


void NoEwald::AllocMem()
{  
   return;     
}


void NoEwald::RecipInit(uint box, BoxDimensions const& boxAxes)
{
   return;   
}


//calculate reciprocate term for a box
void NoEwald::BoxReciprocalSetup(uint box, XYZArray const& molCoords)
{
   return;   
}


//calculate reciprocate term for a box
double NoEwald::BoxReciprocal(uint box) const
{
   return 0.0; 
}


//calculate reciprocate force term for a box
Virial NoEwald::ForceReciprocal(Virial& virial, uint box) const
{
   return virial; 
}

//calculate correction force term for a box
Virial NoEwald::ForceCorrection(Virial& virial, uint box) const
{
   return virial; 
}


//calculate reciprocate term for displacement and rotation move
double NoEwald::MolReciprocal(XYZArray const& molCoords,
			      const uint molIndex, const uint box,
			      XYZ const*const newCOM)
{
   return 0.0; 
}


//calculate self term for a box
double NoEwald::BoxSelf(BoxDimensions const& boxAxes, uint box) const
{
   return 0.0;
}


//calculate correction term for a molecule
double NoEwald::MolCorrection(uint molIndex, uint box) const
{
   return 0.0;
}


//calculate reciprocate term in destination box for swap move
double NoEwald::SwapDestRecip(const cbmc::TrialMol &newMol,
			      const uint box, const int molIndex) 
{
   return 0.0;
}

//calculate reciprocate term for inserting some molecules kind A in destination
// box and removing molecule kind B from dest box
double NoEwald::SwapRecip(const std::vector<cbmc::TrialMol> &newMol,
			  const std::vector<cbmc::TrialMol> &oldMol)
{
  return 0.0;
}


//calculate reciprocate term in source box for swap move
double NoEwald::SwapSourceRecip(const cbmc::TrialMol &oldMol,
				const uint box, const int molIndex) 
{
   return 0.0;
}


//calculate self term after swap move
double NoEwald::SwapSelf(const cbmc::TrialMol& trialMol) const
{
   return 0.0;
}

//calculate correction term after swap move
double NoEwald::SwapCorrection(const cbmc::TrialMol& trialMol) const
{
   return 0.0;
}

//back up reciptocate value to Ref (will be called during initialization)
void NoEwald::SetRecipRef(uint box)
{
   return;
}

//update reciprocate values
void NoEwald::UpdateRecip(uint box)
{
   return;
}

//update the hx,y,z hsqr and prefact
void NoEwald::UpdateRecipVec(uint box)
{
  return;
}

//restore cosMol and sinMol
void NoEwald::RestoreMol(int molIndex)
{
   return;
}

//restore the whole cosMolRef & sinMolRef into cosMolBoxRecip & sinMolBoxRecip
void NoEwald::exgMolCache()
{
   return;
}


