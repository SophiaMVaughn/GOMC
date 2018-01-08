#define _USE_MATH_DEFINES 
#include <math.h> 
#include "DCFreeHedronSeed.h" 
#include "DCData.h" 
#include "TrialMol.h" 
#include "MolSetup.h" 
#include "Forcefield.h" 
#include "PRNG.h" 
#include "NumLib.h" 
 
namespace cbmc 
{ 
 
   DCFreeHedronSeed::DCFreeHedronSeed(DCData* data,
				      const mol_setup::MolKind& kind,
				      uint focus, uint prev) 
      : data(data), hed(data, kind, focus, prev) 
   { 
         using namespace mol_setup; 
         using namespace std; 
         vector<Bond> onFocus = AtomBonds(kind, hed.Focus()); 
         for(uint i = 0; i < onFocus.size(); ++i) { 
            if (onFocus[i].a1 == prev) { 
               anchorBond = data->ff.bonds.Length(onFocus[i].kind); 
	       anchorBondKind = onFocus[i].kind; 
               break; 
            } 
         } 
   } 
 
 
   void DCFreeHedronSeed::PrepareNew(TrialMol& newMol, uint molIndex) 
   { 
      hed.PrepareNew(newMol, molIndex); 
   } 
 
   void DCFreeHedronSeed::PrepareOld(TrialMol& oldMol, uint molIndex) 
   { 
      hed.PrepareOld(oldMol, molIndex); 
      double bondLengthOld = sqrt(oldMol.OldDistSq(hed.Focus(), hed.Prev())); 
      oldBondEnergy =  data->ff.bonds.Calc(anchorBondKind, bondLengthOld); 
   } 
 
 
   void DCFreeHedronSeed::BuildNew(TrialMol& newMol, uint molIndex) 
   { 
      PRNG& prng = data->prng; 
      const CalculateEnergy& calc = data->calc; 
      const EwaldCached *calcEwald = data->calcEwald; 
      const Forcefield& ff = data->ff; 
      uint nLJTrials = data->nLJTrialsNth; 
      double* ljWeights = data->ljWeights; 
      double* inter = data->inter; 
      double* real = data->real; 
 
      std::fill_n(inter, nLJTrials, 0.0); 
      std::fill_n(real, nLJTrials, 0.0); 
      std::fill_n(ljWeights, nLJTrials, 0.0); 
 
      //get info about existing geometry 
      newMol.ShiftBasis(hed.Focus()); 
      const XYZ center = newMol.AtomPosition(hed.Focus()); 
      XYZArray* positions = data->multiPositions; 
 
      for (uint i = 0; i < hed.NumBond(); ++i) 
      { 
         positions[i].Set(0, newMol.RawRectCoords(hed.BondLength(i), 
                                                  hed.Theta(i), hed.Phi(i))); 
      } 
      //add anchor atom 
      positions[hed.NumBond()].Set(0, newMol.RawRectCoords(anchorBond, 0, 0)); 
 
      //counting backward to preserve prototype 
      for (uint lj = nLJTrials; lj-- > 0;) 
      { 
         //convert chosen torsion to 3D positions 
         RotationMatrix spin = 
            RotationMatrix::UniformRandom(prng(), prng(), prng()); 
         for (uint b = 0; b < hed.NumBond() + 1; ++b) 
         { 
               //find positions 
               positions[b].Set(lj, spin.Apply(positions[b][0])); 
               positions[b].Add(lj, center); 
         } 
      } 
 
      for (uint b = 0; b < hed.NumBond() + 1; ++b) 
      { 
         data->axes.WrapPBC(positions[b], newMol.GetBox());  
      } 
 
 
      for (uint b = 0; b < hed.NumBond(); ++b) 
      { 
	 calc.ParticleInter(inter, real, positions[b], hed.Bonded(b),  
			    molIndex, newMol.GetBox(), nLJTrials); 
      } 
      calc.ParticleInter(inter, real, positions[hed.NumBond()], hed.Prev(), 
                         molIndex, newMol.GetBox(), nLJTrials); 
 
      double stepWeight = 0; 
      for (uint lj = 0; lj < nLJTrials; ++lj) 
      { 
	ljWeights[lj] = exp(-ff.beta * (inter[lj] + real[lj])); 
         stepWeight += ljWeights[lj]; 
      } 
      uint winner = prng.PickWeighted(ljWeights, nLJTrials, stepWeight); 
      for(uint b = 0; b < hed.NumBond(); ++b) 
      { 
         newMol.AddAtom(hed.Bonded(b), positions[b][winner]); 
      } 
      newMol.AddAtom(hed.Prev(), positions[hed.NumBond()][winner]); 
      newMol.AddEnergy(Energy(hed.GetEnergy(), hed.GetNonBondedEn(), 
			      inter[winner], real[winner], 
			      0.0, 0.0, 0.0)); 
      newMol.MultWeight(hed.GetWeight()); 
      newMol.MultWeight(stepWeight / nLJTrials); 
   } 
 
   void DCFreeHedronSeed::BuildOld(TrialMol& oldMol, uint molIndex) 
   { 
      PRNG& prng = data->prng; 
      const CalculateEnergy& calc = data->calc; 
      const EwaldCached * calcEwald = data->calcEwald; 
      const Forcefield& ff = data->ff; 
      uint nLJTrials = data->nLJTrialsNth; 
      double* ljWeights = data->ljWeights; 
      double* inter = data->inter; 
      double* real = data->real; 
 
      std::fill_n(inter, nLJTrials, 0.0); 
      std::fill_n(real, nLJTrials, 0.0);  
      std::fill_n(ljWeights, nLJTrials, 0.0); 
 
      //get info about existing geometry 
      oldMol.SetBasis(hed.Focus(), hed.Prev()); 
      //Calculate OldMol Bond Energy & 
      //Calculate phi weight for nTrials using actual theta of OldMol 
      hed.ConstrainedAnglesOld(data->nAngleTrials - 1, oldMol, molIndex); 
      const XYZ center = oldMol.AtomPosition(hed.Focus()); 
      XYZArray* positions = data->multiPositions; 
      double prevPhi[MAX_BONDS]; 
      for (uint i = 0; i < hed.NumBond(); ++i) 
      { 
         //get position and shift to origin 
         positions[i].Set(0, oldMol.AtomPosition(hed.Bonded(i))); 
         data->axes.UnwrapPBC(positions[i], 0, 1, oldMol.GetBox(), center); 
         positions[i].Add(0, -center); 
      } 
      //add anchor atom 
      positions[hed.NumBond()].Set(0, oldMol.AtomPosition(hed.Prev())); 
      data->axes.UnwrapPBC(positions[hed.NumBond()], 0, 1, 
			   oldMol.GetBox(), center); 
      positions[hed.NumBond()].Add(0, -center); 
 
      //counting backward to preserve prototype 
      for (uint lj = nLJTrials; lj-- > 1;) 
      { 
         //convert chosen torsion to 3D positions 
         RotationMatrix spin = 
            RotationMatrix::UniformRandom(prng(), prng(), prng()); 
         for (uint b = 0; b < hed.NumBond() + 1; ++b) 
         { 
            //find positions 
            positions[b].Set(lj, spin.Apply(positions[b][0])); 
            positions[b].Add(lj, center); 
         } 
      } 
 
      for (uint b = 0; b < hed.NumBond() + 1; ++b) 
      { 
         positions[b].Add(0, center); 
         data->axes.WrapPBC(positions[b], oldMol.GetBox()); 
      } 
 
      for (uint b = 0; b < hed.NumBond(); ++b) 
      { 
	calc.ParticleInter(inter, real, positions[b], hed.Bonded(b), 
                            molIndex, oldMol.GetBox(), nLJTrials); 
      } 
      double stepWeight = 0; 
      calc.ParticleInter(inter, real, positions[hed.NumBond()], hed.Prev(), 
                         molIndex, oldMol.GetBox(), nLJTrials);  
 
      for (uint lj = 0; lj < nLJTrials; ++lj) 
      { 
	stepWeight += exp(-ff.beta * (inter[lj] + real[lj])); 
      } 

      for(uint b = 0; b < hed.NumBond(); ++b) 
      { 
         oldMol.ConfirmOldAtom(hed.Bonded(b)); 
      } 

      oldMol.ConfirmOldAtom(hed.Prev()); 
      oldMol.AddEnergy(Energy(hed.GetEnergy() + oldBondEnergy + 
			      hed.GetOldBondEn(), hed.GetNonBondedEn(), 
			      inter[0], real[0], 0.0, 0.0, 0.0)); 
      oldMol.MultWeight(hed.GetWeight()); 
      oldMol.MultWeight(stepWeight / nLJTrials); 
   } 
 
 
}              
