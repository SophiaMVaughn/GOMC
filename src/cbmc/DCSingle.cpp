#include "DCSingle.h" 
#include "TrialMol.h" 
#include "DCData.h" 
#include "PRNG.h" 
#include "CalculateEnergy.h" 
#include "XYZArray.h" 
#include "Forcefield.h" 
 
namespace cbmc 
{ 
 
   void DCSingle::BuildOld(TrialMol& oldMol, uint molIndex) 
   { 
      PRNG& prng = data->prng; 
      XYZArray& positions = data->positions; 
      uint nLJTrials = data->nLJTrialsFirst; 
      double* inter = data->inter; 
      double* real = data->real; 
      double stepWeight = 0; 
       
      std::fill_n(inter, nLJTrials, 0.0);  
      std::fill_n(real, nLJTrials, 0.0); 
 
      prng.FillWithRandom(data->positions, nLJTrials, 
			  data->axes.GetAxis(oldMol.GetBox())); 
      positions.Set(0, oldMol.AtomPosition(atom)); 
      data->calc.ParticleInter(inter, real, positions, atom, molIndex, 
                               oldMol.GetBox(), nLJTrials);  
 
      for (uint trial = 0; trial < nLJTrials; ++trial) 
      { 
         stepWeight += exp(-1 * data->ff.beta * 
			   (inter[trial] + real[trial])); 
      } 
      oldMol.MultWeight(stepWeight); 
      oldMol.AddEnergy(Energy(0.0, 0.0, inter[0], real[0], 
			      0.0, 0.0, 0.0)); 
      oldMol.ConfirmOldAtom(atom); 
   } 
 
   void DCSingle::BuildNew(TrialMol& newMol, uint molIndex) 
   { 
      PRNG& prng = data->prng; 
      XYZArray& positions = data->positions; 
      uint nLJTrials = data->nLJTrialsFirst; 
      double* inter = data->inter; 
      double* real = data->real; 
      double* ljWeights = data->ljWeights; 
 
      std::fill_n(inter, nLJTrials, 0.0);  
      std::fill_n(real, nLJTrials, 0.0); 
      std::fill_n(ljWeights, nLJTrials, 0.0); 
 
      prng.FillWithRandom(positions, nLJTrials, 
			  data->axes.GetAxis(newMol.GetBox())); 
      if(newMol.seedToGrow)
      {
	//set the first atom in position of the other picked molecule.
	positions.Set(0, newMol.sCoords);
      }
      data->calc.ParticleInter(inter, real, positions, atom, molIndex, 
                               newMol.GetBox(), nLJTrials);

      double stepWeight = 0; 
      for (uint trial = 0; trial < nLJTrials; ++trial) 
      { 
	ljWeights[trial] = exp(-1 * data->ff.beta * 
			       (inter[trial] + real[trial])); 
         stepWeight += ljWeights[trial]; 
      } 
      if(newMol.seedToGrow)
      {
	//In IDExchange move we pick the winner.
	uint winner = 0;
      }
      else
      {
	uint winner = prng.PickWeighted(ljWeights, nLJTrials, stepWeight); 
      }
      newMol.MultWeight(stepWeight); 
      newMol.AddEnergy(Energy(0.0, 0.0, inter[winner], real[winner], 
			      0.0, 0.0, 0.0)); 
      newMol.AddAtom(atom, positions[winner]); 
   } 
}            
