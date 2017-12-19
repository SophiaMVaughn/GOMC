#ifndef MOLCULETRANSFER_H
#define MOLCULETRANSFER_H

#if ENSEMBLE==GCMC || ENSEMBLE==GEMC

#include "MoveBase.h"
#include "TrialMol.h"

//#define DEBUG_MOVES

class MoleculeTransfer : public MoveBase
{
 public:

   MoleculeTransfer(System &sys, StaticVals const& statV) : 
   ffRef(statV.forcefield), molLookRef(sys.molLookupRef), 
   MoveBase(sys, statV), perAdjust(statV.GetPerAdjust()) 
   {
     //checking the acceptance statistic for each kind
     trial.resize(BOX_TOTAL);
     accept.resize(BOX_TOTAL);
     for(uint b = 0; b < BOX_TOTAL; b++)
     {
       trial[b].resize(sys.molLookupRef.GetNumKind(), 0);
       accept[b].resize(sys.molLookupRef.GetNumKind(), 0);	    
     }
     //
   }

   virtual uint Prep(const double subDraw, const double movPerc);
   virtual uint Transform();
   virtual void CalcEn();
   virtual void Accept(const uint earlyReject, const uint step);

 private:
   
   void PrintAcceptance(const uint step);
   double GetCoeff() const;
   uint GetBoxPairAndMol(const double subDraw, const double movPerc);
   MolPick molPick;
   uint sourceBox, destBox;
   uint pStart, pLen;
   uint molIndex, kindIndex;
   uint perAdjust;
   //check acceptance
   vector< vector<uint> > trial, accept;

   double W_tc, W_recip;
   double correct_old, correct_new, self_old, self_new;
   cbmc::TrialMol oldMol, newMol;
   Intermolecular tcLose, tcGain, recipLose, recipGain;
   MoleculeLookup & molLookRef;
   Forcefield const& ffRef;
};

inline void MoleculeTransfer::PrintAcceptance(const uint step)
{
  uint t = 0;
  for(uint b = 0; b < BOX_TOTAL; b++)
  {
    subPick = mv::GetMoveSubIndex(mv::MOL_TRANSFER, b);
    t += moveSetRef.GetTrial(subPick);
  }

  if((t + 1) % (10 * perAdjust) == 0)
  {
    for(uint b = 0; b < BOXES_WITH_U_NB; b++)
    {
      for(uint k = 0; k < molLookRef.GetNumKind(); k++)
      {
	if(trial[b][k] != 0)
	{
	  printf("MolTransfer acceptance for kind %d in Box %d: %1.5f \n",
		 k, b, ((double)accept[b][k] / trial[b][k]) * 100);
	}
      }
    }
  }
}

inline uint MoleculeTransfer::GetBoxPairAndMol(const double subDraw,
					       const double movPerc)
{
   // Need to call a function to pick a molecule that is not fixed but cannot be
   // swap between boxes. (beta != 1, beta !=2)
   uint state = prng.PickMolAndBoxPair2(molIndex, kindIndex, sourceBox, destBox,
					subDraw, movPerc);   

   if (state != mv::fail_state::NO_MOL_OF_KIND_IN_BOX)
   {
      pStart = pLen = 0;
      molRef.GetRangeStartLength(pStart, pLen, molIndex);
      trial[destBox][kindIndex]++;
   }
   return state;
}

inline uint MoleculeTransfer::Prep(const double subDraw, const double movPerc)
{
   uint state = GetBoxPairAndMol(subDraw, movPerc);
   newMol = cbmc::TrialMol(molRef.kinds[kindIndex], boxDimRef, destBox);
   oldMol = cbmc::TrialMol(molRef.kinds[kindIndex], boxDimRef, sourceBox);
   oldMol.SetCoords(coordCurrRef, pStart);
   return state;
}


inline uint MoleculeTransfer::Transform()
{
   cellList.RemoveMol(molIndex, sourceBox, coordCurrRef);
   molRef.kinds[kindIndex].Build(oldMol, newMol, molIndex);
   return mv::fail_state::NO_FAIL;
}

inline void MoleculeTransfer::CalcEn()
{
   W_tc = 1.0;
   W_recip = 1.0;
   correct_old = 0.0;
   correct_new = 0.0;
   self_old = 0.0;
   self_new = 0.0;

   if (ffRef.useLRC)
   {
      tcLose = calcEnRef.MoleculeTailChange(sourceBox, kindIndex, false);
      tcGain = calcEnRef.MoleculeTailChange(destBox, kindIndex, true);
      W_tc = exp(-1.0*ffRef.beta*(tcGain.energy + tcLose.energy));
   }
   
   if (newMol.GetWeight() != 0.0)
   {
      correct_new = calcEwald->SwapCorrection(newMol);
      correct_old = calcEwald->SwapCorrection(oldMol);
      self_new = calcEwald->SwapSelf(newMol);
      self_old = calcEwald->SwapSelf(oldMol);
      recipGain.energy =
	calcEwald->SwapDestRecip(newMol, destBox, molIndex);
      recipLose.energy =
	calcEwald->SwapSourceRecip(oldMol, sourceBox, molIndex);
      //need to contribute the self and correction energy 
      W_recip = exp(-1.0 * ffRef.beta * (recipGain.energy + recipLose.energy +
					 correct_new - correct_old +
					 self_new - self_old));
   }

}

inline double MoleculeTransfer::GetCoeff() const
{
#if ENSEMBLE == GEMC
   return (double)(molLookRef.NumKindInBox(kindIndex, sourceBox)) /
      (double)(molLookRef.NumKindInBox(kindIndex, destBox) + 1) *
      boxDimRef.volume[destBox] * boxDimRef.volInv[sourceBox];
#elif ENSEMBLE == GCMC
   if (sourceBox == mv::BOX0) //Delete case
   {
     if(ffRef.isFugacity)
     {
       return (double)(molLookRef.NumKindInBox(kindIndex, sourceBox)) *
         boxDimRef.volInv[sourceBox] /
         (BETA * molRef.kinds[kindIndex].chemPot);
     }
     else
     {
       return (double)(molLookRef.NumKindInBox(kindIndex, sourceBox)) *
         boxDimRef.volInv[sourceBox] *
         exp(-BETA * molRef.kinds[kindIndex].chemPot);
     }
   }
   else //Insertion case
   {
     if(ffRef.isFugacity)
     {
       return boxDimRef.volume[destBox]/
         (double)(molLookRef.NumKindInBox(kindIndex, destBox)+1) *
         (BETA * molRef.kinds[kindIndex].chemPot);
     }
     else
     {
       return boxDimRef.volume[destBox]/
         (double)(molLookRef.NumKindInBox(kindIndex, destBox)+1) *
         exp(BETA * molRef.kinds[kindIndex].chemPot);
     }
   }
#endif
}

inline void MoleculeTransfer::Accept(const uint rejectState, const uint step)
{
   bool result;   
   //print acceptance information
   PrintAcceptance(step);

   //If we didn't skip the move calculation
   if(rejectState == mv::fail_state::NO_FAIL)
   {
      double molTransCoeff = GetCoeff();
      double Wo = oldMol.GetWeight();
      double Wn = newMol.GetWeight();
      double Wrat = Wn / Wo * W_tc * W_recip;

      if((newMol.GetEnergy().real < 1.0e15) &&
	 (oldMol.GetEnergy().real < 1.0e15))
      {
	result = prng() < molTransCoeff * Wrat;
      }
      else
	result = false;

      if(result)
      {
	 //update acceptance
	 accept[destBox][kindIndex]++;
         //Add tail corrections
         sysPotRef.boxEnergy[sourceBox].tc += tcLose.energy;
         sysPotRef.boxEnergy[destBox].tc += tcGain.energy;
         //Add rest of energy.
         sysPotRef.boxEnergy[sourceBox] -= oldMol.GetEnergy();
         sysPotRef.boxEnergy[destBox] += newMol.GetEnergy();
	 //Add Reciprocal energy
	 sysPotRef.boxEnergy[sourceBox].recip += recipLose.energy;
	 sysPotRef.boxEnergy[destBox].recip += recipGain.energy;
	 //Add correction energy
	 sysPotRef.boxEnergy[sourceBox].correction -= correct_old;
	 sysPotRef.boxEnergy[destBox].correction += correct_new;
	 //Add self energy
	 sysPotRef.boxEnergy[sourceBox].self -= self_old;
	 sysPotRef.boxEnergy[destBox].self += self_new;

	 //Set coordinates, new COM; shift index to new box's list
         newMol.GetCoords().CopyRange(coordCurrRef, 0, pStart, pLen);
         comCurrRef.SetNew(molIndex, destBox);
         molLookRef.ShiftMolBox(molIndex, sourceBox, destBox,
				kindIndex);
	 cellList.AddMol(molIndex, destBox, coordCurrRef);


	 //Zero out box energies to prevent small number 
	 //errors in double.
	 if (molLookRef.NumInBox(sourceBox) == 0)
	 {
	    sysPotRef.boxEnergy[sourceBox].Zero();
	    sysPotRef.boxVirial[sourceBox].Zero();
	 }
	 else if (molLookRef.NumInBox(sourceBox) == 1)
	 {
	    sysPotRef.boxEnergy[sourceBox].inter = 0;
	    sysPotRef.boxVirial[sourceBox].inter = 0;
	    sysPotRef.boxEnergy[sourceBox].real = 0;
	    sysPotRef.boxVirial[sourceBox].real = 0;
	 }

	 for (uint b = 0; b < BOX_TOTAL; b++)
	 {
	    calcEwald->UpdateRecip(b);
	 }

	 //Retotal
         sysPotRef.Total();
      }
      else
      {
	 cellList.AddMol(molIndex, sourceBox, coordCurrRef);
	 //when weight is 0, MolDestSwap() will not be executed, thus cos/sin
	 //molRef will not be changed. Also since no memcpy, doing restore
	 //results in memory overwrite
	 if (newMol.GetWeight() != 0.0)
	    calcEwald->RestoreMol(molIndex);
	 
      }

   }
   else  //we didn't even try because we knew it would fail
      result = false;

   subPick = mv::GetMoveSubIndex(mv::MOL_TRANSFER, sourceBox);
   moveSetRef.Update(result, subPick, step);
}

#endif

#endif
