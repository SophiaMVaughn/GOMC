#ifndef IDENTITYEXCHANGE_H
#define IDENTITYEXCHANGE_H

#if ENSEMBLE==GCMC || ENSEMBLE==GEMC

#include "MoveBase.h"
#include "cbmc/TrialMol.h"

using std::vector;

// Identity Exchange Move:
//
// insert the molecule A inside the cavity that has the size of the molecule B
// and vice versa. Only 1 trial must be perform for the first seed.
// Mohammad Soroush Barhaghi. July 2017

class IdentityExchange : public MoveBase
{
 public:

   IdentityExchange(System &sys, StaticVals const& statV) :
    ffRef(statV.forcefield), molLookRef(sys.molLookupRef),
      MoveBase(sys, statV), rmax(statV.mol.kinds[0].rmax / 2.0) 
      {
	enableID = ffRef.enableID;	
	volCav = pow(2.0 * rmax, 3);

	if(enableID)
	{
	  uint kindSnum = 0;
	  uint kindLnum = 0;
	  uint max = 0;
	  for(uint k = 0; k < molLookRef.GetNumKind(); k++)
	  {
	    if(molRef.kinds[k].exchangeRatio != 0)
	    {
	      if(molRef.kinds[k].exchangeRatio > max)
	      {
		if(max == 0)
		{
		  max = molRef.kinds[k].exchangeRatio;
		  kindS = k;
		  kindSnum = molRef.kinds[k].exchangeRatio;
		}
		else
		{
		  kindL = kindS;
		  kindLnum = kindSnum;
		  kindS = k;
		  kindSnum = molRef.kinds[k].exchangeRatio;
		}
	      }
	      else
	      {
		kindL = k;
		kindLnum = molRef.kinds[k].exchangeRatio;
	      }
	    }
	  }

	  if(kindLnum != 0 && kindSnum != 0)
	  { 
	    exchangeRate = kindSnum / kindLnum;
	  }
	  else
	  {
	    printf("Error: Identity Exchange move is valid for exchanging one large molecule with N small molecules. Exchange Ratio of one kind must set to 1.\n");
	    exit(EXIT_FAILURE);
	  }
	}
      }

   virtual uint Prep(const double subDraw, const double movPerc);
   virtual uint Transform();
   virtual void CalcEn();
   virtual void Accept(const uint earlyReject, const uint step);

 private:

   void ShiftMol(const uint n, const uint from, const uint to);
   void RecoverMol(const uint n, const uint from, const uint to);
   uint PickMolInCav();
   uint ReplaceMolecule();
   void CalcTc();
   double GetCoeff() const;
   //calculate factorial
   double Factorial(const uint n) const;
   //calculate ratio of factorial
   double Factorial(const uint n, const uint count) const;
   uint GetBoxPairAndMol(const double subDraw, const double movPerc);
   MolPick molPick;
   uint sourceBox, destBox;
   vector<uint> pStartA, pLenA, pStartB, pLenB;
   vector<uint> molIndexA, kindIndexA, molIndexB, kindIndexB;
   bool insertB, enableID;
   uint numInCavA, numInCavB, exchangeRate, kindS, kindL;

   double rmax, volCav;
   XYZ center;
   double W_tc, W_recip;
   double correct_oldA, correct_newA, self_oldA, self_newA;
   double correct_oldB, correct_newB, self_oldB, self_newB;
   double recipDest, recipSource;
   vector<cbmc::TrialMol> oldMolA, newMolA, oldMolB, newMolB;
   Intermolecular tcNew[BOX_TOTAL];
   vector< vector<uint> > molInCav;
   MoleculeLookup & molLookRef;
   Forcefield const& ffRef;
};

inline uint IdentityExchange::PickMolInCav()
{
   uint state = mv::fail_state::NO_FAIL;
   //pick a random location in dense phase
   XYZ axis = boxDimRef.GetAxis(sourceBox);
   XYZ temp(prng.randExc(axis.x), prng.randExc(axis.y), prng.randExc(axis.z));
   //Use to shift the new insterted molecule
   center = temp;

   //Find the molecule kind 0 in the cavity
   if(calcEnRef.FindMolInCavity(molInCav, center, rmax, sourceBox, kindS,
				exchangeRate))
   {
     molIndexA.clear();
     kindIndexA.clear();
     molIndexB.clear();
     kindIndexB.clear();
     //Find the exchangeRate number of molecules kind 0 in cavity
     numInCavA = exchangeRate;
     for(uint n = 0; n < numInCavA; n++)
     {
       //pick random exchangeRate number of kindS in cavity
       uint mId = molInCav[kindS][prng.randIntExc(molInCav[kindS].size())];
       while(std::find(molIndexA.begin(), molIndexA.end(), mId) !=
	     molIndexA.end())
       {
	 mId = molInCav[kindS][prng.randIntExc(molInCav[kindS].size())];
       }
       molIndexA.push_back(mId);
       kindIndexA.push_back(molRef.GetMolKind(molIndexA[n]));
     } 

     //pick a molecule from Large kind in destBox
     numInCavB = 1;
     state = prng.PickMol(kindS, kindIndexB, molIndexB, numInCavB,
			  destBox);
     
   }
   else
   {
     //reject the move
     state = mv::fail_state::NO_MOL_OF_KIND_IN_BOX;
     //printf("Empty Cavity\n");
   }

   return state;
}


inline uint IdentityExchange::ReplaceMolecule()
{
   uint state = mv::fail_state::NO_FAIL;
   molIndexA.clear();
   kindIndexA.clear();
   molIndexB.clear();
   kindIndexB.clear();
   numInCavA = 1;
   numInCavB = exchangeRate;
   //pick a random molecule of Large kind in dens box
   state = prng.PickMol(kindS, kindIndexA, molIndexA, numInCavA, sourceBox);

   if(state == mv::fail_state::NO_FAIL)
   {
     //Use to shift to the COM of new molecule
     center = comCurrRef.Get(molIndexA[0]);
     //find how many of KindS exist in this center
     calcEnRef.FindMolInCavity(molInCav, center, rmax, sourceBox, kindS,
			       exchangeRate);
     //pick exchangeRate number of Small molecule from dest box
     state = prng.PickMol(kindL, kindIndexB, molIndexB, numInCavB, destBox);  
   }
   return state;
}

inline uint IdentityExchange::GetBoxPairAndMol
(const double subDraw, const double movPerc)
{
   uint state = mv::fail_state::NO_FAIL; 
   //deside to insert or remove the big molecule
   prng.PickBool(insertB, subDraw, movPerc);
   
#if ENSEMBLE == GEMC
   double density;
   double maxDens = 0.0;
   uint densB;
   //choose the sourceBox to be the dense phase
   for(uint b = 0; b < BOX_TOTAL; b++)
   {
     density = 0.0;
     for(uint k = 0; k < molLookRef.GetNumKind(); k++)
     {
       density += molLookRef.NumKindInBox(k, b) * boxDimRef.volInv[b] *
	 molRef.kinds[k].molMass;
     }

     if(density > maxDens)
     {
       maxDens = density;
       densB = b;
     }
   }

   //Pick box in dense phase
   sourceBox = densB; 
   //Pick the destination box
   prng.SetOtherBox(destBox, sourceBox);

#elif ENSEMBLE == GCMC
   sourceBox = 0;
   destBox = 1;
#endif

   if(insertB)
   {
     state = PickMolInCav();
   }
   else
   {
     state = ReplaceMolecule();
   }  
   
   if(state == mv::fail_state::NO_FAIL)
   {
     pStartA.clear();
     pStartB.clear();
     pStartA.resize(numInCavA);
     pStartB.resize(numInCavB);
     pLenA.clear();
     pLenB.clear();
     pLenA.resize(numInCavA);
     pLenB.resize(numInCavB);

     for(uint n = 0; n < numInCavA; n++)
     {
       pStartA[n] = pLenA[n] = 0;
       molRef.GetRangeStartLength(pStartA[n], pLenA[n], molIndexA[n]);
     }

     for(uint n = 0; n < numInCavB; n++)
     {
       pStartB[n] = pLenB[n] = 0;
       molRef.GetRangeStartLength(pStartB[n], pLenB[n], molIndexB[n]);
     }
   }

   return state;
}


inline uint IdentityExchange::Prep(const double subDraw, const double movPerc)
{
   uint state = GetBoxPairAndMol(subDraw, movPerc);
   if(state == mv::fail_state::NO_FAIL)
   {
     newMolA.clear();
     oldMolA.clear();
     newMolB.clear();
     oldMolB.clear();
     //transfering type A from source to dest
     for(uint n = 0; n < numInCavA; n++)
     {
       newMolA.push_back(cbmc::TrialMol(molRef.kinds[kindIndexA[n]], boxDimRef,
					destBox));
       oldMolA.push_back(cbmc::TrialMol(molRef.kinds[kindIndexA[n]], boxDimRef,
					sourceBox));
     }

     for(uint n = 0; n < numInCavB; n++)
     {
       //transfering type B from dest to source
       newMolB.push_back(cbmc::TrialMol(molRef.kinds[kindIndexB[n]], boxDimRef,
					sourceBox));
       oldMolB.push_back(cbmc::TrialMol(molRef.kinds[kindIndexB[n]], boxDimRef,
					destBox));
     }

     //set the old coordinate after unwrap them
     for(uint n = 0; n < numInCavA; n++)
     {
       XYZArray molA(pLenA[n]);
       coordCurrRef.CopyRange(molA, pStartA[n], 0, pLenA[n]);
       boxDimRef.UnwrapPBC(molA, sourceBox, comCurrRef.Get(molIndexA[n]));
       oldMolA[n].SetCoords(molA, 0);
       //set coordinate of moleA to newMolA, later it will shift to center
       newMolA[n].SetCoords(molA, 0);
     }

     for(uint n = 0; n < numInCavB; n++)
     {
       XYZArray molB(pLenB[n]);     
       coordCurrRef.CopyRange(molB, pStartB[n], 0, pLenB[n]);
       boxDimRef.UnwrapPBC(molB, destBox, comCurrRef.Get(molIndexB[n]));
       oldMolB[n].SetCoords(molB, 0);
       //set coordinate of moleB to newMolB, later it will shift to tempD
       newMolB[n].SetCoords(molB, 0);
     }

     XYZ axisD = boxDimRef.GetAxis(destBox);        
     XYZ axisS = boxDimRef.GetAxis(sourceBox);
 
     for(uint n = 0; n < numInCavB; n++)
     {
       if(kindIndexB[0] == kindL)
       {
	 //Inserting molB from destBox to the center of cavity in sourceBox
	 newMolB[n].SetSeed(center, rmax, true, true);
	 //perform rotational trial move in destBox for oldMolB
	 oldMolB[n].SetSeed(false, true);
       }
       else
       {
	 //Inserting molB from destBox to the cavity in sourceBox
	 newMolB[n].SetSeed(center, rmax, true, false);
	 //perform trial move in destBox for oldMolB
	 oldMolB[n].SetSeed(false, false);
       }
     }

     for(uint n = 0; n < numInCavA; n++)
     {
       if(kindIndexA[0] == kindL)
       {
	 //Inserting molA from sourceBox to destBox
	 newMolA[n].SetSeed(false, true);
	 //perform rotational trial move on COM for oldMolA
	 oldMolA[n].SetSeed(center, rmax, true, true);
       }
       else
       {
	 //Inserting molA from sourceBox to destBox
	 newMolA[n].SetSeed(false, false);
	 //perform trial move in cavity in sourceBox for oldMolA
	 oldMolA[n].SetSeed(center, rmax, true, false);
       }
     }
   }

   return state;
}


inline uint IdentityExchange::Transform()
{
  CalcTc();
  //Deleting MoleculeB from sourceBox
  for(uint n = 0; n < numInCavB; n++)
  {
    cellList.RemoveMol(molIndexB[n], destBox, coordCurrRef);
  }
  
  //Insert A to destBox
  for(uint n = 0; n < numInCavA; n++)
  {
    cellList.RemoveMol(molIndexA[n], sourceBox, coordCurrRef);
    //Transfer Type A to destBox
    molRef.kinds[kindIndexA[n]].BuildID(oldMolA[n], newMolA[n], molIndexA[n]);
    ShiftMol(n, sourceBox, destBox);
    cellList.AddMol(molIndexA[n], destBox, coordCurrRef);
  }

  //Deleting MoleculeA before inserting B
  for(uint n = 0; n < numInCavA; n++)
  {
    cellList.RemoveMol(molIndexA[n], destBox, coordCurrRef);
  }
  
  //Add moleculeB back to sourceBox
  for(uint n = 0; n < numInCavB; n++)
  {
    cellList.AddMol(molIndexB[n], destBox, coordCurrRef);
  }

  //Insert B in sourceBox
  for(uint n = 0; n < numInCavB; n++)
  {
    cellList.RemoveMol(molIndexB[n], destBox, coordCurrRef);
    //Transfer Type B to sourceBox
    molRef.kinds[kindIndexB[n]].BuildID(oldMolB[n], newMolB[n], molIndexB[n]);
    ShiftMol(n, destBox, sourceBox);
    cellList.AddMol(molIndexB[n], sourceBox, coordCurrRef);    
  }
  
  //Add moleculeA back to sourcebox
  for(uint n = 0; n < numInCavA; n++)
  {
    cellList.AddMol(molIndexA[n], destBox, coordCurrRef);
  }
  
  

  return mv::fail_state::NO_FAIL;
}

inline void IdentityExchange::CalcTc()
{
  W_tc = 1.0;
  if (ffRef.useLRC)
  {
    double delTC = 0.0;
    for (uint b = 0; b < BOX_TOTAL; ++b)
    {
      uint kCount[molRef.kindsCount];
      for (uint k = 0; k < molRef.kindsCount; ++k)
      {
	kCount[k] = molLookRef.NumKindInBox(k, b);
      }

      if (b == sourceBox)
      {
	kCount[kindIndexA[0]] -= numInCavA;
	kCount[kindIndexB[0]] += numInCavB;	   
      }
      else if (b == destBox)
      {
	kCount[kindIndexA[0]] += numInCavA;
	kCount[kindIndexB[0]] -= numInCavB;
      }
      tcNew[b].energy = calcEnRef.EnergyCorrection(b, kCount);
      delTC += tcNew[b].energy - sysPotRef.boxEnergy[b].tc;
    }
    W_tc = exp(-1.0 * ffRef.beta * delTC); 
  }
}
inline void IdentityExchange::CalcEn()
{   
   W_recip = 1.0;
   correct_oldA = 0.0, correct_newA = 0.0;
   self_oldA = 0.0, self_newA = 0.0;
   correct_oldB = 0.0, correct_newB = 0.0;
   self_oldB = 0.0, self_newB = 0.0;
   recipDest = 0.0, recipSource = 0.0;

   for(uint n = 0; n < numInCavA; n++)
   {
      correct_newA += calcEwald->SwapCorrection(newMolA[n]);
      correct_oldA += calcEwald->SwapCorrection(oldMolA[n]);
      self_newA += calcEwald->SwapSelf(newMolA[n]);
      self_oldA += calcEwald->SwapSelf(oldMolA[n]);
   }
   recipDest = calcEwald->SwapRecip(newMolA, oldMolB);

   for(uint n = 0; n < numInCavB; n++)
   {
     correct_newB += calcEwald->SwapCorrection(newMolB[n]);
     correct_oldB += calcEwald->SwapCorrection(oldMolB[n]);
     self_newB += calcEwald->SwapSelf(newMolB[n]);
     self_oldB += calcEwald->SwapSelf(oldMolB[n]);
   }
   recipSource = calcEwald->SwapRecip(newMolB, oldMolA);

   //need to contribute the self and correction energy 
   W_recip = exp(-1.0 * ffRef.beta * (recipSource + recipDest +
				      correct_newA - correct_oldA +
				      correct_newB - correct_oldB +
				      self_newA - self_oldA +
				      self_newB - self_oldB));
   
}

inline double IdentityExchange::GetCoeff() const
{
  double numTypeASource = molLookRef.NumKindInBox(kindIndexA[0], sourceBox);
  double numTypeADest = molLookRef.NumKindInBox(kindIndexA[0], destBox);
  double numTypeBSource = molLookRef.NumKindInBox(kindIndexB[0], sourceBox);
  double numTypeBDest = molLookRef.NumKindInBox(kindIndexB[0], destBox);
  double volSource = boxDimRef.volume[sourceBox];
  double volDest = boxDimRef.volume[destBox];
#if ENSEMBLE == GEMC
  if(insertB)
  {
    //kindA is the small molecule
    uint totMolInCav = molInCav[kindIndexA[0]].size();
    double ratioF =  Factorial(totMolInCav) /
      (Factorial(totMolInCav - exchangeRate) *
       Factorial(numTypeADest, exchangeRate));

    double ratioV = (volSource / volDest) * pow(volDest / volCav, exchangeRate);
    
    return ratioF * ratioV  * numTypeBDest / (numTypeBSource + 1.0);
  }
  else
  {
    //kindA is the big molecule
    uint totMolInCav = molInCav[kindIndexB[0]].size();
    double ratioF =  Factorial(totMolInCav) *
      Factorial(numTypeBDest - exchangeRate, exchangeRate) /
      Factorial(totMolInCav + exchangeRate);

    double ratioV = (volDest / volSource) * pow(volCav / volDest, exchangeRate);

    return ratioF * ratioV  * numTypeASource / (numTypeADest + 1.0);
  }
#elif ENSEMBLE == GCMC
  if(ffRef.isFugacity)
  {
    double delA = (BETA * molRef.kinds[kindIndexA[0]].chemPot * numInCavA);
    double insB = (BETA * molRef.kinds[kindIndexB[0]].chemPot * numInCavB);
    if(insertB)
    {
      //kindA is the small molecule
      uint totMolInCav = molInCav[kindIndexA[0]].size();
      double ratioF =  Factorial(totMolInCav) /
	Factorial(totMolInCav - exchangeRate);

      double ratioV = volSource / pow(volCav, exchangeRate);
      return (insB / delA) * ratioF * ratioV / (numTypeBSource + 1.0);
    }
    else
    {
      //kindA is the big molecule
      uint totMolInCav = molInCav[kindIndexB[0]].size();
      double ratioF = Factorial(totMolInCav) /
	Factorial(totMolInCav + exchangeRate);

      double ratioV = pow(volCav, exchangeRate) / volSource;
      return (insB / delA) * ratioF * ratioV * numTypeASource;
    }
  }
  else
  {
    double delA = exp(-BETA * molRef.kinds[kindIndexA[0]].chemPot * numInCavA);
    double insB = exp(BETA * molRef.kinds[kindIndexB[0]].chemPot * numInCavB);
    if(insertB)
    {
      //kindA is the small molecule
      uint totMolInCav = molInCav[kindIndexA[0]].size();
      double ratioF =  Factorial(totMolInCav) /
	Factorial(totMolInCav - exchangeRate);

      double ratioV = volSource / pow(volCav, exchangeRate);
      return delA * insB * ratioF * ratioV / (numTypeBSource + 1.0);
    }
    else
    {
      //kindA is the big molecule
      uint totMolInCav = molInCav[kindIndexB[0]].size();
      double ratioF = Factorial(totMolInCav) /
	Factorial(totMolInCav + exchangeRate);

      double ratioV = pow(volCav, exchangeRate) / volSource;
      return delA * insB * ratioF * ratioV * numTypeASource;
    }
  }
#endif
}

inline void IdentityExchange::ShiftMol(const uint n, const uint from,
				       const uint to)
{
  if(from == 0)
  {
    //Add type A to dest box
    newMolA[n].GetCoords().CopyRange(coordCurrRef, 0, pStartA[n], pLenA[n]);
    comCurrRef.SetNew(molIndexA[n], destBox);
    molLookRef.ShiftMolBox(molIndexA[n], sourceBox, destBox, kindIndexA[n]);
  }
  else
  {
    //Add type B to source box
    newMolB[n].GetCoords().CopyRange(coordCurrRef, 0, pStartB[n], pLenB[n]);
    comCurrRef.SetNew(molIndexB[n], sourceBox);
    molLookRef.ShiftMolBox(molIndexB[n], destBox, sourceBox, kindIndexB[n]);
  }
}

inline void IdentityExchange::RecoverMol(const uint n, const uint from,
					 const uint to)
{
  if(from == 0)
  {
    XYZArray molA(pLenA[n]);
    oldMolA[n].GetCoords().CopyRange(molA, 0, 0, pLenA[n]);
    boxDimRef.WrapPBC(molA, sourceBox);

    molA.CopyRange(coordCurrRef, 0, pStartA[n], pLenA[n]);
    comCurrRef.SetNew(molIndexA[n], sourceBox);
    molLookRef.ShiftMolBox(molIndexA[n], destBox, sourceBox, kindIndexA[n]);
  }
  else
  {
    XYZArray molB(pLenB[n]);
    oldMolB[n].GetCoords().CopyRange(molB, 0, 0, pLenB[n]);
    boxDimRef.WrapPBC(molB, destBox);

    molB.CopyRange(coordCurrRef, 0, pStartB[n], pLenB[n]);
    comCurrRef.SetNew(molIndexB[n], destBox);
    molLookRef.ShiftMolBox(molIndexB[n], sourceBox, destBox, kindIndexB[n]);
  }
}

inline double IdentityExchange::Factorial(const uint n) const
{
  double result = 1.0;
  for(uint i = 2; i <= n; i++)
  {
    result *= i;
  }

  return result;
}

inline double IdentityExchange::Factorial(const uint n, const uint count) const
{
  double result = 1.0;
  for(uint i = 1; i <= count; i++)
  {
    result *= n + i;
  }

  return result;
}

inline void IdentityExchange::Accept(const uint rejectState, const uint step)
{
   bool result;
   //If we didn't skip the move calculation
   if(rejectState == mv::fail_state::NO_FAIL)
   {
      double molTransCoeff = GetCoeff();  
      double Wrat = W_tc * W_recip;

      for(uint n = 0; n < numInCavA; n++)
      {
	Wrat *= newMolA[n].GetWeight() / oldMolA[n].GetWeight();
      }

      for(uint n = 0; n < numInCavB; n++)
      {
	Wrat *= newMolB[n].GetWeight() / oldMolB[n].GetWeight();
      }

      //std::cout << "Wrat: " << Wrat << std::endl;
      result = prng() < molTransCoeff * Wrat;
  
      if(result)
      {
         //Add tail corrections
         sysPotRef.boxEnergy[sourceBox].tc = tcNew[sourceBox].energy;
         sysPotRef.boxEnergy[destBox].tc = tcNew[destBox].energy;

         //Add rest of energy.
	 for(uint n = 0; n < numInCavB; n++)
	 {
	   sysPotRef.boxEnergy[sourceBox] += newMolB[n].GetEnergy();
	   sysPotRef.boxEnergy[destBox] -= oldMolB[n].GetEnergy();
	 }

	 for(uint n = 0; n < numInCavA; n++)
	 {
	   sysPotRef.boxEnergy[sourceBox] -= oldMolA[n].GetEnergy();
	   sysPotRef.boxEnergy[destBox] += newMolA[n].GetEnergy();
	 }
	 

	 //Add Reciprocal energy
	 sysPotRef.boxEnergy[sourceBox].recip += recipSource;
	 sysPotRef.boxEnergy[destBox].recip += recipDest;	 
	 //Add correction energy
	 sysPotRef.boxEnergy[sourceBox].correction -= correct_oldA;
	 sysPotRef.boxEnergy[sourceBox].correction += correct_newB;
	 sysPotRef.boxEnergy[destBox].correction += correct_newA;
	 sysPotRef.boxEnergy[destBox].correction -= correct_oldB;	 
	 //Add self energy
	 sysPotRef.boxEnergy[sourceBox].self -= self_oldA;
	 sysPotRef.boxEnergy[sourceBox].self += self_newB;
	 sysPotRef.boxEnergy[destBox].self += self_newA;
	 sysPotRef.boxEnergy[destBox].self -= self_oldB;
	 
	 for (uint b = 0; b < BOX_TOTAL; b++)
	 {
	    calcEwald->UpdateRecip(b);
	 }

	 //molA and molB already transfered to destBox and added to cellist

	 //Retotal
         sysPotRef.Total();
      }
      else
      {
	//transfer molA from destBox to source
	for(uint n = 0; n < numInCavA; n++)
	{
	  cellList.RemoveMol(molIndexA[n], destBox, coordCurrRef);
	  RecoverMol(n, sourceBox, destBox);
	  cellList.AddMol(molIndexA[n], sourceBox, coordCurrRef);
	}
	//transfer molB from sourceBox to dest
	for(uint n = 0; n < numInCavB; n++)
	{
	  cellList.RemoveMol(molIndexB[n], sourceBox, coordCurrRef);
	  RecoverMol(n, destBox, sourceBox);
	  cellList.AddMol(molIndexB[n], destBox, coordCurrRef);
	}
	
      }
   }
   else  //else we didn't even try because we knew it would fail
      result = false;

#if ENSEMBLE == GEMC
   subPick = mv::GetMoveSubIndex(mv::ID_EXCHANGE, sourceBox);
   moveSetRef.Update(result, subPick, step);
   subPick = mv::GetMoveSubIndex(mv::ID_EXCHANGE, destBox);
   moveSetRef.Update(result, subPick, step);
#elif ENSEMBLE == GCMC
   subPick = mv::GetMoveSubIndex(mv::ID_EXCHANGE);
   moveSetRef.Update(result, subPick, step);
#endif
}

#endif

#endif
