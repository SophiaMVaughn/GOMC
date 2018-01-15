#ifndef INTRAIDENTITYEXCHANGE_H
#define INTRAIDENTITYEXCHANGE_H

#include "MoveBase.h"
#include "cbmc/TrialMol.h"

using std::vector;

// Identity Exchange Move:
//
// insert the molecule A inside the cavity that has the size of the molecule B
// and vice versa. Only 1 trial must be perform for the first seed.
// Mohammad Soroush Barhaghi. July 2017

class IntraIdentityExchange : public MoveBase
{
 public:

   IntraIdentityExchange(System &sys, StaticVals const& statV) :
      ffRef(statV.forcefield), molLookRef(sys.molLookupRef),
      MoveBase(sys, statV), rmax(statV.mol.kinds[0].cavDim), cavA(3),
	invCavA(3), cavB(3), invCavB(3) 
      {
	enableID = ffRef.enableID;
	if(rmax.x >= rmax.y)
	  rmax.y = rmax.x;
	else
	  rmax.x = rmax.y;

	volCav = rmax.x * rmax.y * rmax.z; 

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

   void ShiftMol(const uint n, const bool kindA);
   void RecoverMol(const uint n, const bool kindA);
   uint PickMolInCav();
   double GetCoeff() const;
   //calculate factorial
   double Factorial(const uint n) const;
   //calculate ratio of factorial
   double Factorial(const uint n, const uint count) const;
   uint GetBoxPairAndMol(const double subDraw, const double movPerc);
   MolPick molPick;
   uint sourceBox;
   vector<uint> pStartA, pLenA, pStartB, pLenB;
   vector<uint> molIndexA, kindIndexA, molIndexB, kindIndexB;
   bool insertB, enableID;
   uint numInCavA, numInCavB, exchangeRate, kindS, kindL;
   uint numSCavA, numSCavB;

   double volCav;
   XYZ centerA, centerB, rmax;
   XYZArray cavA, invCavA, cavB, invCavB;
   double W_recip;
   double recipDiffA, recipDiffB;
   vector<cbmc::TrialMol> oldMolA, newMolA, oldMolB, newMolB;
   vector< vector<uint> > molInCav;
   MoleculeLookup & molLookRef;
   Forcefield const& ffRef;
};

inline uint IntraIdentityExchange::PickMolInCav()
{
   uint state = mv::fail_state::NO_FAIL;
   //pick a random location in dense phase
   XYZ axis = boxDimRef.GetAxis(sourceBox);
   XYZ temp(prng.randExc(axis.x), prng.randExc(axis.y), prng.randExc(axis.z));
   //Use to shift the new insterted molecule
   centerA = temp;
   //Pick random vector and find two vectors that are perpendicular to V1
   cavA.Set(0, prng.RandomUnitVect());
   cavA.GramSchmidt();
   //Calculate inverse matrix for cav
   cavA.TransposeMatrix(invCavA);
   //Find the molecule kind 0 in the cavity
   if(calcEnRef.FindMolInCavity(molInCav, centerA, rmax, invCavA,
				sourceBox, kindS, exchangeRate))
   {
     molIndexA.clear();
     kindIndexA.clear();
     molIndexB.clear();
     kindIndexB.clear();
     //printf("MolS in cav: %d.\n", molInCav[kindS].size());
     //Find the exchangeRate number of molecules kindS in cavity
     numInCavA = exchangeRate;
     numSCavA = molInCav[kindS].size();
     //printf("water Number: %d \n", numSCavA);
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
     
     //pick a molecule from kindL in same box
     numInCavB = 1;
     state = prng.PickMol(kindS, kindIndexB, molIndexB, numInCavB, sourceBox);
   }
   else
   {
     //reject the move
     state = mv::fail_state::NO_MOL_OF_KIND_IN_BOX;
   }

   //printf("S in Cavity: %d \n",  molInCav[kindS].size());

   if(state == mv::fail_state::NO_FAIL)
   {
     //Use to insert kindS in this centerB
     centerB = comCurrRef.Get(molIndexB[0]);
     //Set the V1 to the vector from first to last atom
     uint pStart = 0;
     uint pLen = 0;
     molRef.GetRangeStartLength(pStart, pLen, molIndexB[0]);
     if(pLen == 1)
       cavB.Set(0, prng.RandomUnitVect());
     else
     {
       uint pEnd = pStart + pLen -1;
       cavB.Set(0, boxDimRef.MinImage(coordCurrRef.Difference(pStart, pEnd),
				      sourceBox));
     }
     cavB.GramSchmidt();
     //Calculate inverse matrix for cav
     cavB.TransposeMatrix(invCavB);
     //find how many of KindS exist in this centerB (COM of kindL)
     calcEnRef.FindMolInCavity(molInCav, centerB, rmax, invCavB,
			       sourceBox, kindS, exchangeRate);
     numSCavB = molInCav[kindS].size();
   }
   return state;
}


inline uint IntraIdentityExchange::GetBoxPairAndMol
(const double subDraw, const double movPerc)
{
   uint state = mv::fail_state::NO_FAIL; 

#if ENSEMBLE == GCMC
   sourceBox = 0;
#else
   prng.PickBox(sourceBox, subDraw, movPerc);
#endif

   //Pick the exchange number of kindS in cavity and a molecule of kindL
   //kindA = kindS, kindB = kindL
   state = PickMolInCav();  
   
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


inline uint IntraIdentityExchange::Prep(const double subDraw,
					const double movPerc)
{
   uint state = GetBoxPairAndMol(subDraw, movPerc);
   if(state == mv::fail_state::NO_FAIL)
   {
     newMolA.clear();
     oldMolA.clear();
     newMolB.clear();
     oldMolB.clear();
     //transfering type A from source 
     for(uint n = 0; n < numInCavA; n++)
     {
       newMolA.push_back(cbmc::TrialMol(molRef.kinds[kindIndexA[n]], boxDimRef,
					sourceBox));
       oldMolA.push_back(cbmc::TrialMol(molRef.kinds[kindIndexA[n]], boxDimRef,
					sourceBox));
     }

     for(uint n = 0; n < numInCavB; n++)
     {
       //transfering type B from source
       newMolB.push_back(cbmc::TrialMol(molRef.kinds[kindIndexB[n]], boxDimRef,
					sourceBox));
       oldMolB.push_back(cbmc::TrialMol(molRef.kinds[kindIndexB[n]], boxDimRef,
					sourceBox));
     }

     //set the old coordinate after unwrap them
     for(uint n = 0; n < numInCavA; n++)
     {
       XYZArray molA(pLenA[n]);
       coordCurrRef.CopyRange(molA, pStartA[n], 0, pLenA[n]);
       boxDimRef.UnwrapPBC(molA, sourceBox, comCurrRef.Get(molIndexA[n]));
       oldMolA[n].SetCoords(molA, 0);
       //copy cavA matrix to slant the old trial of molA
       oldMolA[n].SetCavMatrix(cavA);
       //set coordinate of moleA to newMolA, later it will shift to centerB
       newMolA[n].SetCoords(molA, 0);
       //copy cavB matrix to slant the new trial of molA
       newMolA[n].SetCavMatrix(cavB);
     }

     for(uint n = 0; n < numInCavB; n++)
     {
       XYZArray molB(pLenB[n]);     
       coordCurrRef.CopyRange(molB, pStartB[n], 0, pLenB[n]);
       boxDimRef.UnwrapPBC(molB, sourceBox, comCurrRef.Get(molIndexB[n]));
       oldMolB[n].SetCoords(molB, 0);
       //copy cavB matrix to slant the old trial of molB
       oldMolB[n].SetCavMatrix(cavB);
       //set coordinate of moleB to newMolB, later it will shift to centerA
       newMolB[n].SetCoords(molB, 0);
       //copy cavA matrix to slant the new trial of molB
       newMolB[n].SetCavMatrix(cavA);
     }

     for(uint n = 0; n < numInCavB; n++)
     {
       //Inserting molB from centerB to the centerA 
       newMolB[n].SetSeed(centerA, rmax, true, true, true); 
       //perform rotational trial move for oldMolB
       oldMolB[n].SetSeed(centerB, rmax, true, true, true); 
     }

     for(uint n = 0; n < numInCavA; n++)
     {
       //Inserting molA from cavity(centerA) to the cavityB(centerB)
       newMolA[n].SetSeed(centerB, rmax, true, false, false);
       //perform trial move in cavity in sourceBox for oldMolA
       oldMolA[n].SetSeed(centerA, rmax, true, false, false); 
     }
   }

   return state;
}


inline uint IntraIdentityExchange::Transform()
{
  //Calc old energy before deleting
  for(uint n = 0; n < numInCavA; n++)
  {
    cellList.RemoveMol(molIndexA[n], sourceBox, coordCurrRef);
    molRef.kinds[kindIndexA[n]].BuildIDOld(oldMolA[n], molIndexA[n]);
  }

  //Calc old energy before deleting
  for(uint n = 0; n < numInCavB; n++)
  {
    cellList.RemoveMol(molIndexB[n], sourceBox, coordCurrRef);
    molRef.kinds[kindIndexB[n]].BuildIDOld(oldMolB[n], molIndexB[n]);
  }
  
  //Insert A to cavity of  center B
  for(uint n = 0; n < numInCavA; n++)
  {
    molRef.kinds[kindIndexA[n]].BuildIDNew(newMolA[n], molIndexA[n]);
    ShiftMol(n, true);
    cellList.AddMol(molIndexA[n], sourceBox, coordCurrRef);
  }

  //Insert B to cavity of  center A
  for(uint n = 0; n < numInCavB; n++)
  {
    molRef.kinds[kindIndexB[n]].BuildIDNew(newMolB[n], molIndexB[n]);
    ShiftMol(n, false);
    cellList.AddMol(molIndexB[n], sourceBox, coordCurrRef);    
  }
  
  return mv::fail_state::NO_FAIL;
}


inline void IntraIdentityExchange::CalcEn()
{   
   W_recip = 1.0;
   recipDiffA = 0.0, recipDiffB = 0.0;

   recipDiffA = calcEwald->SwapRecip(newMolA, oldMolA);
   recipDiffB = calcEwald->SwapRecip(newMolB, oldMolB);

   //No need to contribute the self and correction energy 
   W_recip = exp(-1.0 * ffRef.beta * (recipDiffA + recipDiffB));
   
}

inline double IntraIdentityExchange::GetCoeff() const
{
  double ratioF =  double(Factorial(numSCavA) * Factorial(numSCavB)) /
    double((Factorial(numSCavA - exchangeRate) *
	    Factorial(numSCavB + exchangeRate)));

  return ratioF;
}

inline void IntraIdentityExchange::ShiftMol(const uint n, const bool typeA)
{
  if(typeA)
  {
    //update coordinate of molecule typeA
    newMolA[n].GetCoords().CopyRange(coordCurrRef, 0, pStartA[n], pLenA[n]);
    comCurrRef.SetNew(molIndexA[n], sourceBox);
  }
  else
  {
    //update coordinate of molecule typeA
    newMolB[n].GetCoords().CopyRange(coordCurrRef, 0, pStartB[n], pLenB[n]);
    comCurrRef.SetNew(molIndexB[n], sourceBox);
  }
}

inline void IntraIdentityExchange::RecoverMol(const uint n, const bool typeA)
{
  if(typeA)
  {
    XYZArray molA(pLenA[n]);
    oldMolA[n].GetCoords().CopyRange(molA, 0, 0, pLenA[n]);
    boxDimRef.WrapPBC(molA, sourceBox);

    molA.CopyRange(coordCurrRef, 0, pStartA[n], pLenA[n]);
    comCurrRef.SetNew(molIndexA[n], sourceBox);
  }
  else
  {
    XYZArray molB(pLenB[n]);
    oldMolB[n].GetCoords().CopyRange(molB, 0, 0, pLenB[n]);
    boxDimRef.WrapPBC(molB, sourceBox);

    molB.CopyRange(coordCurrRef, 0, pStartB[n], pLenB[n]);
    comCurrRef.SetNew(molIndexB[n], sourceBox);
  }
}

inline double IntraIdentityExchange::Factorial(const uint n) const
{
  double result = 1.0;
  for(uint i = 2; i <= n; i++)
  {
    result *= i;
  }

  return result;
}

inline double IntraIdentityExchange::Factorial(const uint n,
					       const uint count) const
{
  double result = 1.0;
  for(uint i = 1; i <= count; i++)
  {
    result *= n + i;
  }

  return result;
}

inline void IntraIdentityExchange::Accept(const uint rejectState,
					  const uint step)
{
   bool result;
   //If we didn't skip the move calculation
   if(rejectState == mv::fail_state::NO_FAIL)
   {
      double molTransCoeff = GetCoeff();  
      double Wrat = W_recip;

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
         //Add rest of energy.
	 for(uint n = 0; n < numInCavB; n++)
	 {
	   sysPotRef.boxEnergy[sourceBox] += newMolB[n].GetEnergy();
	   sysPotRef.boxEnergy[sourceBox] -= oldMolB[n].GetEnergy();
	 }

	 for(uint n = 0; n < numInCavA; n++)
	 {
	   sysPotRef.boxEnergy[sourceBox] -= oldMolA[n].GetEnergy();
	   sysPotRef.boxEnergy[sourceBox] += newMolA[n].GetEnergy();
	 }
	 
	 //Add Reciprocal energy
	 sysPotRef.boxEnergy[sourceBox].recip += recipDiffA;
	 sysPotRef.boxEnergy[sourceBox].recip += recipDiffB;
	 
	 //Update reciprocal
	 calcEwald->UpdateRecip(sourceBox);

	 //molA and molB already added to cellist

	 //Retotal
         sysPotRef.Total();
      }
      else
      {
	//transfer molA 
	for(uint n = 0; n < numInCavA; n++)
	{
	  cellList.RemoveMol(molIndexA[n], sourceBox, coordCurrRef);
	  RecoverMol(n, true);
	  cellList.AddMol(molIndexA[n], sourceBox, coordCurrRef);
	}
	//transfer molB
	for(uint n = 0; n < numInCavB; n++)
	{
	  cellList.RemoveMol(molIndexB[n], sourceBox, coordCurrRef);
	  RecoverMol(n, false);
	  cellList.AddMol(molIndexB[n], sourceBox, coordCurrRef);
	}
	
      }
   }
   else  //else we didn't even try because we knew it would fail
      result = false;

   subPick = mv::GetMoveSubIndex(mv::INTRA_ID_EXCHANGE, sourceBox);
   moveSetRef.Update(result, subPick, step);

}

#endif
