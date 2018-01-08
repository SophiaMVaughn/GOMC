#define _USE_MATH_DEFINES 
#include <math.h> 
#include "DCRotateCOM.h" 
#include "DCData.h" 
#include "TrialMol.h" 
#include "MolSetup.h" 
#include "Forcefield.h" 
#include "PRNG.h" 
#include "NumLib.h" 
/*
This file is only called from IdentityExchange and IntraIdentityExchange file.
This file depend on ensemble kind, swap the COM of small molecules with one 
large molecule.

GCMC:
Deleting large molecule: Place a cavity on the selected Large Molecule with    
                         center at L-Mol COM and orientation of its backbone
BuildOld: perform nLJ trial around the backbone of Large Molecule.
BuildNew: perform fLJ for COM of S-Mol, for each fLJ, perform nLJ to rotate 
          around COM in sphere.

Insert large molecule: in a random oriented cavity and random center 
BuildOld: perform fLJ for COM of S-Mol inside the cavity, for each fLJ,
          perform nLJ to rotate around its COM in sphere.
BuildNew: insert the COM of L-Mol to the center of cavity and perform nLJ trial
          around the axis(cavity largest length-c)
###############################################################################
GEMC: Each molecule type has NEW and OLD which is similar to GCMC. However, to 
      insert or delete the Lrge molecule from dense box.
Deleting large molecule:Place a cavity on the selected Large Molecule with 
                        center at L-Mol COM and orientation of its backbone
BuildOldL: perform nLJ trial around the backbone of Large Molecule.
BuildNewL: perform fLJ for COM of L-Mol, for each fLJ, perform nLJ to rotate 
           around COM in sphere.
BuildOldS: perform fLJ for COM of S-Mol, for each fLJ, perform nLJ to rotate 
           around its COM in sphere.
BuildNewS: perform fLJ for COM of S-Mol, for each fLJ, perform nLJ to rotate 
           around COM in sphere.

Insert large molecule: In a random oriented cavity and random center.
BuildOldS: perform fLJ for COM of S-Mol inside the cavity, for each fLJ,
           perform nLJ to rotate around its COM in sphere.
BuildNewS: perform fLJ for COM of S-Mol, for each fLJ, perform nLJ to rotate 
           around COM in sphere.
BuildOldL: perform fLJ for COM of L-Mol, for each fLJ, perform nLJ to rotate 
           around COM in sphere.
BuildNewL: Insert the COM of L-Mol to the center of cavity and perform nLJ trial
           around the axis(cavity largest length-c)
*/ 

namespace cbmc 
{ 
 
   DCRotateCOM::DCRotateCOM(DCData* data) 
     : data(data) {} 
 

 
   void DCRotateCOM::PrepareNew(TrialMol& newMol, uint molIndex) 
   { 
     newMol.SetWeight(1.0);
     atomNumber = newMol.GetCoords().Count();
     //old center of mass
     oldCOM = newMol.GetCOM();
   } 

   void DCRotateCOM::PickTransferCOMNew(TrialMol& newMol, uint molIndex)
   {
     PRNG& prng = data->prng;

     if(newMol.SeedFix())
     {
       if(newMol.HasCav())
	 COM = newMol.GetSeed();
       else
	 prng.FillWithRandom(COM, data->axes.GetAxis(newMol.GetBox()));
     }
     else
     {
       //new center of mass that need to be transfered 
       if(newMol.HasCav())
       {
	 //Pick molecule in cav dimension
	 prng.FillWithRandomInCavity(COM, newMol.GetRmax());
	 //rotate using cavity matrix
	 COM = newMol.Transform(COM);
	 //add center 
	 COM += newMol.GetSeed();
       }
       else
	 prng.FillWithRandom(COM, data->axes.GetAxis(newMol.GetBox())); 
     }

     XYZ diff = COM - oldCOM;
     
     for(uint p = 0; p < atomNumber; p++)
     {
       newMol.SetAtomCoords(p, newMol.AtomPosition(p) + diff);
     }

     oldCOM = COM;
   }

 
   void DCRotateCOM::PrepareOld(TrialMol& oldMol, uint molIndex) 
   { 
     oldMol.SetWeight(1.0);
     atomNumber = oldMol.GetCoords().Count();
     //old center of mass
     oldCOM = oldMol.GetCOM();
     COM = oldCOM;
   } 

   void DCRotateCOM::PickTransferCOMOld(TrialMol& oldMol, uint molIndex)
   {
     PRNG& prng = data->prng;
     //new center of mass that need to be transfered
     if(oldMol.HasCav())
     {
       //Pick molecule in cav dimension
       prng.FillWithRandomInCavity(COM, oldMol.GetRmax());
       //rotate using cavity matrix
       COM = oldMol.Transform(COM);
       //add center 
       COM += oldMol.GetSeed();
     }
     else
       prng.FillWithRandom(COM, data->axes.GetAxis(oldMol.GetBox())); 

     XYZ diff = COM - oldCOM;
     
     for(uint p = 0; p < atomNumber; p++)
     {
       oldMol.SetAtomCoords(p, oldMol.AtomPosition(p) + diff);
     }

     oldCOM = COM;
   }
 
   void DCRotateCOM::BuildNew(TrialMol& newMol, uint molIndex) 
   { 
      PRNG& prng = data->prng; 
      const CalculateEnergy& calc = data->calc; 
      const EwaldCached *calcEwald = data->calcEwald; 
      const Forcefield& ff = data->ff; 
      uint nLJTrials = data->nLJTrialsNth; 
      uint fLJTrials = data->nLJTrialsFirst;
      uint totalTrials = data->totalTrials;
      double* ljWeights = data->ljWeightsT; 
      double* inter = data->interT; 
      double* real = data->realT; 
      RotationMatrix spin;
 
      std::fill_n(inter, totalTrials, 0.0); 
      std::fill_n(real, totalTrials, 0.0); 
      std::fill_n(ljWeights, totalTrials, 0.0); 
 
      XYZArray *multiPosRotions;  
      multiPosRotions = new XYZArray[atomNumber];
      for(uint i = 0; i < atomNumber; ++i) 
      { 
        multiPosRotions[i] = XYZArray(totalTrials); 
      }

      if(atomNumber == 1)
      {
	nLJTrials = 1;
	totalTrials = fLJTrials;
      }

      if(newMol.SeedFix())
      {
	fLJTrials = 1;
	totalTrials = nLJTrials;
      }
 

      for (uint p = 0; p < fLJTrials; ++p)
      {
	//Pick a new position for COM and transfer the molecule
	PickTransferCOMNew(newMol, molIndex);
	//get info about existing geometry 
	newMol.ShiftBasis(COM); 
	const XYZ center = COM; 
	uint index = p * nLJTrials;
	
	for (uint a = 0; a < atomNumber; ++a) 
	{ 
	  multiPosRotions[a].Set(index, newMol.AtomPosition(a)); 
	  multiPosRotions[a].Add(index, -center);
	} 
	
	//Rotational trial the molecule around COM
	for (uint r = nLJTrials; r-- > 0;) 
	{ 
	  //convert chosen torsion to 3D positions 
	  spin = RotationMatrix::UniformRandom(prng(), prng(), prng());

	  for (uint a = 0; a < atomNumber; ++a) 
	  { 
	    //find positions 
	    multiPosRotions[a].Set(index + r,
				   spin.Apply(multiPosRotions[a][index]));
	    multiPosRotions[a].Add(index + r, center); 
	  } 
	}
      } 

      for (uint a = 0; a < atomNumber; ++a) 
      { 
         data->axes.WrapPBC(multiPosRotions[a], newMol.GetBox());  
	 calc.ParticleInter(inter, real, multiPosRotions[a], a,  
			    molIndex, newMol.GetBox(), totalTrials); 
      }  
 
      double stepWeight = 0.0; 
      for (uint lj = 0; lj < totalTrials; ++lj) 
      { 
	ljWeights[lj] = exp(-ff.beta * (inter[lj] + real[lj])); 
         stepWeight += ljWeights[lj]; 
      } 
      uint winner = prng.PickWeighted(ljWeights, totalTrials, stepWeight); 

      for(uint a = 0; a < atomNumber; ++a) 
      { 
         newMol.AddAtom(a, multiPosRotions[a][winner]); 
      } 
 
      newMol.AddEnergy(Energy(0.0, 0.0, inter[winner], real[winner], 0.0, 0.0,
			      0.0)); 
      newMol.MultWeight(stepWeight / totalTrials); 

      delete[] multiPosRotions; 
   } 
 
   void DCRotateCOM::BuildOld(TrialMol& oldMol, uint molIndex) 
   { 
      PRNG& prng = data->prng; 
      const CalculateEnergy& calc = data->calc; 
      const EwaldCached * calcEwald = data->calcEwald; 
      const Forcefield& ff = data->ff; 
      uint nLJTrials = data->nLJTrialsNth; 
      uint fLJTrials = data->nLJTrialsFirst;
      uint totalTrials = data->totalTrials;
      double* ljWeights = data->ljWeightsT; 
      double* inter = data->interT; 
      double* real = data->realT; 
      RotationMatrix spin;

      std::fill_n(inter, totalTrials, 0.0); 
      std::fill_n(real, totalTrials, 0.0);  
      std::fill_n(ljWeights, totalTrials, 0.0); 
      
      XYZArray *multiPosRotions;  
      multiPosRotions = new XYZArray[atomNumber];
      for(uint i = 0; i < atomNumber; ++i) 
      { 
        multiPosRotions[i] = XYZArray(totalTrials); 
      }

      if(atomNumber == 1)
      {
	nLJTrials = 1;
	totalTrials = fLJTrials;
      }

      if(oldMol.SeedFix())
      {
	fLJTrials = 1;
	totalTrials = nLJTrials;
      }

      const XYZ orgCenter = COM;

      for (uint p = 0; p < fLJTrials; ++p)
      {
	//First trial is current configuration
	//get info about existing geometry 
	oldMol.ShiftBasis(COM); 
	const XYZ center = COM; 
	uint index = p * nLJTrials;	
 
	for (uint a = 0; a < atomNumber; ++a) 
	{ 
	  //get position and shift to origin 
	  multiPosRotions[a].Set(index, oldMol.AtomPosition(a)); 
	  multiPosRotions[a].Add(index, -center); 
	}  

	//Rotational trial the molecule around COM
	for (uint r = nLJTrials; r-- > 0;)
	{ 
	  if((index + r) == 0)
	    continue;

	  //convert chosen torsion to 3D positions 
	  spin =  RotationMatrix::UniformRandom(prng(), prng(), prng()); 

	  for (uint a = 0; a < atomNumber; ++a) 
	  { 
	    //find positions 
	    multiPosRotions[a].Set(index + r,
				   spin.Apply(multiPosRotions[a][index])); 
            multiPosRotions[a].Add(index + r, center); 
	  } 
	}
	//Pick a new position for COM and transfer the molecule
	PickTransferCOMOld(oldMol, molIndex);
      }
 
      for (uint a = 0; a < atomNumber; ++a) 
      { 
         multiPosRotions[a].Add(0, orgCenter); 
         data->axes.WrapPBC(multiPosRotions[a], oldMol.GetBox()); 
	 calc.ParticleInter(inter, real, multiPosRotions[a], a, 
                            molIndex, oldMol.GetBox(), totalTrials);
      } 

      double stepWeight = 0.0;  
      for (uint lj = 0; lj < totalTrials; ++lj) 
      { 
	stepWeight += exp(-ff.beta * (inter[lj] + real[lj])); 
      } 

      for (uint a = 0; a < atomNumber; ++a)
      { 
         oldMol.AddAtom(a, multiPosRotions[a][0]); 
      } 
 
      oldMol.AddEnergy(Energy(0.0, 0.0, inter[0], real[0], 0.0, 0.0, 0.0)); 
      oldMol.MultWeight(stepWeight / totalTrials); 

      delete[] multiPosRotions; 
   }
 
}              
