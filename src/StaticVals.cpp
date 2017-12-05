#include "StaticVals.h"
#include "ConfigSetup.h" //For types directly read from config. file
#include "Setup.h" //For source of setup data.

void StaticVals::Init(Setup & set, System& sys)
{
   //Standard inits
   simEventFreq.Init(set.config.sys.step);
   forcefield.Init(set);
   mol.Init(set, forcefield, sys);
#ifndef VARIABLE_VOLUME
   boxDimensions.Init(set.config.in.restart, 
		      set.config.sys.volume, set.pdb.cryst, forcefield.rCut, 
		      forcefield.rCutSq);
#endif
#ifndef VARIABLE_PARTICLE_NUMBER
   molLookup.Init(mol, set.pdb.atoms);
#endif
   InitMovePercents(set.config.sys.moves);
#if ENSEMBLE == GEMC || ENSEMBLE == NPT
   kindOfGEMC = set.config.sys.gemc.kind;
   pressure = set.config.sys.gemc.pressure;
   fixVolBox0 = set.config.sys.volume.cstVolBox0;
#endif
}

void StaticVals::InitOver(Setup & set, System& sys)
{
   mol.Init(set, forcefield, sys);
}

void StaticVals::InitMovePercents(config_setup::MovePercents const& perc)
{
   totalPerc = 0.0;
   for (uint m = 0; m < mv::MOVE_KINDS_TOTAL; ++m)
   {
      switch(m)
      {
      case mv::DISPLACE:
	 movePerc[m] = perc.displace; break;
      case mv::ROTATE:
	 movePerc[m] = perc.rotate; break;
      case mv::INTRA_SWAP:
	 movePerc[m] = perc.intraSwap; break;
      case mv::INTRA_ID_EXCHANGE :
	movePerc[m] = perc.intraIdExchange; break;
      case mv::REGROWTH:
	movePerc[m] = perc.regrowth; break;
#ifdef VARIABLE_VOLUME
      case mv::VOL_TRANSFER :
	 movePerc[m] = perc.volume; break;
#endif
#ifdef VARIABLE_PARTICLE_NUMBER
#if ENSEMBLE == GEMC || ENSEMBLE == GCMC
      case mv::MOL_TRANSFER :
	 movePerc[m] = perc.transfer; break;	 
      case mv::ID_EXCHANGE :
	movePerc[m] = perc.idExchange; break; 
#endif
#endif
      default:
	 movePerc[m] = 0.0; break;
      }
      totalPerc += movePerc[m];
   } 
   for (uint m = 0; m < mv::MOVE_KINDS_TOTAL; m++)
      movePerc[m] /= totalPerc;
   totalPerc = 1.0;
}
