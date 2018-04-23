#include <cassert>
#include "DCLinear.h"
#include "DCSingle.h"
#include "DCOnSphere.h"
#include "DCRotateCOM.h"

using namespace cbmc;

DCLinear::DCLinear(System& sys, const Forcefield& ff,
                   const MoleculeKind& kind, const Setup& set) :
  data(sys, ff, set)
{
  mol_setup::MolMap::const_iterator it = set.mol.kindMap.find(kind.name);
  assert(it != set.mol.kindMap.end());
  const mol_setup::MolKind setupKind = it->second;
  uint size = kind.NumAtoms();
  atomSize = size;

  idExchange = new DCRotateCOM(&data);
  //First atom of the molecule
  forward.push_back(new DCSingle(&data, 0));
  backward.push_back(new DCSingle(&data, size - 1));
  //second atom of the molecule
  if(atomSize > 1) {
    forward.push_back(new DCOnSphere(&data, setupKind, 1, 0));
    backward.push_back(new DCOnSphere(&data, setupKind, size - 2, size - 1));
  }
}

DCLinear::~DCLinear()
{
  for(uint i = 0; i < forward.size(); ++i)
  {
    delete forward[i];
    delete backward[i];
  }
  delete idExchange;
}

void DCLinear::Build(TrialMol& oldMol, TrialMol& newMol, uint molIndex)
{
  std::vector<DCComponent*>& comps =data.prng.randInt(1) ? forward : backward;
  for(uint i = 0; i < comps.size(); ++i) {
    comps[i]->PrepareNew(newMol, molIndex);
    comps[i]->BuildNew(newMol, molIndex);
  }

  for(uint i = 0; i < comps.size(); ++i) {
    comps[i]->PrepareOld(oldMol, molIndex);
    comps[i]->BuildOld(oldMol, molIndex);
  }  
}

void DCLinear::BuildIDNew(TrialMol& newMol, uint molIndex)
{
  idExchange->PrepareNew(newMol, molIndex);
  idExchange->BuildNew(newMol, molIndex);
}

void DCLinear::BuildIDOld(TrialMol& oldMol, uint molIndex)
{
  idExchange->PrepareOld(oldMol, molIndex);
  idExchange->BuildOld(oldMol, molIndex);
}

void DCLinear::BuildOld2(TrialMol& oldMol, uint molIndex)
{
  std::vector<DCComponent*>& comps =data.prng.randInt(1) ? forward : backward;
  for(uint i = 0; i < comps.size(); ++i)
  {
    comps[i]->PrepareOld(oldMol, molIndex);
    comps[i]->BuildOld(oldMol, molIndex);
  }
}

void DCLinear::BuildNew2(TrialMol& newMol, uint molIndex)
{
  std::vector<DCComponent*>& comps =data.prng.randInt(1) ? forward : backward;
  for(uint i = 0; i < comps.size(); ++i)
  {
    comps[i]->PrepareNew(newMol, molIndex);
    comps[i]->BuildNew(newMol, molIndex);
  }
}

void DCLinear::BuildGrowOld(TrialMol& oldMol, uint molIndex)
{
  std::vector<DCComponent*>& comps = forward;
  for(uint i = 0; i < comps.size(); ++i)
  {
    comps[i]->PrepareOld(oldMol, molIndex);
    comps[i]->BuildOld(oldMol, molIndex);
  }
}

void DCLinear::BuildGrowNew(TrialMol& newMol, uint molIndex)
{
  std::vector<DCComponent*>& comps = forward;
  for(uint i = 0; i < comps.size(); ++i)
  {
    comps[i]->PrepareNew(newMol, molIndex);
    comps[i]->BuildNew(newMol, molIndex);
  }

}

void DCLinear::Regrowth(TrialMol& oldMol, TrialMol& newMol, uint molIndex)
{
  //perform Intra-Swap move within the same box
  if(atomSize < 2)
  {
    return Build(oldMol, newMol, molIndex);
  }
  else
  {
    //we only have two atoms in molecule: atom 0, 1
    uint fix = data.prng.randInt(1);
    uint grow = 1 - fix;
    //If fix == 0, forward (build atom 1), else backward (build atom 0)
    std::vector<DCComponent*>& comps = fix ? backward : forward;
    
    //copy the coordinate of the fix part
    newMol.AddAtom(fix, oldMol.AtomPosition(fix));
    oldMol.ConfirmOldAtom(fix);
    //build the second atom
    comps[1]->PrepareNew(newMol, molIndex);
    comps[1]->BuildNew(newMol, molIndex);
    comps[1]->PrepareOld(oldMol, molIndex);
    comps[1]->BuildOld(oldMol, molIndex);
  }
}
