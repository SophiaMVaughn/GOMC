#include "CBMC.h"
#include "MoleculeKind.h"
#include "MolSetup.h"
#include "DCLinear.h"
#include "DCGraph.h"
#include <vector>


namespace cbmc
{

CBMC* MakeCBMC(System& sys, const Forcefield& ff,
               const MoleculeKind& kind, const Setup& set)
{
/*
  std::vector<uint> bondCount(kind.NumAtoms(), 0);
  for (uint i = 0; i < kind.bondList.count; ++i)
  {
    bondCount[kind.bondList.part1[i]]++;
    bondCount[kind.bondList.part2[i]]++;
  }
  bool branched = false;
  for (uint i = 0; i < kind.NumAtoms(); ++i)
  {
    if (bondCount[i] > 2)
      branched = true;
  }
  bool cyclic = (kind.NumBonds() > kind.NumAtoms() - 1) ? true : false;
*/

  //Any molecule woth 3 atoms and more will be built in DCGraph
  if (kind.NumAtoms() > 2) {
    return new DCGraph(sys, ff, kind, set);
  } else {
    return new DCLinear(sys, ff, kind, set);
  }

}

}
