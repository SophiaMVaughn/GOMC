#include "DCGraph.h"
#include "DCFreeHedron.h"
#include "DCLinkedHedron.h"
#include "DCRotateCOM.h"
#include "MolSetup.h"
#include "Setup.h"
#include "MoleculeKind.h"
#include <cassert>
#include <map>

namespace cbmc
{
DCGraph::DCGraph(System& sys, const Forcefield& ff,
                 const MoleculeKind& kind, const Setup& set)
  : data(sys, ff, set)
{
  using namespace mol_setup;
  MolMap::const_iterator it = set.mol.kindMap.find(kind.name);
  assert(it != set.mol.kindMap.end());
  const MolKind setupKind = it->second;

  idExchange = new DCRotateCOM(&data);

  std::vector<uint> atomToNode(setupKind.atoms.size(), 0);
  std::vector<uint> bondCount(setupKind.atoms.size(), 0);
  //Count the number of bonds for each atom 
  for (uint b = 0; b < setupKind.bonds.size(); ++b)
  {
    const Bond& bond = setupKind.bonds[b];
    ++bondCount[bond.a0];
    ++bondCount[bond.a1];
  }

  //Find the node (number of bound > 1)
  //Construct the starting node (DCFreeHedron)
  //Construct the Linking node (DCLinkHedron)
  for (uint atom = 0; atom < setupKind.atoms.size(); ++atom)
  {
    if (bondCount[atom] < 2)
    {
      //e.g. H in CH4.
      atomToNode[atom] = -1;
    }
    else
    {
      //Get the information of other Atoms that are bonded to the atom
      std::vector<Bond> bonds = AtomBonds(setupKind, atom);
      atomToNode[atom] = nodes.size();
      //Add atom to the node list and initialize it with DCFreeHedron, atom and
      // the first partner of the atom
      nodes.push_back(Node());
      Node& node = nodes.back();
      //Atoms will be build from prev(a1) to focus (atom)
      node.starting = new DCFreeHedron(&data, setupKind, atom,
                                       bonds[0].a1);
      //Loop through all the bonds
      for (uint i = 0; i < bonds.size(); ++i)
      {
        uint partner = bonds[i].a1;
	//Check if the partner of the atom is a node(has more than 1 bound)
        if(bondCount[partner] == 1)
        {
          continue;
        }
	//Add partner to the edge list of node and initialize it with partner
	//and the atom in DCLinkedHedron
	//Atoms will be build from prev(atom) to focus (partner)
        Edge e = Edge(partner, new DCLinkedHedron(&data, setupKind, partner,
						  atom));
        node.edges.push_back(e);
      }
    }
  }

  //reassign destination values from atom indices to node indices
  for (uint i = 0; i < nodes.size(); ++i)
  {
    for (uint j = 0; j < nodes[i].edges.size(); ++j)
    {
      uint& dest = nodes[i].edges[j].destination;
      dest = atomToNode[dest];
      assert(dest != -1);
    }
  }
}

void DCGraph::Build(TrialMol& oldMol, TrialMol& newMol, uint molIndex)
{
  //Randomely pick a node to call DCFreeHedron on it
  uint current = data.prng.randIntExc(nodes.size());
  visited.assign(nodes.size(), false);
  //Visiting the node
  visited[current] = true;
  //Copy the edges of the node to fringe
  fringe = nodes[current].edges;
  DCComponent* comp = nodes[current].starting;
  //Call DCFreeHedron to build all Atoms connected to the node
  comp->PrepareNew(newMol, molIndex);
  comp->BuildNew(newMol, molIndex);
  comp->PrepareOld(oldMol, molIndex);
  comp->BuildOld(oldMol, molIndex);
  //Advance along edges, building as we go
  while(!fringe.empty())
  {
    //Randomely pick one of the edges connected to node
    uint pick = data.prng.randIntExc(fringe.size());
    comp = fringe[pick].component;
    //Call DCLinkedHedron and build all Atoms connected to selected edge
    comp->PrepareNew(newMol, molIndex);
    comp->BuildNew(newMol, molIndex);
    comp->PrepareOld(oldMol, molIndex);
    comp->BuildOld(oldMol, molIndex);

    //Travel to new node, remove traversed edge
    //Current node is the edge that we picked
    current = fringe[pick].destination;
    //Remove the edge that we visited
    fringe[pick] = fringe.back();
    fringe.pop_back();
    //Visiting the node
    visited[current] = true;

    //Add edges to unvisited nodes
    for(uint i = 0; i < nodes[current].edges.size(); ++i)
    {
      Edge& e = nodes[current].edges[i];
      if(!visited[e.destination])
      {
        fringe.push_back(e);
      }
    }
  }
}

void DCGraph::BuildIDNew(TrialMol& newMol, uint molIndex)
{
  idExchange->PrepareNew(newMol, molIndex);
  idExchange->BuildNew(newMol, molIndex);
}

void DCGraph::BuildIDOld(TrialMol& oldMol, uint molIndex)
{
  idExchange->PrepareOld(oldMol, molIndex);
  idExchange->BuildOld(oldMol, molIndex);
}

void DCGraph::Regrowth(TrialMol& oldMol, TrialMol& newMol, uint molIndex)
{
  Build(oldMol, newMol, molIndex);
}

DCGraph::~DCGraph()
{
  for(uint v = 0; v < nodes.size(); ++v)
  {
    Node& node = nodes[v];
    delete node.starting;
    for(uint e = 0; e < node.edges.size(); ++ e)
    {
      delete node.edges[e].component;
    }
  }

  delete idExchange;
}

}
