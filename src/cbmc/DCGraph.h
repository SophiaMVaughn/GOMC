/*******************************************************************************
GPU OPTIMIZED MONTE CARLO (GOMC) 1.0 (Serial version)
Copyright (C) 2015  GOMC Group
A copy of the GNU General Public License can be found in the COPYRIGHT.txt
along with this program, also can be found at <http://www.gnu.org/licenses/>.
********************************************************************************/
#ifndef DCGRAPH_H
#define DCGRAPH_H
#include "CBMC.h"
#include "DCComponent.h"
#include "DCData.h"
#include <vector>
#include <utility>

/*CBMC graph of a branched/cyclic molecule
* The Decoupled/Coupled CBMC algorithm is represented by
* traversing a spanning tree of the graph.
*/

namespace cbmc
{

   class DCComponent;

   class DCGraph : public CBMC
   {
   public:
      DCGraph(System& sys, const Forcefield& ff,
         const MoleculeKind& kind, const Setup& set);

      void Build(TrialMol& oldMol, TrialMol& newMol, uint molIndex);
      void BuildIDNew(TrialMol& newMol, uint molIndex);
      void BuildIDOld(TrialMol& oldMol, uint molIndex);
      void Regrowth(TrialMol& oldMol, TrialMol& newMol, uint molIndex);
      ~DCGraph();

   private:
      //Store edge's atom that are connected to node and has more than 1 bond
      //Each edge is a node as well
      struct Edge
      {
	//destination is partner of the atom.
	uint destination;
	DCComponent* component;
        Edge(uint d, DCComponent* c) : destination(d), component(c) {}
      };
      
      //Store the branching atom and all Atoms that are connected to this
      //branching atom
      struct Node
      {
	//starting will be initialized with DCFreeHedron
	DCComponent *starting;
	//all the atom that are connected to this node and has more than 1 bond
	//will be in edges and initialized with DCLinkedHedron
	std::vector<Edge> edges;
      };
      
      DCComponent *idExchange;

      DCData data;
      //Branching atoms, C in CH4
      std::vector<Node> nodes;
      std::vector<Edge> fringe;
      std::vector<bool> visited;
   };
}


#endif
