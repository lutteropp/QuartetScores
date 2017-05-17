#pragma once

#include "genesis.hpp"
#include <vector>
#include <cassert>
#include <algorithm>
#include <memory>
#include "TreeInformation.hpp"
#include <unordered_map>
#include <stxxl/vector>
#include "TaxonMapper.hpp"

using namespace genesis;
using namespace tree;
using namespace utils;
using namespace std;

/**
 * Let n be the number of taxa in the reference tree and m be the number of evaluation trees.
 * Count quartet topologies for a small number of evaluation trees, with O(n*m) space and O(m) lookup cost.
 */
class QuartetCounterFewTrees {
public:
	QuartetCounterFewTrees(Tree const &refTree, TreeSet const &evalTrees);
	std::tuple<size_t, size_t, size_t> countQuartetOccurrences(size_t aIdx, size_t bIdx, size_t cIdx, size_t dIdx);
private:
	std::vector<std::unique_ptr<TreeInformation> > infos; /**> helper classes for lca and distance computation in each of the trees */
	std::unique_ptr<TaxonMapper> taxonMapper; /**> helper class for identifying taxa throughout the trees */
	std::vector<size_t> eulerTourLeaves; /**> the node IDs of the taxa, visited in an Eulerian tour order */
	TreeSet const &evaluationTrees; /**> The set of evaluation trees */
};

/**
 * @param refTree the reference tree
 * @param evalTrees the set of evaluation trees
 */
QuartetCounterFewTrees::QuartetCounterFewTrees(Tree const &refTree, TreeSet const &evalTrees) :
		evaluationTrees(evalTrees) {
	size_t n = refTree.node_count();
	std::unordered_map<std::string, size_t> taxonToReferenceID;
#pragma omp parallel for
	for (size_t i = 0; i < n; ++i) {
		taxonToReferenceID[refTree.node_at(i).data<DefaultNodeData>().name] = i;
	}
	for (size_t i = 0; i < evalTrees.size(); ++i) {
		infos.push_back(make_unique<TreeInformation>(evalTrees[i].tree));
	}

	// precompute subtree informations
	for (auto it : eulertour(refTree)) {
		if (it.node().is_leaf()) {
			eulerTourLeaves.push_back(it.node().index());
		}
	}

	// precompute taxon ID mappings
	taxonMapper = make_unique<TaxonMapper>(refTree, evalTrees, eulerTourLeaves);
	//std::cout << "Finished precomputing taxon ID mappings.\n";
}

/**
 * Count the occurrences of the quartet topologies ab|cd, ac|bd, and ad|bc in the evaluation trees.
 * This function uses pairwise distances in order to compute the topology of a quartet in a tree.
 * @param aIdx the ID of the taxon a in the reference tree
 * @param bIdx the ID of the taxon b in the reference tree
 * @param cIdx the ID of the taxon c in the reference tree
 * @param dIdx the ID of the taxon d in the reference tree
 */
std::tuple<size_t, size_t, size_t> QuartetCounterFewTrees::countQuartetOccurrences(size_t aIdx, size_t bIdx,
		size_t cIdx, size_t dIdx) {
	size_t q1, q2, q3;
	q1 = q2 = q3 = 0; // q1 = ab|cd, q2 = ac|bd, q3 = ad|bc
#pragma omp parallel for
	for (size_t i = 0; i < evaluationTrees.size(); ++i) {
		size_t aIdxTransformed = taxonMapper->taxonEvalID(i, aIdx);
		size_t bIdxTransformed = taxonMapper->taxonEvalID(i, bIdx);
		size_t cIdxTransformed = taxonMapper->taxonEvalID(i, cIdx);
		size_t dIdxTransformed = taxonMapper->taxonEvalID(i, dIdx);

		if (aIdxTransformed == std::numeric_limits<size_t>::max()
				|| bIdxTransformed == std::numeric_limits<size_t>::max()
				|| cIdxTransformed == std::numeric_limits<size_t>::max()
				|| dIdxTransformed == std::numeric_limits<size_t>::max()) {
			continue; // not all of a,b,c,d do exist in the current evaluation tree
		}

		// find topology
		size_t rootIdx = evaluationTrees[i].tree.root_node().index();
		size_t lca_ab = infos[i]->lowestCommonAncestorIdx(aIdxTransformed, bIdxTransformed, rootIdx);
		size_t lca_ac = infos[i]->lowestCommonAncestorIdx(aIdxTransformed, cIdxTransformed, rootIdx);
		size_t lca_ad = infos[i]->lowestCommonAncestorIdx(aIdxTransformed, dIdxTransformed, rootIdx);
		size_t lca_bc = infos[i]->lowestCommonAncestorIdx(bIdxTransformed, cIdxTransformed, rootIdx);
		size_t lca_bd = infos[i]->lowestCommonAncestorIdx(bIdxTransformed, dIdxTransformed, rootIdx);
		size_t lca_cd = infos[i]->lowestCommonAncestorIdx(cIdxTransformed, dIdxTransformed, rootIdx);

		if (infos[i]->distanceInEdges(lca_ab, lca_cd) > infos[i]->distanceInEdges(lca_ac, lca_bd)
				&& infos[i]->distanceInEdges(lca_ab, lca_cd) > infos[i]->distanceInEdges(lca_ad, lca_bc)) {
			q1++; // ab|cd
		} else if (infos[i]->distanceInEdges(lca_ac, lca_bd) > infos[i]->distanceInEdges(lca_ab, lca_cd)
				&& infos[i]->distanceInEdges(lca_ac, lca_bd) > infos[i]->distanceInEdges(lca_ad, lca_bc)) {
			q2++; // ac|bd
		} else if (infos[i]->distanceInEdges(lca_ad, lca_bc) > infos[i]->distanceInEdges(lca_ab, lca_cd)
				&& infos[i]->distanceInEdges(lca_ad, lca_bc) > infos[i]->distanceInEdges(lca_ac, lca_bd)) {
			q3++; // ad|bc
		} // else, we have a multifurcation and the quartet has none of these three topologies.
	}
	return std::tuple<size_t, size_t, size_t>(q1, q2, q3);
}
