#pragma once

#include "genesis/genesis.hpp"
#include <vector>

/**
 * This class organizes different taxon IDs (in a set of trees) representing the same taxon.
 * This is needed because the genesis toolkit assigns the IDs separately for each tree.
 */
class TaxonMapper {
public:
	size_t taxonEvalID(size_t treeID, size_t taxonID);
	TaxonMapper(Tree const &refTree, const std::string &evalTreesPath, std::vector<size_t> const &taxonIDs);
private:
	/**
	 * First index represents the evaluation tree, second index represents the taxon index in the reference tree.
	 * Entry represents the taxon index in the respective evaluation tree or std::numeric_limits<size_t>::max()
	 * if the taxon is not present in the evaluation tree.
	 */
	std::vector<std::vector<size_t> > taxonEvaluationMappings;
	size_t numNodesReference; /**< number of nodes in the reference tree */
	std::vector<size_t> numNodes; /**< number of nodes in the evaluation trees */
};

/**
 * Returns the ID of the taxon (which corresponds to the taxonID in the reference tree)
 * in the evaluation tree indexed by treeID.
 * @param treeID the ID of the evaluation tree
 * @param taxonID the ID of the taxon in the reference tree
 */
size_t TaxonMapper::taxonEvalID(size_t treeID, size_t taxonID) {
	if (taxonID > numNodesReference) {
		throw std::runtime_error("taxonID too large!");
	}
	if (treeID > taxonEvaluationMappings.size()) {
		throw std::runtime_error("treeID too large!");
	}
	if (taxonEvaluationMappings[treeID][taxonID] > numNodes[treeID]
			&& taxonEvaluationMappings[treeID][taxonID] < std::numeric_limits<size_t>::max()) {
		throw std::runtime_error(
				"Invalid taxon mapping entry! Got ID " + std::to_string(taxonEvaluationMappings[treeID][taxonID])
						+ " for the reference taxon " + std::to_string(taxonID) + ", but the tree "
						+ std::to_string(treeID) + " has only " + std::to_string(numNodes[treeID]) + " nodes.");
	}
	return taxonEvaluationMappings[treeID][taxonID];
}

/**
 * @param refTree the reference tree
 * @param evalTreesPaths path to the file containing the set of evaluation trees
 * @param taxonIDs the taxon IDs in the reference tree
 */
TaxonMapper::TaxonMapper(Tree const &refTree, const std::string &evalTreesPath, std::vector<size_t> const &taxonIDs) {
	numNodesReference = refTree.node_count();
	std::unordered_map<std::string, size_t> taxonToReferenceID;
	for (size_t i = 0; i < taxonIDs.size(); ++i) {
		taxonToReferenceID[refTree.node_at(taxonIDs[i]).data<DefaultNodeData>().name] = taxonIDs[i];
	}

	utils::InputStream instream(utils::make_unique < utils::FileInputSource > (evalTreesPath));
	auto itTree = NewickInputIterator(instream);
	while (itTree) { // iterate over the set of evaluation trees
		Tree const& tree = *itTree;

		numNodes.push_back(tree.node_count());
		std::vector<size_t> taxonMappings(refTree.node_count(), std::numeric_limits<size_t>::max());
		for (size_t j = 0; j < tree.node_count(); ++j) {
			if (tree.node_at(j).is_leaf()) {
				std::string leafName = tree.node_at(j).data<DefaultNodeData>().name;
				if (taxonToReferenceID.find(leafName) != taxonToReferenceID.end()) {
					taxonMappings[taxonToReferenceID[leafName]] = tree.node_at(j).index();
				}
			}
		}
		taxonEvaluationMappings.push_back(taxonMappings);

		++itTree;
	}
}
