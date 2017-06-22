#pragma once

#include "genesis/genesis.hpp"
#include <vector>
#include <memory>
#include <algorithm>

using namespace genesis;
using namespace tree;
using namespace utils;

/**
 * Helper class that computes the lowest common ancestor of two nodes in respect to any given root node,
 * as well as the distance between any two nodes in number of edges.
 */
class TreeInformation {
public:
	TreeInformation();
	void init(Tree const& tree);
	/*
	 * @brief Returns the index of the lowest common ancestor of the nodes at index uIdx and vIdx with respect to root node at rootIdx.
	 */
	size_t lowestCommonAncestorIdx(size_t uIdx, size_t vIdx, size_t rootIdx);
	unsigned distanceInEdges(size_t uIdx, size_t vIdx);
	size_t getRootIdx();
private:
	size_t myRootIndex; /**< index of the root node */
	std::unique_ptr<RangeMinimumQuery> rrmq; /**< pointer to a RangeMinimumQuery */
	size_t rmqQueryCorrectOrder(size_t i, size_t j);
	std::vector<size_t> eulerTourNodes; /**< the nodes of the tree visited in an Eulerian tour order */
	std::vector<size_t> firstOccurrenceInEulerTour; /**< the nodes' first occurrence in an Eulerian tour */
	std::vector<size_t> dist_to_root; /**< distance of a given node to the root in terms of number of edges */
};

/**
 * Compute the distance in edges between two nodes, using their lowest common ancestor.
 * @param uIdx index of the first node in the tree
 * @param vIdx index of the second node in the tree
 */
unsigned TreeInformation::distanceInEdges(size_t uIdx, size_t vIdx) {
	size_t lcaIdx = lowestCommonAncestorIdx(uIdx, vIdx, myRootIndex);
	return dist_to_root[uIdx] + dist_to_root[vIdx] - 2 * dist_to_root[lcaIdx];
}

/**
 * Query the RangeMinimumQuery to find the index of the smallest entry between positions i and j in the euler tour levels.
 * This returns the ID of the lowest common ancestor of the nodes with IDs i and j, in respect to the tree's root node.
 * @param i the starting position
 * @param j the ending position
 */
size_t TreeInformation::rmqQueryCorrectOrder(size_t i, size_t j) {
	if (i <= j)
		return rrmq->query(i, j);
	else
		return rrmq->query(j, i);
}

/**
 * Return the ID of the root node in the tree.
 */
size_t TreeInformation::getRootIdx() {
	return myRootIndex;
}

/**
 * Return the ID of the lowest common ancestor of the nodes at uIdx and vIdx, using rootIdx as the root node ID.
 * @param uIdx ID of the node u
 * @param vIdx ID of the node v
 * @param rootIdx the ID of the node to be used as the root node
 */
size_t TreeInformation::lowestCommonAncestorIdx(size_t uIdx, size_t vIdx, size_t rootIdx) {
	size_t uEulerIdx = firstOccurrenceInEulerTour[uIdx];
	size_t vEulerIdx = firstOccurrenceInEulerTour[vIdx];
	size_t rootEulerIdx = firstOccurrenceInEulerTour[rootIdx];

	if (rootIdx == myRootIndex) {
		return eulerTourNodes[rmqQueryCorrectOrder(uEulerIdx, vEulerIdx)];
	} else { // take the "odd man out", see http://stackoverflow.com/questions/25371865/find-multiple-lcas-in-unrooted-tree
		size_t candidateOne = eulerTourNodes[rmqQueryCorrectOrder(uEulerIdx, vEulerIdx)];
		size_t candidateTwo = eulerTourNodes[rmqQueryCorrectOrder(uEulerIdx, rootEulerIdx)];
		size_t candidateThree = eulerTourNodes[rmqQueryCorrectOrder(vEulerIdx, rootEulerIdx)];
		if (candidateOne == candidateTwo) {
			return candidateThree;
		} else if (candidateOne == candidateThree) {
			return candidateTwo;
		} else {
			return candidateOne;
		}
	}
}

/**
 * @param tree the tree to build the TreeInformation for.
 */
void TreeInformation::init(Tree const &tree) {
	dist_to_root = node_path_length_vector(tree);
	// Arrays for LCA computation
	std::vector<int> eulerTourLevels;

	firstOccurrenceInEulerTour.resize(tree.node_count());
	std::fill(firstOccurrenceInEulerTour.begin(), firstOccurrenceInEulerTour.end(), std::numeric_limits<size_t>::max());

	for (auto it : eulertour(tree)) {
		eulerTourNodes.push_back(it.node().index());
		eulerTourLevels.push_back(dist_to_root[it.node().index()]);
		if (firstOccurrenceInEulerTour[it.node().index()] == std::numeric_limits<size_t>::max()) {
			firstOccurrenceInEulerTour[it.node().index()] = eulerTourNodes.size() - 1;
		}
	}
	myRootIndex = tree.root_node().index();

	rrmq = make_unique<RangeMinimumQuery>(eulerTourLevels);
}

/**
 * @param tree the tree to build the TreeInformation for.
 */
TreeInformation::TreeInformation() {
	myRootIndex = 0;
}
