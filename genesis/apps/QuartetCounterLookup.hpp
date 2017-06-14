#pragma once

#define USE_STXXL FALSE

#include "genesis.hpp"
#include <vector>
#include <cassert>
#include <algorithm>
#include <memory>
#include "TreeInformation.hpp"
#include <unordered_map>
#include <cstdint>
#if USE_STXXL
#include <stxxl/vector>
#endif

using namespace genesis;
using namespace tree;
using namespace utils;
using namespace std;

#define CO(a,b,c,d) (a) * n_cube + (b) * n_square + (c) * n + (d)

typedef uint32_t cint;

/**
 * Let n be the number of taxa in the reference tree.
 * Count occurrences of quartet topologies in the set of evaluation trees using a O(n^4) lookup table with O(1) lookup cost.
 */

class QuartetCounterLookup {
public:
	QuartetCounterLookup(const Tree &refTree, const TreeSet &evalTrees);
	~QuartetCounterLookup();
	std::tuple<cint, cint, cint> countQuartetOccurrences(size_t aIdx, size_t bIdx, size_t cIdx, size_t dIdx);
private:
	cint lookupQuartetCount(size_t aIdx, size_t bIdx, size_t cIdx, size_t dIdx);
	void countQuartets(const TreeSet &evalTrees, const std::unordered_map<std::string, size_t> &taxonToReferenceID);
	void updateQuartets(const Tree &tree, size_t nodeIdx, std::vector<int> &eulerTourLeaves,
			std::vector<int> &linkToEulerLeafIndex);
	void updateQuartetsThreeLinks(size_t link1, size_t link2, size_t link3, const Tree &tree,
			std::vector<int> &eulerTourLeaves, std::vector<int> &linkToEulerLeafIndex);
	void updateQuartetsThreeClades(size_t startLeafIndexS1, size_t endLeafIndexS1, size_t startLeafIndexS2,
			size_t endLeafIndexS2, size_t startLeafIndexS3, size_t endLeafIndexS3, std::vector<int> &eulerTourLeaves);
	std::pair<size_t, size_t> subtreeLeafIndices(size_t linkIdx, const Tree &tree,
			std::vector<int> &linkToEulerLeafIndex);

#if USE_STXXL
	stxxl::VECTOR_GENERATOR<cint>::result lookupTable; /**> O(n^4) lookup table storing the count of each quartet topology */
#else
	std::vector<cint> lookupTable; /**> O(n^4) lookup table storing the count of each quartet topology */
#endif
	size_t n; /**> number of taxa in the reference tree */
	size_t n_square; /**> n*n */
	size_t n_cube; /**> n*n*n */
	std::vector<size_t> refIdToLookupID;
};

// This is only needed in case the STXXL is used
QuartetCounterLookup::~QuartetCounterLookup() {
#if USE_STXXL
	while (!lookupTable.empty()) {
		lookupTable.pop_back();
	}
#endif
}

/**
 * Update the quartet topology counts for quartets  {a.b.c.d} where a,b \in S_1, c \in S_2, and d \in S_3.
 * @param startLeafIndexS1 the first index in eulerTourLeaves that corresponds to a leaf in subtree S_1
 * @param endLeafIndexS1 the last index in eulerTourLeaves that corresponds to a leaf in subtree S_1
 * @param startLeafIndexS2 the first index in eulerTourLeaves that corresponds to a leaf in subtree S_2
 * @param endLeafIndexS2 the last index in eulerTourLeaves that corresponds to a leaf in subtree S_2
 * @param startLeafIndexS3 the first index in eulerTourLeaves that corresponds to a leaf in subtree S_3
 * @param endLeafIndexS3 the last index in eulerTourLeaves that corresponds to a leaf in subtree S_3
 * @param eulerTourLeaves the leaves' IDs of the tree traversed in an euler tour order
 */
void QuartetCounterLookup::updateQuartetsThreeClades(size_t startLeafIndexS1, size_t endLeafIndexS1,
		size_t startLeafIndexS2, size_t endLeafIndexS2, size_t startLeafIndexS3, size_t endLeafIndexS3,
		std::vector<int> &eulerTourLeaves) {
	size_t aLeafIndex = startLeafIndexS1;
	size_t bLeafIndex = startLeafIndexS2;
	size_t cLeafIndex = startLeafIndexS3;

	while (aLeafIndex != endLeafIndexS1) {
		size_t a = eulerTourLeaves[aLeafIndex];
		size_t a2LeafIndex = (aLeafIndex + 1) % eulerTourLeaves.size();
		while (a2LeafIndex != endLeafIndexS1) {
			size_t a2 = eulerTourLeaves[a2LeafIndex];
			while (bLeafIndex != endLeafIndexS2) {
				size_t b = eulerTourLeaves[bLeafIndex];
				while (cLeafIndex != endLeafIndexS3) {
					size_t c = eulerTourLeaves[cLeafIndex];
#pragma omp atomic
					//lookupTable[CO(refIdToLookupID[a], refIdToLookupID[a2], refIdToLookupID[b], refIdToLookupID[c])]++;
					lookupTable[CO(a, a2, b, c)]++;
					cLeafIndex = (cLeafIndex + 1) % eulerTourLeaves.size();
				}
				bLeafIndex = (bLeafIndex + 1) % eulerTourLeaves.size();
				cLeafIndex = startLeafIndexS3;
			}
			a2LeafIndex = (a2LeafIndex + 1) % eulerTourLeaves.size();
			bLeafIndex = startLeafIndexS2;
			cLeafIndex = startLeafIndexS3;
		}
		aLeafIndex = (aLeafIndex + 1) % eulerTourLeaves.size();
		bLeafIndex = startLeafIndexS2;
		cLeafIndex = startLeafIndexS3;
	}
}

/**
 * Return a pair <start, end> representing the leaf indices in the Euler tour within the subtree induced by the genesis TreeLink with ID linkIdx.
 * The leaf indices are between [start,end), this means they include the start index but not the end index.
 * @param linkIdx the ID of the TreeLink from genesis
 * @param tree the tree
 * @param linkToEulerLeafIndex Mapping of each link in the tree to indices in the euler tour;
 * 	needed for determining first and last index of leaves belonging to a subtree.
 */
std::pair<size_t, size_t> QuartetCounterLookup::subtreeLeafIndices(size_t linkIdx, const Tree &tree,
		std::vector<int> &linkToEulerLeafIndex) {
	size_t outerLinkIdx = tree.link_at(linkIdx).outer().index();
	return {linkToEulerLeafIndex[linkIdx] % linkToEulerLeafIndex.size(), linkToEulerLeafIndex[outerLinkIdx] % linkToEulerLeafIndex.size()};
}

/**
 * Given the genesis links to the tree subtrees induced by an inner node, update the quartet topology counts of all quartets
 * {a,b,c,d} for which a and b are in the same subtree, c is in another subtree, and d is in the remaining subtree.
 * @param link1 link ID to the first subtree
 * @param link2 link ID to the second subtree
 * @param link3 link ID to the third subtree
 * @param tree the evaluation tree
 * @param eulerTourLeaves the leaves' IDs of the tree traversed in an euler tour order
 * @param linkToEulerLeafIndex Mapping of each link in the tree to indices in the euler tour;
 * 	needed for determining first and last index of leaves belonging to a subtree.
 */
void QuartetCounterLookup::updateQuartetsThreeLinks(size_t link1, size_t link2, size_t link3, const Tree &tree,
		std::vector<int> &eulerTourLeaves, std::vector<int> &linkToEulerLeafIndex) {
	std::pair<size_t, size_t> subtree1 = subtreeLeafIndices(link1, tree, linkToEulerLeafIndex);
	std::pair<size_t, size_t> subtree2 = subtreeLeafIndices(link2, tree, linkToEulerLeafIndex);
	std::pair<size_t, size_t> subtree3 = subtreeLeafIndices(link3, tree, linkToEulerLeafIndex);

	size_t startLeafIndexS1 = subtree1.first % eulerTourLeaves.size();
	size_t endLeafIndexS1 = subtree1.second % eulerTourLeaves.size();
	size_t startLeafIndexS2 = subtree2.first % eulerTourLeaves.size();
	size_t endLeafIndexS2 = subtree2.second % eulerTourLeaves.size();
	size_t startLeafIndexS3 = subtree3.first % eulerTourLeaves.size();
	size_t endLeafIndexS3 = subtree3.second % eulerTourLeaves.size();

	updateQuartetsThreeClades(startLeafIndexS1, endLeafIndexS1, startLeafIndexS2, endLeafIndexS2, startLeafIndexS3,
			endLeafIndexS3, eulerTourLeaves);
	updateQuartetsThreeClades(startLeafIndexS2, endLeafIndexS2, startLeafIndexS1, endLeafIndexS1, startLeafIndexS3,
			endLeafIndexS3, eulerTourLeaves);
	updateQuartetsThreeClades(startLeafIndexS3, endLeafIndexS3, startLeafIndexS1, endLeafIndexS1, startLeafIndexS2,
			endLeafIndexS2, eulerTourLeaves);
}

/**
 * An inner node in a bifurcating tree induces three subtrees S_1, S_2, and S_3.
 * Given an evaluation tree and an inner node, update the quartet topology counts of all quartets
 * {a,b,c,d} for which a and b are in the same subtree, c is in another subtree, and d is in the remaining subtree.
 * @param tree the evaluation tree
 * @param nodeIdx ID of an inner node in the evaluation tree
 * @param eulerTourLeaves the leaves' IDs of the tree traversed in an euler tour order
 * @param linkToEulerLeafIndex Mapping of each link in the tree to indices in the euler tour;
 * 	needed for determining first and last index of leaves belonging to a subtree.
 */
void QuartetCounterLookup::updateQuartets(const Tree &tree, size_t nodeIdx, std::vector<int> &eulerTourLeaves,
		std::vector<int> &linkToEulerLeafIndex) {
	// get taxa from subtree clades at nodeIdx
	std::vector<size_t> subtreeLinkIndices;
	const TreeLink* actLinkPtr = &tree.node_at(nodeIdx).link();
	subtreeLinkIndices.push_back(actLinkPtr->index());
	while (subtreeLinkIndices[0] != actLinkPtr->next().index()) {
		actLinkPtr = &actLinkPtr->next();
		subtreeLinkIndices.push_back(actLinkPtr->index());
	}

	for (size_t i = 0; i < subtreeLinkIndices.size(); ++i) {
		for (size_t j = i + 1; j < subtreeLinkIndices.size(); ++j) {
			for (size_t k = j + 1; k < subtreeLinkIndices.size(); ++k) {
				size_t link1 = subtreeLinkIndices[i];
				size_t link2 = subtreeLinkIndices[j];
				size_t link3 = subtreeLinkIndices[k];
				updateQuartetsThreeLinks(link1, link2, link3, tree, eulerTourLeaves, linkToEulerLeafIndex);
			}
		}
	}
}

void QuartetCounterLookup::countQuartets(const TreeSet &evalTrees,
		const std::unordered_map<std::string, size_t> &taxonToReferenceID) {
	unsigned int progress = 1;
	float onePercent = (float) evalTrees.size() / 100;

#pragma omp parallel for
	for (size_t i = 0; i < evalTrees.size(); ++i) {
		size_t nEval = evalTrees[i].tree.node_count();

		// do an euler tour through the tree
		std::vector<int> eulerTourLeaves; // directly containing the mapped IDs from the reference
		std::vector<int> linkToEulerLeafIndex;
		linkToEulerLeafIndex.resize(evalTrees[i].tree.link_count());
		for (auto it : eulertour(evalTrees[i].tree)) {
			if (it.node().is_leaf()) {
				size_t leafIdx = it.node().index();
				/*eulerTourLeaves.push_back(
						taxonToReferenceID.at(evalTrees[i].tree.node_at(leafIdx).data<DefaultNodeData>().name));*/
				eulerTourLeaves.push_back(
										refIdToLookupID[taxonToReferenceID.at(evalTrees[i].tree.node_at(leafIdx).data<DefaultNodeData>().name)]);
			}
			linkToEulerLeafIndex[it.link().index()] = eulerTourLeaves.size();
		}

		for (size_t j = 0; j < nEval; ++j) {
			if (!evalTrees[i].tree.node_at(j).is_leaf()) {
				updateQuartets(evalTrees[i].tree, j, eulerTourLeaves, linkToEulerLeafIndex);
			}
		}

		if (i > progress * onePercent) {
			std::cout << "Counting quartets... " << progress << "%" << std::endl;
#pragma omp atomic
			progress++;
		}
	}
}

/**
 * @param refTree the reference tree
 * @para evalTrees the set of evaluation trees
 */
QuartetCounterLookup::QuartetCounterLookup(Tree const &refTree, TreeSet const &evalTrees) {
	std::unordered_map<std::string, size_t> taxonToReferenceID;
	refIdToLookupID.resize(refTree.node_count());
	n = 0;
	for (auto it : eulertour(refTree)) {
		if (it.node().is_leaf()) {
			taxonToReferenceID[it.node().data<DefaultNodeData>().name] = it.node().index();
			refIdToLookupID[it.node().index()] = n;
			n++;
		}
	}
	n_square = n * n;
	n_cube = n_square * n;
	// initialize the lookup table.
	lookupTable.resize(n * n * n * n);
	//#pragma omp parallel for

	countQuartets(evalTrees, taxonToReferenceID);
	std::cout << "lookup table size: " << lookupTable.size() << "\n";
}

/**
 * Returns the count of the quartet topology ab|cd in the evaluation trees
 * @param aIdx ID of taxon a
 * @param bIdx ID of taxon b
 * @param cIdx ID of taxon c
 * @param dIdx ID of taxon d
 */
cint QuartetCounterLookup::lookupQuartetCount(size_t aIdx, size_t bIdx, size_t cIdx, size_t dIdx) {
	aIdx = refIdToLookupID[aIdx];
	bIdx = refIdToLookupID[bIdx];
	cIdx = refIdToLookupID[cIdx];
	dIdx = refIdToLookupID[dIdx];
	return lookupTable[CO(aIdx, bIdx, cIdx, dIdx)] + lookupTable[CO(aIdx, bIdx, dIdx, cIdx)]
			+ lookupTable[CO(bIdx, aIdx, cIdx, dIdx)] + lookupTable[CO(bIdx, aIdx, dIdx, cIdx)];
}

/**
 * Returns the counts of the quartet topologies ab|cd, ac|bd, and ad|bc in the evaluation trees
 * @param aIdx ID of taxon a
 * @param bIdx ID of taxon b
 * @param cIdx ID of taxon c
 * @param dIdx ID of taxon d
 */
std::tuple<cint, cint, cint> QuartetCounterLookup::countQuartetOccurrences(size_t aIdx, size_t bIdx, size_t cIdx,
		size_t dIdx) {
	cint abCD = lookupQuartetCount(aIdx, bIdx, cIdx, dIdx);
	cint acBD = lookupQuartetCount(aIdx, cIdx, bIdx, dIdx);
	cint adBC = lookupQuartetCount(aIdx, dIdx, bIdx, cIdx);
	return std::tuple<cint, cint, cint>(abCD, acBD, adBC);
}
