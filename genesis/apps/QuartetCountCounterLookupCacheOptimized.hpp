#pragma once

#include "genesis/genesis.hpp"
#include <vector>
#include <cassert>
#include <algorithm>
#include <memory>
#include <unordered_map>
#include <stxxl/vector>

#include "QuartetCountConverter.hpp"
#include "TreeInformation.hpp"

using namespace genesis;
using namespace tree;
using namespace utils;
using namespace std;

// TODO FIXME: This is currently not working, for some reason.

class QuartetCounterLookupCacheOptimized {
public:
	QuartetCounterLookupCacheOptimized(Tree const &refTree, TreeSet const &evalTrees);
	//~QuartetCounterLookupCacheOptimized();
	std::tuple<size_t, size_t, size_t> countQuartetOccurrences(int aIdx, int bIdx, int cIdx, int dIdx);
private:
	void updateQuartets(const Tree &tree, size_t nodeIdx, std::vector<int> &eulerTourLeaves,
			std::vector<int> &linkToEulerLeafIndex);
	void updateQuartetsThreeLinks(size_t link1, size_t link2, size_t link3, const Tree &tree,
			std::vector<int> &eulerTourLeaves, std::vector<int> &linkToEulerLeafIndex);
	void updateQuartetsThreeClades(size_t startLeafIndexS1, size_t endLeafIndexS1, size_t startLeafIndexS2,
			size_t endLeafIndexS2, size_t startLeafIndexS3, size_t endLeafIndexS3, std::vector<int> &eulerTourLeaves);
	std::pair<size_t, size_t> subtreeLeafIndices(size_t linkIdx, const Tree &tree,
			std::vector<int> &linkToEulerLeafIndex);

	QuartetConverter quartetConverter;
	std::vector<QuartetCount> lookupTable;
};

void QuartetCounterLookupCacheOptimized::updateQuartetsThreeClades(size_t startLeafIndexS1, size_t endLeafIndexS1,
		size_t startLeafIndexS2, size_t endLeafIndexS2, size_t startLeafIndexS3, size_t endLeafIndexS3,
		std::vector<int> &eulerTourLeaves) {
	size_t aLeafIndex = startLeafIndexS1;
	size_t bLeafIndex = startLeafIndexS2;
	size_t cLeafIndex = startLeafIndexS3;

	while (aLeafIndex != endLeafIndexS1) {
		int a = eulerTourLeaves[aLeafIndex];
		size_t a2LeafIndex = (aLeafIndex + 1) % eulerTourLeaves.size();
		while (a2LeafIndex != endLeafIndexS1) {
			int a2 = eulerTourLeaves[a2LeafIndex];
			while (bLeafIndex != endLeafIndexS2) {
				int b = eulerTourLeaves[bLeafIndex];
				while (cLeafIndex != endLeafIndexS3) {
					int c = eulerTourLeaves[cLeafIndex];
					//#pragma omp critical
					{
						lookupTable[quartetConverter.quartetToNumber(a, a2, b, c)].addCount(a, a2, b, c, 1);
					}
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

// the leaf indizes are between [start,end)
std::pair<size_t, size_t> QuartetCounterLookupCacheOptimized::subtreeLeafIndices(size_t linkIdx, const Tree &tree,
		std::vector<int> &linkToEulerLeafIndex) {
	size_t outerLinkIdx = tree.link_at(linkIdx).outer().index();
	return {linkToEulerLeafIndex[linkIdx] % linkToEulerLeafIndex.size(), linkToEulerLeafIndex[outerLinkIdx] % linkToEulerLeafIndex.size()};
}

void QuartetCounterLookupCacheOptimized::updateQuartetsThreeLinks(size_t link1, size_t link2, size_t link3,
		const Tree &tree, std::vector<int> &eulerTourLeaves, std::vector<int> &linkToEulerLeafIndex) {
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

void QuartetCounterLookupCacheOptimized::updateQuartets(const Tree &tree, size_t nodeIdx,
		std::vector<int> &eulerTourLeaves, std::vector<int> &linkToEulerLeafIndex) {
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

QuartetCounterLookupCacheOptimized::QuartetCounterLookupCacheOptimized(Tree const &refTree, TreeSet const &evalTrees) :
		quartetConverter(refTree.node_count()) {
	int n = refTree.node_count();
	//quartetConverter.test();

	// initialize the lookup table.
	lookupTable.resize((n * (n - 1) * (n - 2) * (n - 3)) / 24);
	std::unordered_map<std::string, size_t> taxonToReferenceID;
	for (int i = 0; i < n; ++i) {
		taxonToReferenceID[refTree.node_at(i).data<DefaultNodeData>().name] = i;
	}

	unsigned int progress = 1;
	float onePercent = (float) evalTrees.size() / 100;

	for (size_t i = 0; i < evalTrees.size(); ++i) {
		size_t nEval = evalTrees[i].tree.node_count();

		// do an euler tour through the tree
		std::vector<int> eulerTourLeaves; // directly containing the mapped IDs from the reference
		std::vector<int> linkToEulerLeafIndex;
		linkToEulerLeafIndex.resize(evalTrees[i].tree.link_count());
		for (auto it : eulertour(evalTrees[i].tree)) {
			if (it.node().is_leaf()) {
				size_t leafIdx = it.node().index();
				eulerTourLeaves.push_back(
						taxonToReferenceID[evalTrees[i].tree.node_at(leafIdx).data<DefaultNodeData>().name]);
			}
			linkToEulerLeafIndex[it.link().index()] = eulerTourLeaves.size();
		}

#pragma omp parallel for
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

std::tuple<size_t, size_t, size_t> QuartetCounterLookupCacheOptimized::countQuartetOccurrences(int aIdx, int bIdx,
		int cIdx, int dIdx) {
	size_t quartetIndex = quartetConverter.quartetToNumber(aIdx, bIdx, cIdx, dIdx);
	size_t abCD = lookupTable[quartetIndex].getCount(aIdx, bIdx, cIdx, dIdx);
	size_t acBD = lookupTable[quartetIndex].getCount(aIdx, cIdx, bIdx, dIdx);
	size_t adBC = lookupTable[quartetIndex].getCount(aIdx, dIdx, bIdx, cIdx);
	return std::tuple<size_t, size_t, size_t>(abCD, acBD, adBC);
}
