#pragma once

#include "genesis/genesis.hpp"
#include <vector>
#include <cassert>
#include <algorithm>
#include <memory>
#include <utility>
#include <unordered_map>
#include "TreeInformation.hpp"
#include <unordered_map>
#include "QuartetCounterLookup.hpp"
#include "QuartetCounterFewTrees.hpp"
#include <unistd.h>
#include "QuartetCountCounterLookupCacheOptimized.hpp"

using namespace genesis;
using namespace tree;
using namespace utils;

// TODO: FIXME: The QuartetCounterLookup solution works, the QuartetCounterLookupCacheOptimized solution does not work.
#define CACHE_OPTIMIZED_LOOKUP FALSE

/**
 * Hash a std::pair.
 */
struct pairhash {
public:
	template<typename T, typename U>
	std::size_t operator()(const std::pair<T, U> &x) const {
		return std::hash<T>()(x.first) ^ std::hash<U>()(x.second);
	}
};

/**
 * Compute LQ-, QP-, and EQP-IC support scores for quartets.
 */
class QuartetScoreComputer {
public:
	QuartetScoreComputer(Tree const &refTree, const std::string &evalTreesPath, bool verboseOutput);
	std::vector<double> getLQICScores();
	std::vector<double> getQPICScores();
	std::vector<double> getEQPICScores();
private:
	size_t countEvalTrees(const std::string& evalTreesPath);

	double log_score(size_t q1, size_t q2, size_t q3);
	void computeQuartetScoresBifurcating();

	void computeQuartetScoresBifurcatingQuartets();
	std::unordered_map<std::pair<size_t, size_t>, std::tuple<size_t, size_t, size_t>, pairhash> countBuffer;

	void computeQuartetScoresMultifurcating();
	std::pair<size_t, size_t> nodePairForQuartet(size_t aIdx, size_t bIdx, size_t cIdx, size_t dIdx);
	void processNodePair(size_t uIdx, size_t vIdx);
	std::tuple<size_t, size_t, size_t> countQuartetOccurrences(size_t aIdx, size_t bIdx, size_t cIdx, size_t dIdx);

	std::pair<size_t, size_t> subtreeLeafIndices(size_t linkIdx);

	Tree referenceTree; /**< the reference tree */
	size_t rootIdx; /**< ID of the genesis root node in the reference tree */
	//TreeSet evaluationTrees; /**< a set of evaluation trees */
	bool useLookupTable; /**< toggles usage of a lookup table versus usage of the distance-based method for few trees */

	bool verbose;

	std::unique_ptr<TreeInformation> informationReferenceTree;
	std::vector<double> LQICScores;
	std::vector<double> QPICScores;
	std::vector<double> EQPICScores;

	std::vector<size_t> eulerTourLeaves;
	std::vector<size_t> linkToEulerLeafIndex;

#if CACHE_OPTIMIZED_LOOKUP
	std::unique_ptr<QuartetCounterLookupCacheOptimized> quartetCounterLookup;
#else
	std::unique_ptr<QuartetCounterLookup> quartetCounterLookup;
#endif
	std::unique_ptr<QuartetCounterFewTrees> quartetCounterFewTrees;
};

/**
 * Retrieve the IDs of the node pair {u,v} that induces the metaquartet which induces the quartet {a,b,c,d}.
 * @param aIdx ID of the taxon a in the reference tree
 * @param bIdx ID of the taxon b in the reference tree
 * @param cIdx ID of the taxon c in the reference tree
 * @param dIdx ID of the taxon d in the reference tree
 */
std::pair<size_t, size_t> QuartetScoreComputer::nodePairForQuartet(size_t aIdx, size_t bIdx, size_t cIdx, size_t dIdx) {
	size_t uIdx = informationReferenceTree->lowestCommonAncestorIdx(aIdx, bIdx, rootIdx);
	size_t vIdx = informationReferenceTree->lowestCommonAncestorIdx(cIdx, dIdx, rootIdx);

	if (uIdx == informationReferenceTree->lowestCommonAncestorIdx(aIdx, bIdx, vIdx)) {
		vIdx = informationReferenceTree->lowestCommonAncestorIdx(cIdx, dIdx, uIdx);
	} else {
		uIdx = informationReferenceTree->lowestCommonAncestorIdx(aIdx, bIdx, vIdx);
	}
	return std::pair<size_t, size_t>(uIdx, vIdx);
}

/**
 * Return the computed LQ-IC support scores.
 */
std::vector<double> QuartetScoreComputer::getLQICScores() {
	return LQICScores;
}

/**
 * Return the computed QP-IC support scores.
 */
std::vector<double> QuartetScoreComputer::getQPICScores() {
	return QPICScores;
}

/**
 * Return the computed EQP-IC support scores.
 */
std::vector<double> QuartetScoreComputer::getEQPICScores() {
	return EQPICScores;
}

/**
 * Compute qic = 1 + p_q1 * log(p_q1) + p_q2 * log(p_q2) + p_q3 * log(p_q3)
 * Take into account corner cases when one or more values are zero.
 * Let ab|cd be the topology of the quartet {a,b,c,d} in the reference tree.
 * @param q1 count of the quartet topology ab|cd
 * @param q2 count of the quartet topology ac|bd
 * @param q3 count of the quartet topology ad|bc
 */
double QuartetScoreComputer::log_score(size_t q1, size_t q2, size_t q3) {
	if (q1 == 0 && q2 == 0 && q3 == 0)
		return 0;

	size_t sum = q1 + q2 + q3;
	double p_q1 = (double) q1 / sum;
	double p_q2 = (double) q2 / sum;
	double p_q3 = (double) q3 / sum;

	// qic = 1 + p_q1 * log(p_q1) + p_q2 * log(p_q2) + p_q3 * log(p_q3);
	double qic = 1;
	if (p_q1 != 0)
		qic += p_q1 * log(p_q1) / log(3);
	if (p_q2 != 0)
		qic += p_q2 * log(p_q2) / log(3);
	if (p_q3 != 0)
		qic += p_q3 * log(p_q3) / log(3);

	if (q1 < q2 || q1 < q3) {
		return qic * -1;
	} else {
		return qic;
	}
}

/**
 * @brief For a given pair of nodes (and their LCA), return the indices of the two links of those
 * nodes that point towards the respective other node.
 *
 * For example, given node A and node B, the first index in the returned std::pair is the index
 * of the link of node A that points in direction of node B, and the second index in the std::pair
 * is the index of the link of node B that points towards node A.
 */
std::pair<size_t, size_t> get_path_inner_links(TreeNode const& u, TreeNode const& v, TreeNode const& lca) {
	if (&u == &v) {
		throw std::runtime_error("No inner links on a path between a node and itself.");
	}
	if (&u == &lca) {
		TreeLink const* u_link = nullptr;
		for (auto it : path_set(v, u, lca)) {
			if (it.is_lca())
				continue;
			u_link = &it.link().outer();
		}
		assert(u_link != nullptr);
		return {u_link->index(), v.primary_link().index()};
	}
	if (&v == &lca) {
		TreeLink const* v_link = nullptr;
		for (auto it : path_set(u, v, lca)) {
			if (it.is_lca())
				continue;
			v_link = &it.link().outer();
		}
		assert(v_link != nullptr);
		return {u.primary_link().index(), v_link->index()};
	}
	return {u.primary_link().index(), v.primary_link().index()};
}

/**
 * Count the occurrences of the quartet topologies ab|cd, ac|bd, and ad|bc in the evaluation trees.
 * @param aIdx ID of the taxon a
 * @param bIdx ID of the taxon b
 * @param cIdx ID of the taxon c
 * @param dIdx ID of the taxon d
 */
std::tuple<size_t, size_t, size_t> QuartetScoreComputer::countQuartetOccurrences(size_t aIdx, size_t bIdx, size_t cIdx,
		size_t dIdx) {
	if (useLookupTable) {
		return quartetCounterLookup->countQuartetOccurrences(aIdx, bIdx, cIdx, dIdx);
	} else {
		return quartetCounterFewTrees->countQuartetOccurrences(aIdx, bIdx, cIdx, dIdx);
	}
}

/**
 * Compute the LQ-IC, QP-IC, and EQP-IC scores for a bifurcating reference tree, iterating over quartets instead of node pairs.
 */
void QuartetScoreComputer::computeQuartetScoresBifurcatingQuartets() {
	// Process all quartets
	//#pragma omp parallel for schedule(dynamic)
	for (size_t uLeafIdx = 0; uLeafIdx < eulerTourLeaves.size(); uLeafIdx++) {
		for (size_t vLeafIdx = uLeafIdx + 1; vLeafIdx < eulerTourLeaves.size(); vLeafIdx++) {
			for (size_t wLeafIdx = vLeafIdx + 1; wLeafIdx < eulerTourLeaves.size(); wLeafIdx++) {
				for (size_t zLeafIdx = wLeafIdx + 1; zLeafIdx < eulerTourLeaves.size(); zLeafIdx++) {
					size_t uIdx = eulerTourLeaves[uLeafIdx];
					size_t vIdx = eulerTourLeaves[vLeafIdx];
					size_t wIdx = eulerTourLeaves[wLeafIdx];
					size_t zIdx = eulerTourLeaves[zLeafIdx];
					// find topology ab|cd of {u,v,w,z}
					size_t lca_uv = informationReferenceTree->lowestCommonAncestorIdx(uIdx, vIdx, rootIdx);
					size_t lca_uw = informationReferenceTree->lowestCommonAncestorIdx(uIdx, wIdx, rootIdx);
					size_t lca_uz = informationReferenceTree->lowestCommonAncestorIdx(uIdx, zIdx, rootIdx);
					size_t lca_vw = informationReferenceTree->lowestCommonAncestorIdx(vIdx, wIdx, rootIdx);
					size_t lca_vz = informationReferenceTree->lowestCommonAncestorIdx(vIdx, zIdx, rootIdx);
					size_t lca_wz = informationReferenceTree->lowestCommonAncestorIdx(wIdx, zIdx, rootIdx);

					size_t aIdx, bIdx, cIdx, dIdx;

					if (informationReferenceTree->distanceInEdges(lca_uv, lca_wz)
							> informationReferenceTree->distanceInEdges(lca_uw, lca_vz)
							&& informationReferenceTree->distanceInEdges(lca_uv, lca_wz)
									> informationReferenceTree->distanceInEdges(lca_uz, lca_vw)) {
						aIdx = uIdx;
						bIdx = vIdx;
						cIdx = wIdx;
						dIdx = zIdx; // ab|cd = uv|wz
					} else if (informationReferenceTree->distanceInEdges(lca_uw, lca_vz)
							> informationReferenceTree->distanceInEdges(lca_uv, lca_wz)
							&& informationReferenceTree->distanceInEdges(lca_uw, lca_vz)
									> informationReferenceTree->distanceInEdges(lca_uz, lca_vw)) {
						aIdx = uIdx;
						bIdx = wIdx;
						cIdx = vIdx;
						dIdx = zIdx; // ab|cd = uw|vz
					} else if (informationReferenceTree->distanceInEdges(lca_uz, lca_vw)
							> informationReferenceTree->distanceInEdges(lca_uv, lca_wz)
							&& informationReferenceTree->distanceInEdges(lca_uz, lca_vw)
									> informationReferenceTree->distanceInEdges(lca_uw, lca_vz)) {
						aIdx = uIdx;
						bIdx = zIdx;
						cIdx = vIdx;
						dIdx = wIdx; // ab|cd = uz|vw
					} else {
						// else, we have a multifurcation and the quartet has none of these three topologies. In this case, ignore the quartet.
						continue;
					}

					std::tuple<size_t, size_t, size_t> quartetOccurrences = countQuartetOccurrences(aIdx, bIdx, cIdx,
							dIdx);

					{ //***** Code for QP-IC and EQP-IC scores start
						std::pair<size_t, size_t> nodePair = nodePairForQuartet(aIdx, bIdx, cIdx, dIdx);
						std::pair<size_t, size_t> nodePairSorted = std::pair<size_t, size_t>(
								std::min(nodePair.first, nodePair.second), std::max(nodePair.first, nodePair.second));

						if (countBuffer.find(nodePairSorted) == countBuffer.end()) {
							std::tuple<size_t, size_t, size_t> emptyTuple(0, 0, 0);
							countBuffer[nodePairSorted] = emptyTuple;
						}

						size_t lcaIdx = informationReferenceTree->lowestCommonAncestorIdx(nodePairSorted.first,
								nodePairSorted.second, rootIdx);
						std::pair<size_t, size_t> innerLinks = get_path_inner_links(
								referenceTree.node_at(nodePairSorted.first),
								referenceTree.node_at(nodePairSorted.second), referenceTree.node_at(lcaIdx));
						size_t linkSubtree1 = referenceTree.link_at(innerLinks.first).next().index();
						size_t linkSubtree3 = referenceTree.link_at(innerLinks.second).next().index();

						size_t aLinkIdx = get_path_inner_links(referenceTree.node_at(nodePair.first),
								referenceTree.node_at(aIdx),
								referenceTree.node_at(
										informationReferenceTree->lowestCommonAncestorIdx(nodePair.first, aIdx,
												rootIdx))).first;
						size_t cLinkIdx = get_path_inner_links(referenceTree.node_at(nodePair.second),
								referenceTree.node_at(cIdx),
								referenceTree.node_at(
										informationReferenceTree->lowestCommonAncestorIdx(nodePair.second, cIdx,
												rootIdx))).first;

						size_t count_S1_S2_S3_S4;
						size_t count_S1_S3_S2_S4;
						size_t count_S1_S4_S2_S3;
						if ((aLinkIdx == linkSubtree1 && cLinkIdx == linkSubtree3)
								|| (cLinkIdx == linkSubtree1 && aLinkIdx == linkSubtree3)) {
							count_S1_S2_S3_S4 = std::get<0>(quartetOccurrences);
							count_S1_S3_S2_S4 = std::get<1>(quartetOccurrences);
							count_S1_S4_S2_S3 = std::get<2>(quartetOccurrences);
						} else {
							count_S1_S2_S3_S4 = std::get<0>(quartetOccurrences);
							count_S1_S3_S2_S4 = std::get<2>(quartetOccurrences);
							count_S1_S4_S2_S3 = std::get<1>(quartetOccurrences);
						}

						std::get<0>(countBuffer[nodePairSorted]) += count_S1_S2_S3_S4;
						std::get<1>(countBuffer[nodePairSorted]) += count_S1_S3_S2_S4;
						std::get<2>(countBuffer[nodePairSorted]) += count_S1_S4_S2_S3;
					} //***** Code for QP-IC and EQP-IC scores end

					double qic = log_score(std::get<0>(quartetOccurrences), std::get<1>(quartetOccurrences),
							std::get<2>(quartetOccurrences));

					// find path ends
					size_t lca_ab = informationReferenceTree->lowestCommonAncestorIdx(aIdx, bIdx, rootIdx);
					size_t lca_cd = informationReferenceTree->lowestCommonAncestorIdx(cIdx, dIdx, rootIdx);
					size_t fromIdx, toIdx;
					if (lca_cd == informationReferenceTree->lowestCommonAncestorIdx(cIdx, dIdx, lca_ab)) {
						fromIdx = informationReferenceTree->lowestCommonAncestorIdx(aIdx, bIdx, lca_cd);
						toIdx = lca_cd;
					} else {
						fromIdx = lca_ab;
						toIdx = informationReferenceTree->lowestCommonAncestorIdx(cIdx, dIdx, lca_ab);
					}
					// update the LQ-IC scores of the edges from fromIdx to toIdx
					size_t lcaFromToIdx = informationReferenceTree->lowestCommonAncestorIdx(fromIdx, toIdx, rootIdx);
					for (auto it : path_set(referenceTree.node_at(fromIdx), referenceTree.node_at(toIdx),
							referenceTree.node_at(lcaFromToIdx))) {
						if (it.is_lca())
							continue;
#pragma omp critical
						LQICScores[it.edge().index()] = std::min(LQICScores[it.edge().index()], qic);
					}
				}
			}
		}
	}

	{ // ***** Code for QP-IC and EQP-IC scores, finalizing, start
		for (auto kv : countBuffer) {
			std::pair<size_t, size_t> nodePair = kv.first;
			size_t uIdx = nodePair.first;
			size_t vIdx = nodePair.second;
			std::tuple<size_t, size_t, size_t> counts = kv.second;

			// compute the QP-IC score of the current metaquartet
			double qpic = log_score(std::get<0>(counts), std::get<1>(counts), std::get<2>(counts));

			// check if uIdx and vIdx are neighbors; if so, set QP-IC score of the edge connecting u and v
			auto const& u_link = referenceTree.node_at(uIdx).link();
			auto const& v_link = referenceTree.node_at(vIdx).link();
			if (u_link.outer().node().index() == vIdx) {
				QPICScores[u_link.edge().index()] = qpic;
			} else if (v_link.outer().node().index() == uIdx) {
				QPICScores[v_link.edge().index()] = qpic;
			}

			size_t lcaIdx = informationReferenceTree->lowestCommonAncestorIdx(uIdx, vIdx, rootIdx);
			// update the EQP-IC scores of the edges from uIdx to vIdx
			for (auto it : path_set(referenceTree.node_at(uIdx), referenceTree.node_at(vIdx),
					referenceTree.node_at(lcaIdx))) {
				if (it.is_lca())
					continue;
#pragma omp critical
				EQPICScores[it.edge().index()] = std::min(EQPICScores[it.edge().index()], qpic);
			}
		}
	} // ***** Code for QP-IC and EQP-IC scores, finalizing, end
}

/**
 * Update LQ-IC, QP-IC, and EQP-IC scores of all internodes influenced by the pair of inner nodes {u,v}.
 * As each quartet contributes to only one metaquartet induced by a node pair, our code visits each quartet exactly once.
 * @param uIdx the ID of the inner node u
 * @param vIdx the ID of the inner node v
 */
void QuartetScoreComputer::processNodePair(size_t uIdx, size_t vIdx) {
	// occurrences of topologies of the the metaquartet induced by {u,v} in the evaluation trees
	unsigned p1, p2, p3;
	p1 = 0;
	p2 = 0;
	p3 = 0;
	// find metaquartet indices by {u,v}

	size_t lcaIdx = informationReferenceTree->lowestCommonAncestorIdx(uIdx, vIdx, rootIdx);

	std::pair<size_t, size_t> innerLinks = get_path_inner_links(referenceTree.node_at(uIdx),
			referenceTree.node_at(vIdx), referenceTree.node_at(lcaIdx));

	size_t linkSubtree1 = referenceTree.link_at(innerLinks.first).next().index();
	size_t linkSubtree2 = referenceTree.link_at(innerLinks.first).next().next().index();
	size_t linkSubtree3 = referenceTree.link_at(innerLinks.second).next().index();
	size_t linkSubtree4 = referenceTree.link_at(innerLinks.second).next().next().index();
	// iterate over all quartets from {subtree1} x {subtree2} x {subtree3} x {subtree4}.
	// These are the quartets relevant to the current metaquartet.
	size_t startLeafIndexS1 = linkToEulerLeafIndex[linkSubtree1] % eulerTourLeaves.size();
	size_t endLeafIndexS1 = linkToEulerLeafIndex[referenceTree.link_at(linkSubtree1).outer().index()]
			% eulerTourLeaves.size();
	size_t startLeafIndexS2 = linkToEulerLeafIndex[linkSubtree2] % eulerTourLeaves.size();
	size_t endLeafIndexS2 = linkToEulerLeafIndex[referenceTree.link_at(linkSubtree2).outer().index()]
			% eulerTourLeaves.size();
	size_t startLeafIndexS3 = linkToEulerLeafIndex[linkSubtree3] % eulerTourLeaves.size();
	size_t endLeafIndexS3 = linkToEulerLeafIndex[referenceTree.link_at(linkSubtree3).outer().index()]
			% eulerTourLeaves.size();
	size_t startLeafIndexS4 = linkToEulerLeafIndex[linkSubtree4] % eulerTourLeaves.size();
	size_t endLeafIndexS4 = linkToEulerLeafIndex[referenceTree.link_at(linkSubtree4).outer().index()]
			% eulerTourLeaves.size();

	size_t aLeafIndex = startLeafIndexS1;
	size_t bLeafIndex = startLeafIndexS2;
	size_t cLeafIndex = startLeafIndexS3;
	size_t dLeafIndex = startLeafIndexS4;

	while (aLeafIndex != endLeafIndexS1) {
		while (bLeafIndex != endLeafIndexS2) {
			while (cLeafIndex != endLeafIndexS3) {
				while (dLeafIndex != endLeafIndexS4) {
					size_t aIdx = eulerTourLeaves[aLeafIndex];
					size_t bIdx = eulerTourLeaves[bLeafIndex];
					size_t cIdx = eulerTourLeaves[cLeafIndex];
					size_t dIdx = eulerTourLeaves[dLeafIndex];

					// process the quartet (a,b,c,d)
					// We already know by the way we defined S1,S2,S3,S4 that the reference tree has the quartet topology ab|cd
					std::tuple<size_t, size_t, size_t> quartetOccurrences = countQuartetOccurrences(aIdx, bIdx, cIdx,
							dIdx);
					p1 += std::get<0>(quartetOccurrences);
					p2 += std::get<1>(quartetOccurrences);
					p3 += std::get<2>(quartetOccurrences);
					double qic = log_score(std::get<0>(quartetOccurrences), std::get<1>(quartetOccurrences),
							std::get<2>(quartetOccurrences));

					// find path ends
					size_t lca_ab = informationReferenceTree->lowestCommonAncestorIdx(aIdx, bIdx, rootIdx);
					size_t lca_cd = informationReferenceTree->lowestCommonAncestorIdx(cIdx, dIdx, rootIdx);
					size_t fromIdx, toIdx;
					if (lca_cd == informationReferenceTree->lowestCommonAncestorIdx(cIdx, dIdx, lca_ab)) {
						fromIdx = informationReferenceTree->lowestCommonAncestorIdx(aIdx, bIdx, lca_cd);
						toIdx = lca_cd;
					} else {
						fromIdx = lca_ab;
						toIdx = informationReferenceTree->lowestCommonAncestorIdx(cIdx, dIdx, lca_ab);
					}
					// update the LQ-IC scores of the edges from fromIdx to toIdx
					size_t lcaFromToIdx = informationReferenceTree->lowestCommonAncestorIdx(fromIdx, toIdx, rootIdx);
					for (auto it : path_set(referenceTree.node_at(fromIdx), referenceTree.node_at(toIdx),
							referenceTree.node_at(lcaFromToIdx))) {
						if (it.is_lca())
							continue;
#pragma omp critical
						LQICScores[it.edge().index()] = std::min(LQICScores[it.edge().index()], qic);
					}

					dLeafIndex = (dLeafIndex + 1) % eulerTourLeaves.size();
				}
				cLeafIndex = (cLeafIndex + 1) % eulerTourLeaves.size();
				dLeafIndex = startLeafIndexS4;
			}
			bLeafIndex = (bLeafIndex + 1) % eulerTourLeaves.size();
			cLeafIndex = startLeafIndexS3;
			dLeafIndex = startLeafIndexS4;
		}
		aLeafIndex = (aLeafIndex + 1) % eulerTourLeaves.size();
		bLeafIndex = startLeafIndexS2;
		cLeafIndex = startLeafIndexS3;
		dLeafIndex = startLeafIndexS4;
	}

	// compute the QP-IC score of the current metaquartet
	double qpic = log_score(p1, p2, p3);

	// check if uIdx and vIdx are neighbors; if so, set QP-IC score of the edge connecting u and v
	auto const& u_link = referenceTree.node_at(uIdx).link();
	auto const& v_link = referenceTree.node_at(vIdx).link();
	if (u_link.outer().node().index() == vIdx) {
		QPICScores[u_link.edge().index()] = qpic;
	} else if (v_link.outer().node().index() == uIdx) {
		QPICScores[v_link.edge().index()] = qpic;
	}

	// update the EQP-IC scores of the edges from uIdx to vIdx
	for (auto it : path_set(referenceTree.node_at(uIdx), referenceTree.node_at(vIdx), referenceTree.node_at(lcaIdx))) {
		if (it.is_lca())
			continue;
#pragma omp critical
		EQPICScores[it.edge().index()] = std::min(EQPICScores[it.edge().index()], qpic);
	}
}

/**
 * Compute the LQ-IC, QP-IC, and EQP-IC support scores, iterating over node pairs.
 */
void QuartetScoreComputer::computeQuartetScoresBifurcating() {
	// Process all pairs of inner nodes
#pragma omp parallel for schedule(dynamic)
	for (size_t i = 0; i < referenceTree.node_count(); ++i) {
		if (!referenceTree.node_at(i).is_inner())
			continue;
		for (size_t j = i + 1; j < referenceTree.node_count(); ++j) {
			if (!referenceTree.node_at(j).is_inner())
				continue;
			processNodePair(i, j);
		}
	}
}

/**
 * Compute LQ-IC scores for a (possibly multifurcating) reference tree. Iterate over quartets.
 */
void QuartetScoreComputer::computeQuartetScoresMultifurcating() {
	// Process all quartets
#pragma omp parallel for schedule(dynamic)
	for (size_t uLeafIdx = 0; uLeafIdx < eulerTourLeaves.size(); uLeafIdx++) {
		for (size_t vLeafIdx = uLeafIdx + 1; vLeafIdx < eulerTourLeaves.size(); vLeafIdx++) {
			for (size_t wLeafIdx = vLeafIdx + 1; wLeafIdx < eulerTourLeaves.size(); wLeafIdx++) {
				for (size_t zLeafIdx = wLeafIdx + 1; zLeafIdx < eulerTourLeaves.size(); zLeafIdx++) {
					size_t uIdx = eulerTourLeaves[uLeafIdx];
					size_t vIdx = eulerTourLeaves[vLeafIdx];
					size_t wIdx = eulerTourLeaves[wLeafIdx];
					size_t zIdx = eulerTourLeaves[zLeafIdx];
					// find topology ab|cd of {u,v,w,z}
					size_t lca_uv = informationReferenceTree->lowestCommonAncestorIdx(uIdx, vIdx, rootIdx);
					size_t lca_uw = informationReferenceTree->lowestCommonAncestorIdx(uIdx, wIdx, rootIdx);
					size_t lca_uz = informationReferenceTree->lowestCommonAncestorIdx(uIdx, zIdx, rootIdx);
					size_t lca_vw = informationReferenceTree->lowestCommonAncestorIdx(vIdx, wIdx, rootIdx);
					size_t lca_vz = informationReferenceTree->lowestCommonAncestorIdx(vIdx, zIdx, rootIdx);
					size_t lca_wz = informationReferenceTree->lowestCommonAncestorIdx(wIdx, zIdx, rootIdx);

					size_t aIdx, bIdx, cIdx, dIdx;

					if (informationReferenceTree->distanceInEdges(lca_uv, lca_wz)
							> informationReferenceTree->distanceInEdges(lca_uw, lca_vz)
							&& informationReferenceTree->distanceInEdges(lca_uv, lca_wz)
									> informationReferenceTree->distanceInEdges(lca_uz, lca_vw)) {
						aIdx = uIdx;
						bIdx = vIdx;
						cIdx = wIdx;
						dIdx = zIdx; // ab|cd = uv|wz
					} else if (informationReferenceTree->distanceInEdges(lca_uw, lca_vz)
							> informationReferenceTree->distanceInEdges(lca_uv, lca_wz)
							&& informationReferenceTree->distanceInEdges(lca_uw, lca_vz)
									> informationReferenceTree->distanceInEdges(lca_uz, lca_vw)) {
						aIdx = uIdx;
						bIdx = wIdx;
						cIdx = vIdx;
						dIdx = zIdx; // ab|cd = uw|vz
					} else if (informationReferenceTree->distanceInEdges(lca_uz, lca_vw)
							> informationReferenceTree->distanceInEdges(lca_uv, lca_wz)
							&& informationReferenceTree->distanceInEdges(lca_uz, lca_vw)
									> informationReferenceTree->distanceInEdges(lca_uw, lca_vz)) {
						aIdx = uIdx;
						bIdx = zIdx;
						cIdx = vIdx;
						dIdx = wIdx; // ab|cd = uz|vw
					} else {
						// else, we have a multifurcation and the quartet has none of these three topologies. In this case, ignore the quartet.
						continue;
					}

					std::tuple<size_t, size_t, size_t> quartetOccurrences = countQuartetOccurrences(aIdx, bIdx, cIdx,
							dIdx);

					double qic = log_score(std::get<0>(quartetOccurrences), std::get<1>(quartetOccurrences),
							std::get<2>(quartetOccurrences));

					// find path ends
					size_t lca_ab = informationReferenceTree->lowestCommonAncestorIdx(aIdx, bIdx, rootIdx);
					size_t lca_cd = informationReferenceTree->lowestCommonAncestorIdx(cIdx, dIdx, rootIdx);
					size_t fromIdx, toIdx;
					if (lca_cd == informationReferenceTree->lowestCommonAncestorIdx(cIdx, dIdx, lca_ab)) {
						fromIdx = informationReferenceTree->lowestCommonAncestorIdx(aIdx, bIdx, lca_cd);
						toIdx = lca_cd;
					} else {
						fromIdx = lca_ab;
						toIdx = informationReferenceTree->lowestCommonAncestorIdx(cIdx, dIdx, lca_ab);
					}
					// update the LQ-IC scores of the edges from fromIdx to toIdx
					size_t lcaFromToIdx = informationReferenceTree->lowestCommonAncestorIdx(fromIdx, toIdx, rootIdx);
					for (auto it : path_set(referenceTree.node_at(fromIdx), referenceTree.node_at(toIdx),
							referenceTree.node_at(lcaFromToIdx))) {
						if (it.is_lca())
							continue;
#pragma omp critical
						LQICScores[it.edge().index()] = std::min(LQICScores[it.edge().index()], qic);
					}
				}
			}
		}
	}
}

/**
 * Get the start and stop indices of the leaves belonging to the subtree induced by the given link.
 * The leaf indices are between [start,stop).
 * @param linkIdx the genesis ID of the link inducing a subtree
 */
std::pair<size_t, size_t> QuartetScoreComputer::subtreeLeafIndices(size_t linkIdx) {
	size_t outerLinkIdx = referenceTree.link_at(linkIdx).outer().index();
	return {linkToEulerLeafIndex[linkIdx] % linkToEulerLeafIndex.size(), linkToEulerLeafIndex[outerLinkIdx] % linkToEulerLeafIndex.size()};
}

/**
 * Get the amount of total system RAM memory available.
 */
// Code adapted from http://stackoverflow.com/questions/2513505/how-to-get-available-memory-c-g
inline size_t getTotalSystemMemory() {
	size_t pages = sysconf(_SC_PHYS_PAGES);
	size_t page_size = sysconf(_SC_PAGE_SIZE);
	return pages * page_size;
}

/**
 * Count the number of evaluation trees.
 * @param evalTreesPath path to the file containing the set of evaluation trees
 */
size_t QuartetScoreComputer::countEvalTrees(const std::string &evalTreesPath) {
	size_t count = 0;
	utils::InputStream instream(utils::make_unique<utils::FileInputSource>(evalTreesPath));
	auto it = NewickInputIterator(instream);
	while (it) {
		count++;
		++it;
	}
	return count;
}

/**
 * @param refTree the reference tree
 * @param evalTrees path to the file containing the set of evaluation trees
 * @param verboseOutput print some additional (debug) information
 */
QuartetScoreComputer::QuartetScoreComputer(Tree const &refTree, const std::string &evalTreesPath, bool verboseOutput) {
	referenceTree = refTree;
	rootIdx = referenceTree.root_node().index();

	verbose = verboseOutput;
	size_t m = countEvalTrees(evalTreesPath);

	std::cout << "There are " << m << " evaluation trees.\n";

	std::cout << "Building subtree informations for reference tree..." << std::endl;
	informationReferenceTree = make_unique<TreeInformation>(referenceTree);
	// precompute subtree informations
	linkToEulerLeafIndex.resize(referenceTree.link_count());
	for (auto it : eulertour(referenceTree)) {
		if (it.node().is_leaf()) {
			eulerTourLeaves.push_back(it.node().index());
		}
		linkToEulerLeafIndex[it.link().index()] = eulerTourLeaves.size();
	}
	size_t n = eulerTourLeaves.size();

	std::cout << "Finished precomputing subtree informations in reference tree.\n";
	std::cout << "The reference tree has " << n << " taxa.\n";

	//estimate memory requirements for lookup table
	size_t memoryLookup = n * n * n * n * sizeof(cint);
	//estimate memory requirements for few trees
	size_t memoryFewTrees = m * n * 4 * sizeof(size_t);
	size_t estimatedMemory = getTotalSystemMemory();

	std::cout << "Estimated memory usages:" << std::endl;
	std::cout << "  Lookup table: " << memoryLookup << std::endl;
	std::cout << "  Distance and LCA approach: " << memoryFewTrees << std::endl;
	std::cout << "  Estimated available memory: " << estimatedMemory << std::endl;

	std::chrono::steady_clock::time_point begin = std::chrono::steady_clock::now();

	if (memoryLookup < estimatedMemory / 1.5) {
		std::cout << "Using lookup table.\n";
		quartetCounterLookup = make_unique<QuartetCounterLookup>(refTree, evalTreesPath, m);
		//quartetCounterLookup = make_unique<QuartetCounterLookupCacheOptimized>(refTree, evalTrees);
		useLookupTable = true;
	} else {
		std::cout << "Using distance and LCA approach.\n";
		quartetCounterFewTrees = make_unique<QuartetCounterFewTrees>(refTree, evalTreesPath);
		useLookupTable = false;
	}

	std::chrono::steady_clock::time_point end = std::chrono::steady_clock::now();

	std::cout << "Finished counting quartets.\n";
	std::cout << "It took: " << std::chrono::duration_cast<std::chrono::microseconds>(end - begin).count()
			<< " microseconds." << std::endl;

	begin = std::chrono::steady_clock::now();

	// precompute taxon ID mappings
	// this is commented out as it was not needed anywhere in the code.
	//taxonMapper = make_unique<TaxonMapper>(referenceTree, evaluationTrees, eulerTourLeaves);
	//std::cout << "Finished precomputing taxon ID mappings.\n";

	if (!is_bifurcating(refTree)) {
		std::cout << "The reference tree is multifurcating.\n";
		LQICScores.resize(referenceTree.edge_count());
		std::fill(LQICScores.begin(), LQICScores.end(), std::numeric_limits<double>::infinity());
		// compute only LQ-IC scores
		computeQuartetScoresMultifurcating();
	} else {
		std::cout << "The reference tree is bifurcating.\n";
		LQICScores.resize(referenceTree.edge_count());
		QPICScores.resize(referenceTree.edge_count());
		EQPICScores.resize(referenceTree.edge_count());
		// initialize the LQ-IC, QP-IC and EQP-IC scores of all internodes (edges) to INFINITY
		std::fill(LQICScores.begin(), LQICScores.end(), std::numeric_limits<double>::infinity());
		std::fill(QPICScores.begin(), QPICScores.end(), std::numeric_limits<double>::infinity());
		std::fill(EQPICScores.begin(), EQPICScores.end(), std::numeric_limits<double>::infinity());
		// compute LQ-IC, QP-IC and EQP-IC scores
		computeQuartetScoresBifurcating();
		//computeQuartetScoresBifurcatingQuartets();
	}

	end = std::chrono::steady_clock::now();

	std::cout << "Finished computing scores.\n";
	std::cout << "It took: " << std::chrono::duration_cast<std::chrono::microseconds>(end - begin).count()
			<< " microseconds." << std::endl;
}
