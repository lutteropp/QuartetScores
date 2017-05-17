#pragma once

#include "genesis.hpp"
#include <vector>
#include <cassert>
#include <algorithm>
#include <memory>
#include "TreeInformation.hpp"
#include <unordered_map>
#include "QuartetCounterLookup.hpp"
#include "TaxonMapper.hpp"

using namespace genesis;
using namespace tree;
using namespace utils;

/**
 * Compute LQ-IC, QP-IC, and EQP-IC scores by directly applying the scores' definitions.
 * This is, quartets will be visited more than once as we compute the scores in an internode-based approach.
 */
class QuartetScoreComputerNaive {
public:
	QuartetScoreComputerNaive(Tree const &refTree, TreeSet const &evalTrees, bool printDebug);
	std::vector<double> getLQICScores();
	std::vector<double> getQPICScores();
	std::vector<double> getEQPICScores();
	void computeLQICScores();
	void computeQPICScores();
	void computeEQPICScores();
private:
	double log_score(size_t q1, size_t q2, size_t q3);
	void computeQuartetScoresBifurcating();
	void computeQuartetScoresMultifurcating();
	void processNodePair(size_t uIdx, size_t vIdx);
	std::tuple<size_t, size_t, size_t> countQuartetOccurrences(size_t aIdx, size_t bIdx, size_t cIdx, size_t dIdx);
	Tree referenceTree;
	TreeSet evaluationTrees;

	bool verbose;

	std::unique_ptr<TreeInformation> informationReferenceTree;
	std::vector<double> LQICScores;
	std::vector<double> QPICScores;
	std::vector<double> EQPICScores;

	std::vector<std::string> LQIC_Origins; // TODO: Just for debug
	std::vector<std::string> EQPIC_Origins; // TODO: Just for debug

	std::pair<size_t, size_t> subtreeLeafIndices(size_t linkIdx);
	std::vector<size_t> eulerTourLeaves;
	std::vector<size_t> linkToEulerLeafIndex;

	std::unique_ptr<TaxonMapper> taxonMapper; // TODO: this seems to be unnecessary
	std::unique_ptr<QuartetCounterLookup> quartetCounter;
};

/**
 * Return the computed LQ-IC support scores.
 */
std::vector<double> QuartetScoreComputerNaive::getLQICScores() {
	return LQICScores;
}

/**
 * Return the computed QP-IC support scores.
 */
std::vector<double> QuartetScoreComputerNaive::getQPICScores() {
	return QPICScores;
}

/**
 * Return the computed EQP-IC support scores.
 */
std::vector<double> QuartetScoreComputerNaive::getEQPICScores() {
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
double QuartetScoreComputerNaive::log_score(size_t q1, size_t q2, size_t q3) {
	if (q1 == 0 && q2 == 0 && q3 == 0)
		return 0;

	size_t sum = q1 + q2 + q3;
	double p_q1 = (double) q1 / sum;
	double p_q2 = (double) q2 / sum;
	double p_q3 = (double) q3 / sum;

	double qic = 1;
	if (p_q1 != 0)
		qic += p_q1 * log(p_q1) / log(3);
	if (p_q2 != 0)
		qic += p_q2 * log(p_q2) / log(3);
	if (p_q3 != 0)
		qic += p_q3 * log(p_q3) / log(3);

	//double qic = 1 + p_q1 * log(p_q1) + p_q2 * log(p_q2) + p_q3 * log(p_q3);

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
std::tuple<size_t, size_t, size_t> QuartetScoreComputerNaive::countQuartetOccurrences(size_t aIdx, size_t bIdx,
		size_t cIdx, size_t dIdx) {
	return quartetCounter->countQuartetOccurrences(aIdx, bIdx, cIdx, dIdx);
}

/**
 * Compute the LQ-IC scores by iterating over the internodes (this is, naively applying the definition). Quartets are visited more than once.
 */
void QuartetScoreComputerNaive::computeLQICScores() {
#pragma omp parallel for
	for (size_t i = 0; i < referenceTree.edge_count(); ++i) { // for each edge
		// find left and right subtree, i.e. bipartition indexed with i
		std::pair<size_t, size_t> subtreeLeaves = subtreeLeafIndices(referenceTree.edge_at(i).primary_link().index());
		std::vector<size_t> nodeIndicesLeft;
		std::vector<size_t> nodeIndicesRight;
		size_t j = subtreeLeaves.first % eulerTourLeaves.size();
		nodeIndicesLeft.push_back(referenceTree.node_at(eulerTourLeaves[j]).index());
		j = (j + 1) % eulerTourLeaves.size();
		while (j != subtreeLeaves.second % eulerTourLeaves.size()) {
			nodeIndicesLeft.push_back(referenceTree.node_at(eulerTourLeaves[j]).index());
			j = (j + 1) % eulerTourLeaves.size();
		}
		nodeIndicesRight.push_back(referenceTree.node_at(eulerTourLeaves[j]).index());
		j = (j + 1) % eulerTourLeaves.size();
		while (j != subtreeLeaves.first % eulerTourLeaves.size()) {
			nodeIndicesRight.push_back(referenceTree.node_at(eulerTourLeaves[j]).index());
			j = (j + 1) % eulerTourLeaves.size();
		}

		for (size_t i1 = 0; i1 < nodeIndicesLeft.size(); ++i1) {
			size_t aIdx = nodeIndicesLeft[i1];
			for (size_t i2 = i1 + 1; i2 < nodeIndicesLeft.size(); ++i2) {
				size_t bIdx = nodeIndicesLeft[i2];
				for (size_t j1 = 0; j1 < nodeIndicesRight.size(); ++j1) {
					size_t cIdx = nodeIndicesRight[j1];
					for (size_t j2 = j1 + 1; j2 < nodeIndicesRight.size(); ++j2) {
						size_t dIdx = nodeIndicesRight[j2];
						std::tuple<size_t, size_t, size_t> quartetOccurrences = countQuartetOccurrences(aIdx, bIdx,
								cIdx, dIdx);
						double qic = log_score(std::get<0>(quartetOccurrences), std::get<1>(quartetOccurrences),
								std::get<2>(quartetOccurrences));
						LQICScores[i] = std::min(LQICScores[i], qic);
					}
				}
			}
		}
	}
}

/**
 * Compute the QP-IC scores by iterating over the internodes (this is, naively applying the definition). Quartets are visited more than once.
 */
void QuartetScoreComputerNaive::computeQPICScores() {
	size_t rootIdx = referenceTree.root_node().index();
#pragma omp parallel for
	for (size_t i = 0; i < referenceTree.edge_count(); ++i) { // for each edge
		size_t uIdx = referenceTree.edge_at(i).primary_node().index();
		size_t vIdx = referenceTree.edge_at(i).secondary_node().index();
		// occurrences of topologies of the the metaquartet induced by {u,v} in the evaluation trees
		unsigned p1 = 0;
		unsigned p2 = 0;
		unsigned p3 = 0;
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
						std::tuple<size_t, size_t, size_t> quartetOccurrences = countQuartetOccurrences(aIdx, bIdx,
								cIdx, dIdx);
						p1 += std::get<0>(quartetOccurrences);
						p2 += std::get<1>(quartetOccurrences);
						p3 += std::get<2>(quartetOccurrences);

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
		QPICScores[i] = qpic;
	}
}

/**
 * Compute the EQP-IC scores by iterating over the internodes (this is, naively applying the definition). Quartets are visited more than once.
 */
void QuartetScoreComputerNaive::computeEQPICScores() {
	// TODO: FIXME: Implement me.
	throw std::runtime_error("not implemented yet!");
}

void QuartetScoreComputerNaive::processNodePair(size_t uIdx, size_t vIdx) {
	size_t rootIdx = referenceTree.root_node().index();
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

	// TODO: This is just for debug/ verbose mode!
	std::string metaS1;
	std::string metaS2;
	std::string metaS3;
	std::string metaS4;
	if (verbose) {
		size_t aLeafIndexDebug = startLeafIndexS1;
		size_t bLeafIndexDebug = startLeafIndexS2;
		size_t cLeafIndexDebug = startLeafIndexS3;
		size_t dLeafIndexDebug = startLeafIndexS4;

		metaS1 = metaS1 + referenceTree.node_at(eulerTourLeaves[aLeafIndexDebug]).data<DefaultNodeData>().name;
		aLeafIndexDebug = (aLeafIndexDebug + 1) % eulerTourLeaves.size();
		while (aLeafIndexDebug != endLeafIndexS1) {
			metaS1 = metaS1 + ","
					+ referenceTree.node_at(eulerTourLeaves[aLeafIndexDebug]).data<DefaultNodeData>().name;
			aLeafIndexDebug = (aLeafIndexDebug + 1) % eulerTourLeaves.size();
		}
		metaS2 = metaS2 + referenceTree.node_at(eulerTourLeaves[bLeafIndexDebug]).data<DefaultNodeData>().name;
		bLeafIndexDebug = (bLeafIndexDebug + 1) % eulerTourLeaves.size();
		while (bLeafIndexDebug != endLeafIndexS2) {
			metaS2 = metaS2 + ","
					+ referenceTree.node_at(eulerTourLeaves[bLeafIndexDebug]).data<DefaultNodeData>().name;
			bLeafIndexDebug = (bLeafIndexDebug + 1) % eulerTourLeaves.size();
		}
		metaS3 = metaS3 + referenceTree.node_at(eulerTourLeaves[cLeafIndexDebug]).data<DefaultNodeData>().name;
		cLeafIndexDebug = (cLeafIndexDebug + 1) % eulerTourLeaves.size();
		while (cLeafIndexDebug != endLeafIndexS3) {
			metaS3 = metaS3 + ","
					+ referenceTree.node_at(eulerTourLeaves[cLeafIndexDebug]).data<DefaultNodeData>().name;
			cLeafIndexDebug = (cLeafIndexDebug + 1) % eulerTourLeaves.size();
		}
		metaS4 = metaS4 + referenceTree.node_at(eulerTourLeaves[dLeafIndexDebug]).data<DefaultNodeData>().name;
		dLeafIndexDebug = (dLeafIndexDebug + 1) % eulerTourLeaves.size();
		while (dLeafIndexDebug != endLeafIndexS4) {
			metaS4 = metaS4 + ","
					+ referenceTree.node_at(eulerTourLeaves[dLeafIndexDebug]).data<DefaultNodeData>().name;
			dLeafIndexDebug = (dLeafIndexDebug + 1) % eulerTourLeaves.size();
		}

		//#pragma omp critical
		//std::cout << "Node pair (" << uIdx << "," << vIdx << ") has metaquartet {" << metaS1 << "}, {" << metaS2 << "}, {" << metaS3 << "}, {" << metaS4 << "}\n";
	}

	//TODO: Remove again, this is just for debug
	size_t visitedQuartets = 0;

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

					//TODO: remove again.
					if (verbose) {
						visitedQuartets++;
					}

					if (verbose) {
#pragma omp critical
						{
							std::cout << "Quartet (" << referenceTree.node_at(aIdx).data<DefaultNodeData>().name << ","
									<< referenceTree.node_at(bIdx).data<DefaultNodeData>().name << "|"
									<< referenceTree.node_at(cIdx).data<DefaultNodeData>().name << ","
									<< referenceTree.node_at(dIdx).data<DefaultNodeData>().name << ") has occurrences: "
									<< std::get<0>(quartetOccurrences) << ", " << std::get<1>(quartetOccurrences)
									<< ", " << std::get<2>(quartetOccurrences) << std::endl;

							std::cout << "Quartet (" << referenceTree.node_at(aIdx).data<DefaultNodeData>().name << ","
									<< referenceTree.node_at(bIdx).data<DefaultNodeData>().name << "|"
									<< referenceTree.node_at(cIdx).data<DefaultNodeData>().name << ","
									<< referenceTree.node_at(dIdx).data<DefaultNodeData>().name << ") has QIC: " << qic
									<< std::endl;
						}
					}

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
						//std::cout << "Updating LQICScores[" << it.edge().index() << "] with " << qic << "\n";
#pragma omp critical
						{
							//LQICScores[it.edge().index()] = std::min(LQICScores[it.edge().index()], qic);
							if (qic < LQICScores[it.edge().index()]) {
								LQICScores[it.edge().index()] = qic;
								if (verbose) {
									std::string origin = "";
									origin = origin + referenceTree.node_at(aIdx).data<DefaultNodeData>().name + ","
											+ referenceTree.node_at(bIdx).data<DefaultNodeData>().name + "|"
											+ referenceTree.node_at(cIdx).data<DefaultNodeData>().name + ","
											+ referenceTree.node_at(dIdx).data<DefaultNodeData>().name;
									LQIC_Origins[it.edge().index()] = origin;
								}
							}
						}
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

	// TODO: This is just for debugging/ verbose mode
	std::string metaquartet;
	if (verbose) {
		metaquartet = "({" + metaS1 + "}, {" + metaS2 + "}, {" + metaS3 + "}, {" + metaS4 + "})";
#pragma omp critical
		std::cout << "Metaquartet " << metaquartet << " has counts p1=" << p1 << ",p2=" << p2 << ",p3=" << p3
				<< " and QPIC-score " << qpic << ". It consists of " << visitedQuartets << " quartets." << std::endl;
	}

	// update the EQP-IC scores of the edges from uIdx to vIdx
	for (auto it : path_set(referenceTree.node_at(uIdx), referenceTree.node_at(vIdx), referenceTree.node_at(lcaIdx))) {
		if (it.is_lca())
			continue;
#pragma omp critical
		{
			//EQPICScores[it.edge().index()] = std::min(EQPICScores[it.edge().index()], qpic);
			if (qpic < EQPICScores[it.edge().index()]) {
				EQPICScores[it.edge().index()] = qpic;
				if (verbose) {
					EQPIC_Origins[it.edge().index()] = metaquartet;
				}
			}
		}
	}
}

void QuartetScoreComputerNaive::computeQuartetScoresBifurcating() {
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

void QuartetScoreComputerNaive::computeQuartetScoresMultifurcating() {
	size_t rootIdx = referenceTree.root_node().index();

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

					if (verbose) {
#pragma omp critical
						{
							std::cout << "Quartet (" << referenceTree.node_at(aIdx).data<DefaultNodeData>().name << ","
									<< referenceTree.node_at(bIdx).data<DefaultNodeData>().name << "|"
									<< referenceTree.node_at(cIdx).data<DefaultNodeData>().name << ","
									<< referenceTree.node_at(dIdx).data<DefaultNodeData>().name << ") has occurrences: "
									<< std::get<0>(quartetOccurrences) << ", " << std::get<1>(quartetOccurrences)
									<< ", " << std::get<2>(quartetOccurrences) << std::endl;

							std::cout << "Quartet (" << referenceTree.node_at(aIdx).data<DefaultNodeData>().name << ","
									<< referenceTree.node_at(bIdx).data<DefaultNodeData>().name << "|"
									<< referenceTree.node_at(cIdx).data<DefaultNodeData>().name << ","
									<< referenceTree.node_at(dIdx).data<DefaultNodeData>().name << ") has QIC: " << qic
									<< std::endl;
						}
					}

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
						{
							//LQICScores[it.edge().index()] = std::min(LQICScores[it.edge().index()], qic);
							if (qic < LQICScores[it.edge().index()]) {
								LQICScores[it.edge().index()] = qic;
								if (verbose) {
									std::string origin = "";
									origin = origin + referenceTree.node_at(aIdx).data<DefaultNodeData>().name + ","
											+ referenceTree.node_at(bIdx).data<DefaultNodeData>().name + "|"
											+ referenceTree.node_at(cIdx).data<DefaultNodeData>().name + ","
											+ referenceTree.node_at(dIdx).data<DefaultNodeData>().name;
									LQIC_Origins[it.edge().index()] = origin;
								}
							}
						}
					}
				}
			}
		}
	}
}

// the leaf indizes are between [start,end)
std::pair<size_t, size_t> QuartetScoreComputerNaive::subtreeLeafIndices(size_t linkIdx) {
	size_t outerLinkIdx = referenceTree.link_at(linkIdx).outer().index();
	return {linkToEulerLeafIndex[linkIdx] % linkToEulerLeafIndex.size(), linkToEulerLeafIndex[outerLinkIdx] % linkToEulerLeafIndex.size()};
}

QuartetScoreComputerNaive::QuartetScoreComputerNaive(Tree const &refTree, TreeSet const &evalTrees, bool printDebug) {
	referenceTree = refTree;
	evaluationTrees = evalTrees;

	quartetCounter = make_unique<QuartetCounterLookup>(refTree, evalTrees);

	verbose = printDebug;

	std::cout << "There are " << evaluationTrees.size() << " evaluation trees.\n";

	informationReferenceTree = make_unique<TreeInformation>(referenceTree);

	std::cout << "Finished building tree informations.\n";

	// precompute subtree informations
	linkToEulerLeafIndex.resize(referenceTree.link_count());
	for (auto it : eulertour(referenceTree)) {
		if (it.node().is_leaf()) {
			eulerTourLeaves.push_back(it.node().index());
		}
		linkToEulerLeafIndex[it.link().index()] = eulerTourLeaves.size();
	}

	std::cout << "Finished precomputing subtree informations.\n";
	std::cout << "The reference tree has " << eulerTourLeaves.size() << " taxa.\n";

	// precompute taxon ID mappings
	// Commented out because it was not used anywhere in the code
	//taxonMapper = make_unique<TaxonMapper>(referenceTree, evaluationTrees, eulerTourLeaves);
	//std::cout << "Finished precomputing taxon ID mappings.\n";

	if (!is_bifurcating(refTree)) {
		std::cout << "The reference tree is multifurcating.\n";
		LQICScores.resize(referenceTree.edge_count());
		if (verbose) {
			LQIC_Origins.resize(referenceTree.edge_count());
		}
		std::fill(LQICScores.begin(), LQICScores.end(), std::numeric_limits<double>::infinity());
		// compute only LQ-IC scores
		//computeQuartetScoresMultifurcating();
	} else {
		std::cout << "The reference tree is bifurcating.\n";
		LQICScores.resize(referenceTree.edge_count());
		QPICScores.resize(referenceTree.edge_count());
		EQPICScores.resize(referenceTree.edge_count());
		if (verbose) {
			LQIC_Origins.resize(referenceTree.edge_count());
			EQPIC_Origins.resize(referenceTree.edge_count());
		}
		// initialize the LQ-IC, QP-IC and EQP-IC scores of all internodes (edges) to INFINITY
		std::fill(LQICScores.begin(), LQICScores.end(), std::numeric_limits<double>::infinity());
		std::fill(QPICScores.begin(), QPICScores.end(), std::numeric_limits<double>::infinity());
		std::fill(EQPICScores.begin(), EQPICScores.end(), std::numeric_limits<double>::infinity());
		// compute LQ-IC, QP-IC and EQP-IC scores
		//computeQuartetScoresBifurcating();
	}
	std::cout << "Finished computing scores.\n";
}
