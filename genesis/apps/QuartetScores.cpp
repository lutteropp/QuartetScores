#include "genesis/genesis.hpp"
#include "QuartetScoreComputer.hpp"
#include <string>
#include <vector>
#include <iostream>
#include <fstream>

#include <chrono>

using namespace genesis;
using namespace tree;

/**
 * Count the number of evaluation trees.
 * @param evalTreesPath path to the file containing the set of evaluation trees
 */
size_t countEvalTrees(const std::string &evalTreesPath) {
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
 * The main method. Compute quartet scores and store the result in a tree file.
 */
int main(int argc, char* argv[]) {
	std::chrono::steady_clock::time_point begin = std::chrono::steady_clock::now();

	// TODO: Change argument list to also support --fewtrees, --lookup, and --auto (default).
	// TODO: Something like --nthreads would also be nice to have.

	if (argc < 4 || argc > 5) {
		std::cout << "Usage: " << argv[0] << " <Optional:-v> <refTreeFilePath> <evalTreesFilePath> <outputFilePath>"
				<< std::endl;
		return 1;
	}

	bool verbose = false;
	if (argc == 5)
		verbose = true;

	std::string pathToReferenceTree;
	std::string pathToEvaluationTrees;
	std::string outputFilePath;

	if (argc == 4) {
		pathToReferenceTree = std::string(argv[1]);
		pathToEvaluationTrees = std::string(argv[2]);
		outputFilePath = std::string(argv[3]);
	} else {
		pathToReferenceTree = std::string(argv[2]);
		pathToEvaluationTrees = std::string(argv[3]);
		outputFilePath = std::string(argv[4]);
	}

	std::ifstream infile(outputFilePath);
	if (infile.good()) {
		std::cout << "ERROR: The specified output file already exists.\n";
		return 1;
	}

	//read trees
	DefaultTreeNewickReader reader;
	Tree referenceTree = reader.from_file(pathToReferenceTree);

	/*for( auto const& tree : tree_iter(...)) {

	 }*/

	if (verbose) {
		auto tp = PrinterCompact();
		auto res = tp.print(referenceTree, []( TreeNode const& node, TreeEdge const& edge ) {
			//return node.data<DefaultNodeData>().name + " edge: " + std::to_string(edge.index()) + " vertex: " + std::to_string(node.index());
				return node.data<DefaultNodeData>().name + " " + std::to_string( node.index() );
			});

		std::cout << res << std::endl;
	}

	std::vector<double> lqic;
	std::vector<double> qpic;
	std::vector<double> eqpic;
	size_t m = countEvalTrees(pathToEvaluationTrees);
	if (m < (size_t(1) << 8)) {
		QuartetScoreComputer<uint8_t> qsc(referenceTree, pathToEvaluationTrees, m, verbose);
		lqic = qsc.getLQICScores();
		qpic = qsc.getQPICScores();
		eqpic = qsc.getEQPICScores();
	} else if (m < (size_t(1) << 16)) {
		QuartetScoreComputer<uint16_t> qsc(referenceTree, pathToEvaluationTrees, m, verbose);
		lqic = qsc.getLQICScores();
		qpic = qsc.getQPICScores();
		eqpic = qsc.getEQPICScores();
	} else if (m < (size_t(1) << 32)) {
		QuartetScoreComputer<uint32_t> qsc(referenceTree, pathToEvaluationTrees, m, verbose);
		lqic = qsc.getLQICScores();
		qpic = qsc.getQPICScores();
		eqpic = qsc.getEQPICScores();
	} else {
		QuartetScoreComputer<uint64_t> qsc(referenceTree, pathToEvaluationTrees, m, verbose);
		lqic = qsc.getLQICScores();
		qpic = qsc.getQPICScores();
		eqpic = qsc.getEQPICScores();
	}
	// Create the writer and assign values.
	auto writer = QuartetTreeNewickWriter();
	writer.set_lq_ic_scores(lqic);
	if (!eqpic.empty()) { // bifurcating tree
		writer.set_eqp_ic_scores(eqpic);
	}
	if (!qpic.empty()) { // bifurcating tree
		writer.set_qp_ic_scores(qpic);
	}

	writer.to_file(referenceTree, outputFilePath);

	std::chrono::steady_clock::time_point end = std::chrono::steady_clock::now();

	std::cout << "Elapsed time: " << std::chrono::duration_cast<std::chrono::microseconds>(end - begin).count()
			<< " microseconds." << std::endl;
	return 0;
}
