#include "genesis.hpp"
#include "QuartetScoreComputer.hpp"
#include <string>
#include <vector>
#include <iostream>
#include <fstream>

#include <chrono>

#include "tree/formats/quartet_newick_writer.hpp"

using namespace genesis;
using namespace tree;

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
	Tree referenceTree;
	DefaultTreeNewickReader reader;
	reader.from_file(pathToReferenceTree, referenceTree);
	TreeSet evaluationTrees;
	reader.from_file(pathToEvaluationTrees, evaluationTrees);

	if (verbose) {
		auto tp = PrinterCompact();
		auto res = tp.print(referenceTree, []( TreeNode const& node, TreeEdge const& edge ) {
			//return node.data<DefaultNodeData>().name + " edge: " + std::to_string(edge.index()) + " vertex: " + std::to_string(node.index());
				return node.data<DefaultNodeData>().name + " " + std::to_string( node.index() );
			});

		std::cout << res << std::endl;
	}

	QuartetScoreComputer qsc(referenceTree, evaluationTrees, verbose);
	std::vector<double> lqic = qsc.getLQICScores();
	std::vector<double> qpic = qsc.getQPICScores();
	std::vector<double> eqpic = qsc.getEQPICScores();

	// Create the writer and assign values.
	auto writer = QuartetNewickWriter();
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
