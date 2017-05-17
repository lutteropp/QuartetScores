#include "genesis.hpp"
#include "QuartetScoreComputerNaive.hpp"
#include <string>
#include <vector>
#include <iostream>
#include <fstream>

#include "tree/formats/quartet_newick_writer.hpp"

using namespace genesis;
using namespace tree;

int main(int argc, char* argv[]) {
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

				return node.data<DefaultNodeData>().name + " " + std::to_string( edge.index() );
			});

		std::cout << res << std::endl;
	}

	QuartetScoreComputerNaive qsc(referenceTree, evaluationTrees, verbose);
	qsc.computeLQICScores();
	qsc.computeQPICScores();
	qsc.computeEQPICScores();
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
	return 0;
}
