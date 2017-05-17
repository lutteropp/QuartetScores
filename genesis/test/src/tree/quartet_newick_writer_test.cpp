/*
    Genesis - A toolkit for working with phylogenetic data.
    Copyright (C) 2014-2016 Lucas Czech

    This program is free software: you can redistribute it and/or modify
    it under the terms of the GNU General Public License as published by
    the Free Software Foundation, either version 3 of the License, or
    (at your option) any later version.

    This program is distributed in the hope that it will be useful,
    but WITHOUT ANY WARRANTY; without even the implied warranty of
    MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
    GNU General Public License for more details.

    You should have received a copy of the GNU General Public License
    along with this program.  If not, see <http://www.gnu.org/licenses/>.

    Contact:
    Lucas Czech <lucas.czech@h-its.org>
    Exelixis Lab, Heidelberg Institute for Theoretical Studies
    Schloss-Wolfsbrunnenweg 35, D-69118 Heidelberg, Germany
*/

/**
 * @brief Testing Newick class.
 *
 * @file
 * @ingroup test
 */

#include "common.hpp"

// Path to the quartet writer header.
#include "lib/tree/formats/quartet_newick_writer.hpp"

#include "lib/tree/default/newick_reader.hpp"
#include "lib/tree/default/newick_writer.hpp"
#include "lib/tree/formats/newick/reader.hpp"
#include "lib/tree/formats/newick/writer.hpp"
#include "lib/tree/tree.hpp"

using namespace genesis;
using namespace tree;

TEST( Tree, QuartetNewickWriter )
{
    // Create a test tree.
    std::string input = "((A,(B,C)D)E,((F,(G,H)I)J,K)L)R;";
    Tree tree;
    EXPECT_TRUE( DefaultTreeNewickReader().from_string( input, tree ));

    // Create some test vectors, simply using the index as value.
    auto qp_vector  = std::vector<double>( tree.edge_count() );
    auto eqp_vector = std::vector<double>( tree.edge_count() );
    for( size_t i = 0; i < tree.edge_count(); ++i ) {
        qp_vector[i]  = static_cast<double>( i );
        eqp_vector[i] = static_cast<double>( i );
    }

    // Create the writer and assign values.
    auto writer = QuartetNewickWriter();
    writer.set_qp_ic_scores(  qp_vector  );
    writer.set_eqp_ic_scores( eqp_vector );

    // For testing.
    std::string output = writer.to_string( tree );
    LOG_DBG << output;

    // For production.
    // writer.to_file( tree, filename );
}
