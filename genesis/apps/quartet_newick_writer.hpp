#ifndef QUARTET_NEWICK_WRITER_H_
#define QUARTET_NEWICK_WRITER_H_

/*
    Genesis - A toolkit for working with phylogenetic data.
    Copyright (C) 2014-2017 Lucas Czech

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
 * @brief
 *
 * @file
 * @ingroup tree
 */

#include <assert.h>
#include <limits>
#include <stdexcept>
#include <string>
#include <vector>

#include "genesis/tree/default/newick_writer.hpp"
#include "genesis/tree/formats/newick/writer.hpp"
#include "genesis/utils/text/string.hpp"

namespace genesis {
namespace tree {

// =================================================================================================
//     Quartet Tree Newick Writer Mixin
// =================================================================================================

/**
 * @brief
 */
class QuartetNewickWriterPlugin
{
public:

    // -------------------------------------------------------------------------
    //     Typedefs and Enums
    // -------------------------------------------------------------------------

    using self_type = QuartetNewickWriterPlugin;

    // -------------------------------------------------------------------------
    //     Constructor and Rule of Five
    // -------------------------------------------------------------------------

    QuartetNewickWriterPlugin() = default;
    virtual ~QuartetNewickWriterPlugin() = default;

    QuartetNewickWriterPlugin(QuartetNewickWriterPlugin const&) = default;
    QuartetNewickWriterPlugin(QuartetNewickWriterPlugin&&)      = default;

    QuartetNewickWriterPlugin& operator= (QuartetNewickWriterPlugin const&) = default;
    QuartetNewickWriterPlugin& operator= (QuartetNewickWriterPlugin&&)      = default;

    // -------------------------------------------------------------------------
    //     Properties
    // -------------------------------------------------------------------------

    bool enable_lq_ic_scores() const
    {
        return enable_lq_ic_scores_;
    }

    bool enable_qp_ic_scores() const
    {
        return enable_qp_ic_scores_;
    }

    void enable_lq_ic_scores(bool value)
    {
        enable_lq_ic_scores_ = value;
    }

    void enable_qp_ic_scores(bool value)
    {
        enable_qp_ic_scores_ = value;
    }

    bool enable_eqp_ic_scores() const
    {
        return enable_eqp_ic_scores_;
    }

    void enable_eqp_ic_scores(bool value)
    {
        enable_eqp_ic_scores_ = value;
    }

    void set_lq_ic_scores( std::vector<double> const& values )
    {
        enable_lq_ic_scores_ = true;
        lq_ic_scores_ = values;
    }

    void set_qp_ic_scores( std::vector<double> const& values )
    {
        enable_qp_ic_scores_ = true;
        qp_ic_scores_ = values;
    }

    void set_eqp_ic_scores( std::vector<double> const& values )
    {
        enable_eqp_ic_scores_ = true;
        eqp_ic_scores_ = values;
    }


    // -------------------------------------------------------------------------
    //     Plugin Functions
    // -------------------------------------------------------------------------

    void prepare_writing( Tree const& tree, NewickBroker& broker ) const
    {
        (void) broker;
        if( enable_lq_ic_scores_ && tree.edge_count() != lq_ic_scores_.size() ) {
            throw std::runtime_error(
                "Size of QP-IC-Scores (" + std::to_string( lq_ic_scores_.size() ) +
                ") does not match numer of edges of the tree (" +
                std::to_string( tree.edge_count() ) + ")."
            );
        }
        if( enable_eqp_ic_scores_ && tree.edge_count() != eqp_ic_scores_.size() ) {
            throw std::runtime_error(
                "Size of EQP-IC-Scores (" + std::to_string( eqp_ic_scores_.size() ) +
                ") does not match numer of edges of the tree (" +
                std::to_string( tree.edge_count() ) + ")."
            );
        }
        if( enable_qp_ic_scores_ && tree.edge_count() != qp_ic_scores_.size() ) {
            throw std::runtime_error(
                "Size of QP-IC-Scores (" + std::to_string( qp_ic_scores_.size() ) +
                ") does not match numer of edges of the tree (" +
                std::to_string( tree.edge_count() ) + ")."
            );
        }
    }

    void edge_to_element( TreeEdge const& edge, NewickBrokerElement& element ) const
    {
        std::vector<std::string> edge_comments;

        if (enable_qp_ic_scores_) {
            assert(edge.index() < qp_ic_scores_.size());
            if (lq_ic_scores_[edge.index()] != std::numeric_limits<double>::infinity()) {
                edge_comments.push_back( "qp-ic:" + std::to_string(qp_ic_scores_[edge.index()]) );
            }
        }

        if (enable_lq_ic_scores_) {
            assert(edge.index() < lq_ic_scores_.size());
            if (lq_ic_scores_[edge.index()] != std::numeric_limits<double>::infinity()) {
                edge_comments.push_back( "lq-ic:" + std::to_string(lq_ic_scores_[edge.index()]) );
            }
        }

        if (enable_eqp_ic_scores_) {
            assert(edge.index() < eqp_ic_scores_.size());
            if (eqp_ic_scores_[edge.index()] != std::numeric_limits<double>::infinity()) {
                edge_comments.push_back( "eqp-ic:" + std::to_string(eqp_ic_scores_[edge.index()]) );
            }
        }

        if( ! edge_comments.empty() ) {
            element.comments.push_back( utils::join( edge_comments, ";" ));
        }
    }

    void register_with( NewickWriter& writer ) const
    {
        // Add node functions.
        writer.prepare_writing_plugins.push_back(
            [&]( tree::Tree const& tree, tree::NewickBroker& broker ) {
                prepare_writing( tree, broker );
            }
        );

        // Add edge functions.
        writer.edge_to_element_plugins.push_back(
            [&]( TreeEdge const& edge, NewickBrokerElement& element ) {
                edge_to_element( edge, element );
            }
        );
    }

    // -------------------------------------------------------------------------
    //     Data Members
    // -------------------------------------------------------------------------

private:

    bool enable_lq_ic_scores_  = false;
    bool enable_qp_ic_scores_  = false;
    bool enable_eqp_ic_scores_ = false;

    std::vector<double> lq_ic_scores_;
    std::vector<double> qp_ic_scores_;
    std::vector<double> eqp_ic_scores_;

};

// =================================================================================================
//     Default Tree Newick Writer
// =================================================================================================

/**
 * @brief
 */
class QuartetTreeNewickWriter
    : public NewickWriter
    , public tree::DefaultTreeNewickWriterPlugin
    , public QuartetNewickWriterPlugin
{
public:

    // -------------------------------------------------------------------------
    //     Constructor and Rule of Five
    // -------------------------------------------------------------------------

    QuartetTreeNewickWriter()
    {
        DefaultTreeNewickWriterPlugin::register_with( *this );
        QuartetNewickWriterPlugin::register_with( *this );
    }
};

} // namespace tree
} // namespace genesis

#endif // include guard
