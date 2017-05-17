#ifndef QUARTET_NEWICK_WRITER_H_
#define QUARTET_NEWICK_WRITER_H_

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

#include "tree/default/newick_writer.hpp"
#include "tree/formats/newick/writer.hpp"

namespace genesis {
namespace tree {

// =================================================================================================
//     Quartet Tree Newick Writer Mixin
// =================================================================================================

/**
 * @brief
 */
template <typename Base>
class QuartetNewickWriterMixin : public Base
{
    // -------------------------------------------------------------------------
    //     Properties
    // -------------------------------------------------------------------------

public:

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
    //     Overridden Member Functions
    // -------------------------------------------------------------------------

protected:

    virtual void prepare_writing( tree::Tree const& tree, tree::NewickBroker& broker ) override
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

    virtual void edge_to_element( tree::TreeEdge const& edge, tree::NewickBrokerElement& element ) override
    {
        Base::edge_to_element( edge, element );
        
        //std::string edgeComment = "edge-id:" + std::to_string(edge.index());
        
        std::string edgeComment = "";
        
        bool addSemicolon = false;
        bool edgeHasInfo = false;
        
        if (enable_qp_ic_scores_) {
          assert(edge.index() < qp_ic_scores_.size());
          if (lq_ic_scores_[edge.index()] != std::numeric_limits<double>::infinity()) {
              edgeComment = edgeComment + "qp-ic:" + std::to_string(qp_ic_scores_[edge.index()]);
              edgeHasInfo = true;
          }
          addSemicolon = true;
        }
        
        if (enable_lq_ic_scores_) {
            assert(edge.index() < lq_ic_scores_.size());
            if (lq_ic_scores_[edge.index()] != std::numeric_limits<double>::infinity()) {
              if (addSemicolon) {
                edgeComment = edgeComment + ";lq-ic:" + std::to_string(lq_ic_scores_[edge.index()]);
              } else {
                edgeComment = edgeComment + "lq-ic:" + std::to_string(lq_ic_scores_[edge.index()]);
              }
              edgeHasInfo = true;
            }
            addSemicolon = true;
        }
        
        if (enable_eqp_ic_scores_) {
            assert(edge.index() < eqp_ic_scores_.size());
            if (eqp_ic_scores_[edge.index()] != std::numeric_limits<double>::infinity()) {
              if (addSemicolon) {
                edgeComment = edgeComment + ";eqp-ic:" + std::to_string(eqp_ic_scores_[edge.index()]);
              } else {
                edgeComment = edgeComment + "eqp-ic:" + std::to_string(eqp_ic_scores_[edge.index()]);
              }
              edgeHasInfo = true;
            }
            addSemicolon = true;
        }
        
        if (edgeHasInfo) {
          element.comments.push_back(edgeComment);
        }

        /*if (enable_lq_ic_scores_) {
            assert(edge.index() < lq_ic_scores_.size());
            if (lq_ic_scores_[edge.index()] != std::numeric_limits<double>::infinity()) {
                element.comments.push_back(std::to_string(lq_ic_scores_[edge.index()]));
            }
        }
        if (enable_eqp_ic_scores_) {
            assert(edge.index() < eqp_ic_scores_.size());
            if (eqp_ic_scores_[edge.index()] != std::numeric_limits<double>::infinity()) {
                element.comments.push_back(std::to_string(eqp_ic_scores_[edge.index()]));
            }
        }
        
        if (enable_qp_ic_scores_) {
            assert(edge.index() < qp_ic_scores_.size());
            if (qp_ic_scores_[edge.index()] != std::numeric_limits<double>::infinity()) {
                element.comments.push_back(std::to_string(eqp_ic_scores_[edge.index()]));
            }
        }*/
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
//     Quartet Tree Newick Writer
// =================================================================================================

using QuartetNewickWriter =
    QuartetNewickWriterMixin< tree::DefaultTreeNewickWriterMixin< tree::NewickWriter >>;

} // namespace tree
} // namespace genesis

#endif // include guard
