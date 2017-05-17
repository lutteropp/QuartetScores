#ifndef RANGE_MINIMUM_QUERY_H_
#define RANGE_MINIMUM_QUERY_H_

#include "utils/math/matrix.hpp"

#include <math.h>
#include <stddef.h>
#include <stdint.h>
#include <vector>

namespace external {

using namespace genesis;
using namespace utils;

/**
 * @brief Class that allows to find the index of the minimum element within an interval of an array.
 *
 * The implementation is based on the Succinct RMQ implementation
 * (https://www.bio.ifi.lmu.de/forschung/succinct/#software) by Johanens Fischer
 * (http://www.cs.tu-dortmund.de/nps/de/Home/Personen/F/Fischer__Johannes.html).
 * Most of the implementation is used as-is, with only some minor changes to the data types in order
 * to make it work with genesis.
 *
 * Currently, the original implementation is licensed under GPL v2, while genesis used GPL v3.
 * Thus, they are incompatible. As a consequence, this code currently cannot be distributed with
 * genesis. This includes the code in the corresponding implementation file `range_minimum_query.cpp`
 * in this directory. Both files and their code are for internal testing purposes only.
 */
class RangeMinimumQuery
{
public:
    // -------------------------------------------------------------------------
    //     Typedefs and Enums
    // -------------------------------------------------------------------------

    /**
     * @brief Data type of the array for which we want to run RMQ queries.
     *
     * Currently, this is fixed to a signed 32-bit integer. If a wider range is needed,
     * many internal functions need to be adapted first.
     */
    using IntType       = int32_t;

    using SuccinctType  = unsigned char;
    using BlockTypeType = unsigned short;

    // -------------------------------------------------------------------------
    //     Constructors and Rule of Five
    // -------------------------------------------------------------------------

    RangeMinimumQuery( std::vector<IntType> const& array );
    RangeMinimumQuery( std::vector<IntType>&&      array );

    ~RangeMinimumQuery() = default;

    RangeMinimumQuery( RangeMinimumQuery const& ) = default;
    RangeMinimumQuery( RangeMinimumQuery&& )      = default;

    RangeMinimumQuery& operator= ( RangeMinimumQuery const& ) = default;
    RangeMinimumQuery& operator= ( RangeMinimumQuery&& )      = default;

    // -------------------------------------------------------------------------
    //     Member Functions
    // -------------------------------------------------------------------------

    size_t query( size_t i, size_t j ) const;

    // -------------------------------------------------------------------------
    //     Internal Inline Functions
    // -------------------------------------------------------------------------

private:

    /**
     * @brief Return the microblock-number of entry i.
     */
    inline size_t microblock_( size_t i ) const
    {
        return i / micro_size_;
    }

    /**
     * @brief Return the block-number of entry i.
     */
    inline size_t block_( size_t i ) const
    {
        return i / block_size_;
    }

    /**
     * @brief Return the superblock-number of entry i.
     */
    inline size_t superblock_( size_t i ) const
    {
        return i / super_size_;
    }

    /**
     * @brief Return least signigicant bit in constant time (change for 64bit version).
     */
    inline size_t lsb_( SuccinctType v ) const
    {
        return lsb_table_256_[v];
    }

    /**
     * @brief Because M just stores offsets (rel. to start of block),
     * this method re-calculates the true index:
     */
    inline size_t m_(size_t k, size_t block) const
    {
        return m_matrix_(k, block)+(block*block_size_);
    }

    // -------------------------------------------------------------------------
    //     Internal Functions
    // -------------------------------------------------------------------------

    void init_();

    size_t log2fast_( size_t v ) const;
    SuccinctType clearbits_( SuccinctType n, size_t x ) const;

    // -------------------------------------------------------------------------
    //     Lookup Tables
    // -------------------------------------------------------------------------

private:

    static const size_t       catalan_numbers_[17][17];
    static const char         log_table_256_[256];
    static const char         lsb_table_256_[256];
    static const SuccinctType highest_bits_set_[8];

    // -------------------------------------------------------------------------
    //     Internal Members
    // -------------------------------------------------------------------------

private:

    // data
    std::vector<IntType> array_;

    // table M for the out-of-block queries (contains indices of block-minima)
    Matrix<SuccinctType> m_matrix_;

    // table M' for superblock-queries (contains indices of block-minima)
    Matrix<size_t> m_prime_;

    // type of blocks
    std::vector<BlockTypeType> block_types_;

    // precomputed in-block queries
    Matrix<SuccinctType> precomputed_queries_;

    // microblock size
    size_t micro_size_;

    // block size
    size_t block_size_;

    // superblock size
    size_t super_size_;

    // If the data array is too small, we use the naive approach instead.
    bool naive_ = false;

};

} // namespace external

#endif // include guard
