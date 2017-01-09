/* atds - affix tree data structures
    Copyright (C)2015-2016 Rodrigo Canovas

    This program is free software: you can redistribute it and/or modify
    it under the terms of the GNU General Public License as published by
    the Free Software Foundation, either version 3 of the License, or
    (at your option) any later version.

    This program is distributed in the hope that it will be useful,
    but WITHOUT ANY WARRANTY; without even the implied warranty of
    MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
    GNU General Public License for more details.

    You should have received a copy of the GNU General Public License
    along with this program.  If not, see http://www.gnu.org/licenses/ .
*/

#ifndef AFFIXTREE_AFFIX_TREE_HPP
#define AFFIXTREE_AFFIX_TREE_HPP

#include <utility>
#include <sdsl/construct.hpp>
#include <sdsl/suffix_trees.hpp>
#include <string>

struct affixtree_tag {};
struct cat_tag {};
struct bd_tag{};

const char KEY_SLi[] 		= "sli";
const char KEY_SLj[]        = "slj";

namespace atds {

    // Declaration of the Asynchronoues AffixTree's node type
    template<class node_type, class t_int = sdsl::int_vector<>::size_type>
    struct async_node {
        node_type node;     //!< The explicit node.
        t_int l;            //
        t_int k;            // value of the k-context of the node
        t_int dir;          // direction of the current tree (reverse = 1 or not = 0)

        //! Constructor
        async_node(node_type node=0, t_int l=0, t_int k=0, t_int dir=0):node(node),l(l),k(k),dir(dir) {};

        //! Copy constructor
        async_node(const async_node& iv) = default;
        //! Move copy constructor
        async_node(async_node&& iv) = default;

        //! Equality operator.
        /*! Two nodes are equal if and only if all their corresponding member variables have the same values.
         * */
        bool operator==(const async_node& v)const {
            return node==v.node; //and l==v.l and k==v.k and dir==v.dir;
        }

        //! Inequality operator.
        /*! Two nodes are not equal if and only if not all their corresponding member variables have the same values.
         * */
        bool operator!=(const async_node& v)const {
            return !(*this==v);
        }

        //! Assignment operator.
        async_node& operator=(const async_node& v) = default;
        //! Move assignment
        async_node& operator=(async_node&& v) = default;
    };

    template<class node_type, class t_int>
    inline std::ostream& operator<<(std::ostream& os, const async_node<node_type, t_int>& interval) {
        os << interval.node << " l: " << interval.l << " k: " << interval.k << " dir: " << interval.dir;
        return os;
    }


    template<class t_index>
    void
    construct(t_index &idx, std::string file, uint8_t num_bytes = 0) {
        sdsl::tMSS file_map;
        sdsl::cache_config config;
        if (sdsl::is_ram_file(file)) {
            config.dir = "@";
        }
        typename t_index::index_category index_tag;
        construct(idx, file, config, num_bytes, index_tag);
    }


//!Special case for Affix Tree structures
// Constructs an index object of type t_index for a text stored on disk.
/*!
 * \param idx       	t_index object.  Any assembly graph structure.
 * \param file      	Name of the text file. The representation of the file
 *                  	is dependent on the next parameter.
 * \
 * \param num_bytes 	If `num_bytes` equals 0, the file format is a serialized
 *				    	int_vector<>. Otherwise the file is interpreted as sequence
 *                  	of `num_bytes`-byte integer stored in big endian order.
 */
    template<class t_index>
    void
    construct(t_index &idx, const std::string &file, sdsl::cache_config &config, uint8_t num_bytes, affixtree_tag) {
        auto event = sdsl::memory_monitor::event("construct Affix Tree");
        const char *KEY_TEXT = sdsl::key_text_trait<t_index::alphabet_category::WIDTH>::KEY_TEXT;
        typedef sdsl::int_vector<t_index::alphabet_category::WIDTH> text_type;
        sdsl::cache_config rev_config;
        { // (1) create temporal direct structure needed
            {
                auto event = sdsl::memory_monitor::event("parse input text");
                // (1.1) check, if the text is cached
                if (!cache_file_exists(KEY_TEXT, config)) {
                    text_type text;
                    load_vector_from_file(text, file, num_bytes);
                    if (contains_no_zero_symbol(text, file)) {
                        append_zero_symbol(text);
                        store_to_cache(text, KEY_TEXT, config);
                    }
                }
                register_cache_file(KEY_TEXT, config);
            }
            {
                // (1.2) check, if the suffix array is cached
                auto event = sdsl::memory_monitor::event("SA");
                if (!cache_file_exists(sdsl::conf::KEY_SA, config)) {
                    sdsl::construct_sa<t_index::alphabet_category::WIDTH>(config);
                }
                register_cache_file(sdsl::conf::KEY_SA, config);
            }
            {
                // (1.3) check, if the lcp array is cached
                auto event = sdsl::memory_monitor::event("LCP");
                if (!cache_file_exists(sdsl::conf::KEY_LCP, config)) {
                    if (t_index::alphabet_category::WIDTH==8)
                        sdsl::construct_lcp_semi_extern_PHI(config);
                    else
                        sdsl::construct_lcp_PHI<t_index::alphabet_category::WIDTH>(config);
                }
                register_cache_file(sdsl::conf::KEY_LCP, config);
            }
        }
        {
            // (2) create the reverse text temporal structure needed
            rev_config = sdsl::cache_config(config.delete_files, config.dir, "rev_" + config.id);
            {
                 auto event = sdsl::memory_monitor::event("parse reverse text");
                 // (2.1) check, if the text is cached
                 if (!cache_file_exists(KEY_TEXT, rev_config)) {
                     text_type text;
                     load_vector_from_file(text, file, num_bytes);
                     std::reverse(text.begin(), text.end());
                     if (contains_no_zero_symbol(text, file)) {
                         append_zero_symbol(text);
                         store_to_cache(text, KEY_TEXT, rev_config);
                     }
                 }
                 register_cache_file(KEY_TEXT, rev_config);
            }
            {
                // (2.2) check, if the suffix array is cached
                auto event = sdsl::memory_monitor::event("rev_SA");
                if (!cache_file_exists(sdsl::conf::KEY_SA, rev_config)) {
                    sdsl::construct_sa<t_index::alphabet_category::WIDTH>(rev_config);
                }
                register_cache_file(sdsl::conf::KEY_SA, rev_config);
            }
            {
                // (2.3) check, if the lcp array is cached
                auto event = sdsl::memory_monitor::event("rev_LCP");
                if (!cache_file_exists(sdsl::conf::KEY_LCP, rev_config)) {
                    if (t_index::alphabet_category::WIDTH==8)
                        sdsl::construct_lcp_semi_extern_PHI(rev_config);
                    else
                        sdsl::construct_lcp_PHI<t_index::alphabet_category::WIDTH>(rev_config);
                }
                register_cache_file(sdsl::conf::KEY_LCP, rev_config);
            }
        }
        {
            // (3) create AffixTree
            auto event = sdsl::memory_monitor::event("AffixTree");
            t_index tmp(config, rev_config);  //call constructor of the index used
            tmp.swap(idx);
        }
        if(rev_config.delete_files){
            auto event = sdsl::memory_monitor::event("delete temporary reverse files");
            sdsl::util::delete_all_files(rev_config.file_map);
        }
        if (config.delete_files) {
            auto event = sdsl::memory_monitor::event("delete temporary files");
            sdsl::util::delete_all_files(config.file_map);
        }
    }

    //!Special case for Compressed Affix Tree structures
// Constructs an index object of type t_index for a text stored on disk.
/*!
 * \param idx       	t_index object.  Any assembly graph structure.
 * \param file      	Name of the text file. The representation of the file
 *                  	is dependent on the next parameter.
 * \
 * \param num_bytes 	If `num_bytes` equals 0, the file format is a serialized
 *				    	int_vector<>. Otherwise the file is interpreted as sequence
 *                  	of `num_bytes`-byte integer stored in big endian order.
 */
    template<class t_index>
    void
    construct(t_index &idx, const std::string &file, sdsl::cache_config &config, uint8_t num_bytes, cat_tag) {
        auto event = sdsl::memory_monitor::event("construct CAT");
        const char* KEY_TEXT = sdsl::key_text_trait<t_index::alphabet_category::WIDTH>::KEY_TEXT;
        typedef sdsl::int_vector<t_index::alphabet_category::WIDTH> text_type;
        sdsl::cst_tag cst_t;
        sdsl::cache_config rev_config;
        {  // (1) create temporal direct structure needed
            typename t_index::cst_type cst;
            if (!cache_file_exists(std::string(sdsl::conf::KEY_CST)+"_"+sdsl::util::class_to_hash(cst), config)) {
                sdsl::cache_config cst_config(false, config.dir, config.id, config.file_map);
                construct(cst, file, cst_config, num_bytes, cst_t);
                auto event = sdsl::memory_monitor::event("store CST");
                config.file_map = cst_config.file_map;
                store_to_cache(cst,std::string(sdsl::conf::KEY_CST)+"_"+sdsl::util::class_to_hash(cst), config);
            }
            register_cache_file(std::string(sdsl::conf::KEY_CST)+"_"+sdsl::util::class_to_hash(cst), config);
        }
        {  // (2) create the reverse text temporal structure needed
            rev_config = sdsl::cache_config(config.delete_files, config.dir, "rev_" + config.id);
            {
                auto event = sdsl::memory_monitor::event("parse reverse text");
                // (2.1) check, if the text is cached
                if (!cache_file_exists(KEY_TEXT, rev_config)) {
                    text_type text;
                    load_vector_from_file(text, file, num_bytes);
                    std::reverse(text.begin(), text.end());
                    if (contains_no_zero_symbol(text, file)) {
                        append_zero_symbol(text);
                        store_to_cache(text, KEY_TEXT, rev_config);
                    }
                }
                register_cache_file(KEY_TEXT, rev_config);
            }
            {
                typename t_index::cst_type cst;
                if (!cache_file_exists(std::string(sdsl::conf::KEY_CST)+"_"+sdsl::util::class_to_hash(cst), rev_config)) {
                    sdsl::cache_config cst_config(false, rev_config.dir, rev_config.id, rev_config.file_map);
                    construct(cst, file, cst_config, num_bytes, cst_t);
                    auto event = sdsl::memory_monitor::event("store reverse CST");
                    rev_config.file_map = cst_config.file_map;
                    store_to_cache(cst, std::string(sdsl::conf::KEY_CST)+"_"+sdsl::util::class_to_hash(cst), rev_config);
                }
                register_cache_file(std::string(sdsl::conf::KEY_CST)+"_"+sdsl::util::class_to_hash(cst), rev_config);
            }
        }
        {
            auto event = sdsl::memory_monitor::event("CAT");
            t_index tmp(config, rev_config);
            tmp.swap(idx);
        }
        if (config.delete_files) {
            auto event = sdsl::memory_monitor::event("delete temporary files");
            sdsl::util::delete_all_files(config.file_map);
            sdsl::util::delete_all_files(rev_config.file_map);
        }
    }

    //!Special case for Bidirectional WT structures
// Constructs an index object of type t_index for a text stored on disk.
/*!
 * \param idx       	t_index object.  Any assembly graph structure.
 * \param file      	Name of the text file. The representation of the file
 *                  	is dependent on the next parameter.
 * \
 * \param num_bytes 	If `num_bytes` equals 0, the file format is a serialized
 *				    	int_vector<>. Otherwise the file is interpreted as sequence
 *                  	of `num_bytes`-byte integer stored in big endian order.
 */

    template<class t_index>
    void
    construct(t_index &idx, const std::string &file, sdsl::cache_config &config, uint8_t num_bytes, bd_tag) {
        auto event = sdsl::memory_monitor::event("construct Bidirection WT");
        const char* KEY_TEXT = sdsl::key_text_trait<t_index::alphabet_category::WIDTH>::KEY_TEXT;
        const char* KEY_BWT  = sdsl::key_bwt_trait<t_index::alphabet_category::WIDTH>::KEY_BWT;

        typedef sdsl::int_vector<t_index::alphabet_category::WIDTH> text_type;
        sdsl::cache_config rev_config;
        sdsl::csa_tag csa_t;
        { // (1) create temporal direct structure needed
            typename t_index::csa_type csa;
            if (!cache_file_exists(std::string(sdsl::conf::KEY_CSA) + "_" + sdsl::util::class_to_hash(csa), config)) {
                sdsl::cache_config csa_config(false, config.dir, config.id, config.file_map);
                construct(csa, file, csa_config, num_bytes, csa_t);
                auto event = sdsl::memory_monitor::event("store CSA");
                config.file_map = csa_config.file_map;
                store_to_cache(csa, std::string(sdsl::conf::KEY_CSA) + "_" + sdsl::util::class_to_hash(csa), config);
            }
            register_cache_file(std::string(sdsl::conf::KEY_CSA) + "_" + sdsl::util::class_to_hash(csa), config);
        }
        // (2) create the reverse text temporal structure needed
        rev_config = sdsl::cache_config(config.delete_files, config.dir, "rev_" + config.id);
        {
            auto event = sdsl::memory_monitor::event("parse reverse text");
            // (2.1) check, if the text is cached
            if (!cache_file_exists(KEY_TEXT, rev_config)) {
                text_type text;
                load_vector_from_file(text, file, num_bytes);
                std::reverse(text.begin(), text.end());
                if (contains_no_zero_symbol(text, file)) {
                    append_zero_symbol(text);
                    store_to_cache(text, KEY_TEXT, rev_config);
                }
            }
            register_cache_file(KEY_TEXT, rev_config);
        }
        {
            typename t_index::csa_type csa;
            if (!cache_file_exists(std::string(sdsl::conf::KEY_CSA) + "_" + sdsl::util::class_to_hash(csa),
                                   rev_config)) {
                sdsl::cache_config csa_config(false, rev_config.dir, rev_config.id, rev_config.file_map);
                construct(csa, file, csa_config, num_bytes, csa_t);
                auto event = sdsl::memory_monitor::event("store reverse CSA");
                rev_config.file_map = csa_config.file_map;
                store_to_cache(csa, std::string(sdsl::conf::KEY_CSA) + "_" + sdsl::util::class_to_hash(csa),
                               rev_config);
            }
            register_cache_file(std::string(sdsl::conf::KEY_CSA) + "_" + sdsl::util::class_to_hash(csa), rev_config);
        }
        {   // (3) create the Bidirection WT
            auto event = sdsl::memory_monitor::event("Bidirection Wavelet Tree");
            t_index tmp(config, rev_config);  //call constructor of the index used
            tmp.swap(idx);
        }
        if(rev_config.delete_files){
            auto event = sdsl::memory_monitor::event("delete temporary reverse files");
            sdsl::util::delete_all_files(rev_config.file_map);
        }
        if (config.delete_files) {
            auto event = sdsl::memory_monitor::event("delete temporary files");
            sdsl::util::delete_all_files(config.file_map);
        }
    }

} //end namespace

#include "BidWT.hpp"
#include "AFA.hpp"
#include "ACAT.hpp"
#include "ACATS.hpp"
#include "ACATN.hpp"
#include "RACATN.hpp"
#include "SCAT.hpp"

#endif //AFFIXTREE_AFFIX_TREE_HPP
