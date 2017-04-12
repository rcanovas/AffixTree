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

#ifndef AFFIXTREE_AFFIX_TREE_BD_WT_HPP
#define AFFIXTREE_AFFIX_TREE_BD_WT_HPP


//#include "affix_tree.hpp"
#include <sdsl/int_vector.hpp>
#include <sdsl/rank_support.hpp>
#include <sdsl/wavelet_trees.hpp>


/*! This is an implementation of the bidirectional wavelet index.
 * It is parameterized by a wavelet tree (supporting extract ([]),
 * occ (rank(i,c)) and getBounds(lex_count) operations) and the number SampleDens specifying
 * how many samples of SA are explicitly stored.
 * The code presented is simply a mask over the original code of Thomas Schnattinger 
 * supplied within the sdsl library
 */

namespace atds {

    template<class t_int = sdsl::int_vector<>::size_type>
    struct searchstate {
        t_int occ_begin;
        t_int occ_end;
        t_int occ_forward_begin;
        t_int occ_forward_end;
        t_int length;

        searchstate() : occ_begin(0), occ_end(0), occ_forward_begin(0), occ_forward_end(0), length(0) { }

        searchstate(t_int i1, t_int j1, t_int i2, t_int j2, t_int l) : occ_begin(i1), occ_end(i2),
                                                                       occ_forward_begin(j1), occ_forward_end(j2),
                                                                       length(l) {}
        //! Copy constructor
        searchstate(const searchstate & s) = default;
        //! Move copy constructor
        searchstate(searchstate && s) = default;
        //! Assignment operator.
        searchstate & operator=(const searchstate & interval) = default;
        //! Move assignment
        searchstate & operator=(searchstate && interval) = default;

        //! Equality operator.
        /*! Two searchstate are equal if and only if all their corresponding member variables have the same values.
         * */
        bool operator==(const searchstate & interval)const {
            return occ_begin==interval.occ_begin and occ_end==interval.occ_end and
                    occ_forward_end==interval.occ_forward_end and occ_forward_begin==interval.occ_forward_begin
                   and length==interval.length;
        }

        //! Inequality operator.
        bool operator!=(const searchstate & interval)const {
            return !(*this==interval);
        }
    };


    template<class t_csa = sdsl::csa_wt< sdsl::wt_blcd<>,
                                         32,
                                         64,
                                         sdsl::sa_order_sa_sampling<>,
                                         sdsl::isa_sampling<>,
                                         sdsl::byte_alphabet > >
    class BidWT {


    public: //rename some variable types
        typedef typename t_csa::alphabet_type               alphabet_type;
        typedef typename alphabet_type::alphabet_category   alphabet_category;
        typedef typename alphabet_type::char_type           char_type;
        typedef typename t_csa::size_type                   size_type;
        typedef t_csa                                       csa_type;
        typedef bd_tag                                      index_category;
        typedef searchstate<size_type>                      node_type;


    private: //variables
        csa_type backward_index;
		csa_type forward_index;

    public: //methods and constructors
        const csa_type&             bw_index = backward_index;
        const csa_type&             fw_index = forward_index;
        const typename alphabet_type::char2comp_type& char2comp    = backward_index.char2comp;
        const typename alphabet_type::comp2char_type& comp2char    = backward_index.comp2char;
        const typename alphabet_type::C_type&         C            = backward_index.C;
        const typename alphabet_type::sigma_type&     sigma        = backward_index.sigma;

        //! Default constructor
        BidWT() {}

        //! Copy constructor
        BidWT(const BidWT& wt) {
            copy(wt);
        }

        //! Move constructor
        BidWT(BidWT&& wt) {
            *this = std::move(wt);
        }

        //! Construct from cache config
        BidWT(sdsl::cache_config& config, sdsl::cache_config &rev_config);

        void
        swap(BidWT &wt) {
            if (this != &wt) {
                backward_index.swap(wt.backward_index);
                forward_index.swap(wt.forward_index);
            }
        }

        //Size of the text being represented
        size_type
        size() const {
            return backward_index.size();
        }

        //! Calculate the number of leaves in the subtree rooted at node v.
        /*! \param v A valid node.
         *  \return The number of leaves in the subtree rooted at node v.
         */
        size_type
        size(const node_type& v) const {
            if (v.occ_end > v.occ_begin)
                return v.occ_end - v.occ_begin + 1;
            return 0;
        }

        node_type
        root() {
            node_type r;
            init_search_state(r);
            return r;
        }

        bool
        is_leaf(const node_type &v) const {
            return v.occ_begin == v.occ_end;
        }

        size_type
        degree(const node_type &v, const bool dir = 0) {
            char_type symb;
            size_type len = 0, d = 0;
            for (size_type c = sigma - 1; c > 0; c--) {
                node_type t = v;
                symb = comp2char[c];
                if (dir)
                    len = backward_search(symb, t);
                else
                    len = forward_search(symb, t);
                if (len != 0)
                    d++;
            }
            return d;
        }

        node_type
        parent(const node_type& v) {
            if (v == root())
                return root();
            else
                return root(); //not supported unless parent tables are stored;
        }

        node_type
        sl(const node_type& v) {
            if (v == root())
                return root();
            return root(); //not supported unless slink tables are stored;
        }


        /*!
		 * Initializes a searchstate to the searchstate of the empty word
		 * \param[out] s searchstate to be initialized
		 */
		void
        init_search_state(node_type &s) {
            s.occ_begin = 0;
            s.occ_end = backward_index.size() - 1;
            s.occ_forward_begin = 0;
            s.occ_forward_end = s.occ_end;
            s.length = 0;
        }

        std::vector<node_type>
        suffix_children(const node_type &v) {
            std::vector<node_type> suff_ch;
            //not implemented
            char_type symb;
            size_type len = 0, d = 0;
            for (size_type c = sigma - 1; c > 0; c--) {  //this could be improve by using the information of the bwt
                node_type t = v;
                symb = comp2char[c];
                len = forward_search(symb, t);
                if (len != 0)
                    suff_ch.push_back(t);
            }
            return suff_ch;
        }

        std::vector<node_type>
        prefix_children(const node_type &v) {
            std::vector<node_type> pref_ch;
            char_type symb;
            size_type len = 0, d = 0;
            for (size_type c = sigma - 1; c > 0; c--) {  //this could be improve by using the information of the bwt
                node_type t = v;
                symb = comp2char[c];
                len = backward_search(symb, t);
                if (len != 0)
                    pref_ch.push_back(t);
            }
            return pref_ch;
        }


        /*!
		 * Search for character c in forward direction, starting with a given node (searchstate)
		 * \param[in] pattern string of character to search for
		 * \param[in,out] s searchstate (will be overwritten)
		 * \return size of the new interval
		 */
        template<class t_pat_iter>
		size_type
        forward_search(t_pat_iter begin, t_pat_iter end, node_type &s) {
            size_type find;
            find = sdsl::bidirectional_search_forward(backward_index,
                                                          forward_index,
                                                          s.occ_begin,
                                                          s.occ_end,
                                                          s.occ_forward_begin,
                                                          s.occ_forward_end,
                                                          begin,
                                                          end,
                                                          s.occ_begin,
                                                          s.occ_end,
                                                          s.occ_forward_begin,
                                                          s.occ_forward_end);
            if (find == 0)
                s.length = 0;
            else
                s.length += std::distance(begin, end);
            return find;
        }

        /*!
		 * Search for pattern in forward direction, starting with a given node (searchstate)
		 * \param[in] character c to search for
		 * \param[in,out] s searchstate (will be overwritten)
		 * \return size of the new interval
		 */
		size_type
        forward_search(char_type c, node_type &s) {
            size_type find;
            find = sdsl::bidirectional_search(forward_index,
                                              s.occ_forward_begin,
                                              s.occ_forward_end,
                                              s.occ_begin,
                                              s.occ_end,
                                              c,
                                              s.occ_forward_begin,
                                              s.occ_forward_end,
                                              s.occ_begin,
                                              s.occ_end);
            if (find == 0)
                s.length = 0;
            else
                s.length += 1;
            return find;
        }

        /*!
		 * Search for pattern in backward direction, starting with a given node (searchstate)
		 * \param[in] pattern string of character to search for
		 * \param[in,out] s searchstate (will be overwritten)
		 * \return size of the new interval
		 */
        template<class t_pat_iter>
		size_type
        backward_search(t_pat_iter begin, t_pat_iter end, node_type &s) {
            size_type find;
            find = sdsl::bidirectional_search_backward(backward_index,
                                                           forward_index,
                                                           s.occ_begin,
                                                           s.occ_end,
                                                           s.occ_forward_begin,
                                                           s.occ_forward_end,
                                                           begin,
                                                           end,
                                                           s.occ_begin,
                                                           s.occ_end,
                                                           s.occ_forward_begin,
                                                           s.occ_forward_end);
            if (find == 0)
                s.length = 0;
            else
                s.length += std::distance(begin, end);
            return find;
        }

        /*!
		 * Search for pattern in backward direction, starting with a given node (searchstate)
		 * \param[in] character c to search for
		 * \param[in,out] s searchstate (will be overwritten)
		 * \return size of the new interval
		 */
		size_type
        backward_search(char_type c, node_type &s) {
            size_type find;
            find = sdsl::bidirectional_search(backward_index,
                                              s.occ_begin,
                                              s.occ_end,
                                              s.occ_forward_begin,
                                              s.occ_forward_end,
                                              c,
                                              s.occ_begin,
                                              s.occ_end,
                                              s.occ_forward_begin,
                                              s.occ_forward_end);
            if (find == 0)
                s.length = 0;
            else
                s.length += 1;
            return find;
        }


        /*!
		 * Extract the SA-values of the interval $[i..j]$
		 * \param[in] i begin of the interval (including)
		 * \param[in] j end of the interval (excluding)
		 * \param[out] occ result
		 */
		void
        extract_sa(size_t i, size_t j, sdsl::int_vector<> &occ) {
            occ.resize(j - i);
            for (size_type k = 0; k < j - i; k++) {
                occ[k] = backward_index[k];
            }
        }

        //! Assignment Operator.
        /*!
         *    Required for the Assignable Concept of the STL.
         */
        BidWT& operator=(const BidWT& wt) {
            if (this != &wt) {
                copy(wt);
            }
            return *this;
        }

        //! Assignment Move Operator.
        /*!
         *    Required for the Assignable Concept of the STL.
         */
        BidWT& operator=(BidWT&& wt) {
            if (this != &wt) {
                backward_index = std::move(wt.backward_index);
                forward_index = std::move(wt.forward_index);
            }
            return *this;
        }

        /*!
		 * Serializes this data structure to an output stream
		 * \param[in,out] out output stream to write to
		 * \return number of bytes written to out
		 */
		size_t
        serialize(std::ostream& out, sdsl::structure_tree_node* v=nullptr, std::string name="") const {
            sdsl::structure_tree_node* child = sdsl::structure_tree::add_child(v, name, sdsl::util::class_name(*this));
            size_type written_bytes = 0;
            written_bytes += backward_index.serialize(out, child, "backward index");
            written_bytes += forward_index.serialize(out, child, "forward index");
            sdsl::structure_tree::add_size(child, written_bytes);
            return written_bytes;
        }

		/*!
		 * Load this data structure from an input stream.
		 * \param[in,out] in input stream to read from
		 */
		void
        load(std::istream &in) {
            backward_index.load(in);
            forward_index.load(in);
        }

        void
        print_node_path(const node_type &v) {
            size_type pos = backward_index[v.occ_begin];
            std::string path = sdsl::extract(backward_index, pos, pos + v.length - 1);
            std::cout << path << std::endl;
        }

        size_type
        path_length(const node_type &v) {
            return v.length;
        }

    private: //local methods

        void
        copy(const BidWT& wt) {
            backward_index = wt.backward_index;
            forward_index = wt.forward_index;
        }



    };

    template<class t_csa>
    BidWT<t_csa>::BidWT(sdsl::cache_config& config, sdsl::cache_config &rev_config) {
        //assuming that t_csa is base in the BWT, if not, need to change directions
        {
            auto event = sdsl::memory_monitor::event("load csa");
            load_from_cache(backward_index, std::string(sdsl::conf::KEY_CSA) + "_" + sdsl::util::class_to_hash(backward_index), config);
        }
        {
            auto event = sdsl::memory_monitor::event("load reverse csa");
            load_from_cache(forward_index, std::string(sdsl::conf::KEY_CSA) + "_" + sdsl::util::class_to_hash(forward_index), rev_config);
        }
    }

template<class t_int>
inline std::ostream& operator<<(std::ostream& os, const searchstate<t_int>& interval)
{
    os<<"-["<<interval.occ_begin<<","<<interval.occ_end<<"] ["<<interval.occ_forward_begin<<
            ","<<interval.occ_forward_end<<"]  length: "<<interval.length;
    return os;
}


}

#endif //AFFIXTREE_AFFIX_TREE_BD_WT_HPP
