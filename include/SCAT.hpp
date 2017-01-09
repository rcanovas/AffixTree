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

/*! \file SCAT.hpp
		\brief SCAT.hpp contains an implementation of the Synchronous Compressed Affix Array.
		\author Rodrigo Canovas
*/

#ifndef AFFIXTREE_AFFIX_TREE_CAT_SYN_HPP
#define AFFIXTREE_AFFIX_TREE_CAT_SYN_HPP

#include <stack>
#include <vector>
#include <sdsl/cst_sct3.hpp>
#include <sdsl/cst_fully.hpp>
#include <sdsl/cst_sada.hpp>
#include <sdsl/csa_alphabet_strategy.hpp>
#include <sdsl/bit_vectors.hpp>
#include <sdsl/suffix_tree_helper.hpp>
//#include "affix_tree.hpp"

namespace atds {

    // Declaration of the AffixTree's node type
    template<class node_type, class t_int = sdsl::int_vector<>::size_type>
    struct syn_node {
        node_type node[2];     // The explicit node [0] and reverse node [1]
        t_int l;            // Length of the path read

        //! Constructor
        syn_node(node_type node_t[2], t_int l_t = 0){
            node[0] = node_t[0]; node[1] = node_t[1];
            l =l_t;
        }

        //! Copy constructor
        syn_node(const syn_node &iv) = default;

        //! Move copy constructor
        syn_node(syn_node &&iv) = default;

        //! Equality operator.
        /*! Two syn-nodes are equal if and only if all their corresponding member variables have the same values.
         * */
        bool operator==(const syn_node &v) const {
            return node == v.node and l == v.l;
        }

        //! Inequality operator.
        /*! Two affix-nodes are not equal if and only if not all their corresponding member variables have the same values.
         * */
        bool operator!=(const syn_node &v) const {
            return !(*this == v);
        }

        //! Assignment operator.
        syn_node &operator=(const syn_node &v) = default;

        //! Move assignment
        syn_node &operator=(syn_node &&v) = default;
    };


    template<class t_cst = sdsl::cst_sct3<>,
            uint8_t t_width = 0>
    class SCAT {

    public: //rename some variable types
        typedef typename t_cst::size_type                       size_type;
        typedef typename t_cst::alphabet_category               alphabet_category;
        typedef t_cst                                           cst_type;
        typedef typename t_cst::char_type                       char_type;
        typedef typename t_cst::node_type                       s_node; //single node
        typedef syn_node<typename t_cst::node_type, size_type>  node_type;
        typedef cat_tag                                         index_category;

    private: //variables
        cst_type CST[2];
        size_type n_text;

    public: //methods and constructors
        const typename cst_type::csa_type::alphabet_type::char2comp_type &char2comp = CST[0].csa.char2comp;
        const typename cst_type::csa_type::alphabet_type::comp2char_type &comp2char = CST[0].csa.comp2char;
        const typename cst_type::csa_type::alphabet_type::C_type &C = CST[0].csa.C;
        const typename cst_type::csa_type::alphabet_type::sigma_type &sigma = CST[0].csa.sigma;


        /*default constructor*/
        SCAT() { }

        SCAT(sdsl::cache_config &config, sdsl::cache_config &rev_config);

        void
        swap(SCAT &aff) {
            if (this != &aff) {
                for (uint8_t i = 0; i < 2; ++i)
                    CST[i].swap(aff.CST[i]);
                std::swap(n_text, aff.n_text);
            }
        }

        void
        load(std::istream &in) {
            for (uint8_t i = 0; i < 2; ++i)
                CST[i].load(in);
            n_text = CST[0].size();
        }

        size_type
        serialize(std::ostream &out, sdsl::structure_tree_node *v = NULL, std::string name = "") const {
            sdsl::structure_tree_node *child = sdsl::structure_tree::add_child(v, name, sdsl::util::class_name(*this));
            size_type written_bytes = 0;
            written_bytes += CST[0].serialize(out, child, "cst");
            written_bytes += CST[1].serialize(out, child, "rev_cst");
            sdsl::structure_tree::add_size(child, written_bytes);
            return written_bytes;
        }

        void
        copy(const SCAT &aff) {
            for (uint8_t i = 0; i < 2; ++i)
                CST[i] = aff.CST[i];
            n_text = CST[0].size();
        }

        //Size of the text being represented
        size_type
        size() const {
            return n_text;
        }

        //! Calculate the number of leaves in the subtree rooted at node v.
        /*! \param v A valid node.
         *  \return The number of leaves in the subtree rooted at node v.
         */
        size_type
        size(const node_type &v) const {
            return CST[0].size(v.node); //should be the same number of leaves in both trees
        }

        node_type
        root() const {
            typename cst_type::node_type node[2];
            node[0] = CST[0].root(); node[1] = CST[1].root();
            return node_type(node, 0);
        }

        bool
        is_leaf(node_type v) const {
            return CST[0].is_leaf(v.node[0]); //if v.node[0] is a leaf then v.node[1] is a leaf
        }

        //!Calculates the index of the leftmost leaf in the corresponding suffix array.
        size_type
        lb(const node_type &v, const bool dir) const {
            return CST[dir].lb(v.node[dir]);
        }

        //! Calculates the index of the rightmost leaf in the corresponding suffix array.
        size_type
        rb(const node_type &v, const bool dir = 0) const {
            return CST[dir].rb(v.node[dir]);
        }

        //! Returns the string depth of node v (LCP value) in the corresponding tree (note that this could
        // be different than v.l).
        size_type
        depth(const node_type &v, const bool dir = 0) const {
            return CST[dir].depth(v.node[dir]);
        }

        //! Returns the d-th character (1-based indexing) of the edge-label pointing to v.
        char_type
        edge(node_type v, size_type d, const bool dir = 0) const {
            return CST[dir].edge(v.node[dir], d);
        }

        //! Returns the node depth of node v
        size_type
        node_depth(node_type v, const bool dir = 0) const {
            return CST[dir].node_depth(v.node[dir]);
        }

        //! return the node resulting after moving to the parent of the corresponding node of v.
        node_type
        parent(const node_type& v, const bool dir) const {
            if (v == root())
                return root();
            node_type p = v;
            size_type aux_l = 0;
            p.node[dir] = CST[dir].parent(p.node[dir]);
            p.l = depth(p, dir);
            aux_l = v.l - p.l;
            while (aux_l > 0) {
                p.node[1 - dir] = CST[1 - dir].sl(p.node[1 - dir]);
                --aux_l;
            }
            return p;
        }

        size_type
        degree(const node_type &v, const bool dir = 0) const {
            size_type aux_l = depth(v,dir) - v.l;
            if (aux_l > 0)
                return 1;
            else
                return CST[dir].degree(v.node[dir]);
        }

        //! Compute the suffix link of node v.
        node_type
        sl(const node_type &v, const bool dir = 0) const {
            s_node p;
            node_type s = v;
            s.node[dir] = CST[dir].sl(s.node[dir]);
            if (s.node[dir] == CST[dir].root())
                return root();
            size_type l_aux = v.l;
            p = CST[1 - dir].parent(v.node[1 - dir]);

            l_aux -= CST[1 - dir].depth(p);
            if (l_aux == 1)
                s.node[1 - dir] = p;
            s.l -= 1;
            return s;
        }

        //! Get the child w of corresponding node of v which edge label (v.node(or _rev),w) starts with character c.
        node_type
        child(const node_type &v, const char_type c, const bool dir) const {
            size_type char_pos;
            return child(v, c, char_pos, dir);
        }

        node_type
        child(const node_type &v, const char_type c, size_type &char_pos, const bool dir) const {
            if (is_leaf(v))
                return root();
            //two cases: l == string depth of the queried node or l is lower and there is only one possible child
            node_type c_node = v;
            size_type aux_l = depth(v, dir) - v.l; //if aux_l>0 then there is only one possible child for v
            size_type shift = 0;
            char_type letter = c;
            if (aux_l == 0) { //we are in an explicit node
                c_node.node[dir] = CST[dir].child(c_node.node[dir], c, char_pos);
                if (c_node.node[dir] == CST[dir].root()) //the child does not exist
                    return root();
                aux_l = depth(c_node, dir) - v.l;
            }
            else { //only one possible child, we are in an implicit node
                letter = edge(c_node, v.l, dir);
                if (letter != c)
                    return root();
            }
            while (aux_l != shift) { //length of the path between the node and the child
                letter = edge(c_node, v.l + shift, dir);
                c_node.node[1 - dir] = wl_all(c_node.node[1 - dir], letter, 1 - dir);
                //c_node.node[1 - dir] = CST[1 - dir].wl(c_node.node[1 - dir], letter);
                shift++;
            }
            c_node.l += aux_l;
            return c_node;
        }

        //! Compute the Weiner link of the corresponding node of v using character c.
        node_type
        wl(node_type v, const char_type c, const bool dir) const {
            node_type wl_node = v;
            size_type shift = 0, aux_l = 0;
            char_type letter = c;
            wl_node.node[dir] = CST[dir].wl(wl_node.node[dir], c);
            if (wl_node.node[dir] == CST[dir].root()) //the wl does not exist
                return root();
            //move in the reverse node
            aux_l = depth(wl_node, (bool)(1 - dir)) - v.l;
            if (aux_l == 0)
                wl_node.node[1 - dir] = CST[1 - dir].child(wl_node.node[1 - dir], c);
            wl_node.l += 1;
            return wl_node;
        }

        //! Compute the suffix number of a leaf node v (or the value of the left most sa).
        size_type
        sn(node_type v, const bool dir) const {
            typename cst_type::node_type leftmost;
            leftmost = CST[dir].leftmost_leaf(v.node[dir]);
            return CST[dir].sn(leftmost);
        }

        //! Return all possible node extension of one character for v, returning the resulting nodes
        std::vector<node_type>
        suffix_children(const node_type &v) {
            std::vector<node_type> suff_ch;
            node_type c_node = v;
            char_type letter;
            size_type l_v = depth(v,0), aux_l, shift;
            if (v.l < l_v) {
                letter = edge(v, v.l + 1, 0);
                //c_node.node[1] = CST[1].wl(c_node.node[1], letter);
                c_node.node[1] = wl_all(c_node.node[1], letter, 1);
                c_node.l += 1; //=aux_l; //update length
                suff_ch.push_back(c_node);
            }
            else { // compute all the children of v
                for (auto &n_child: CST[0].children(v.node[0])) {
                    c_node.node[0] = n_child;
                    //compute the reverse node
                    letter = edge(c_node, v.l + 1, 0);
                    //c_node.node[1] = CST[1].wl(c_node.node[1], letter);
                    c_node.node[1] = wl_all(c_node.node[1], letter, 1);
                    c_node.l += 1;//=aux_l;
                    suff_ch.push_back(c_node);
                }
            }
            return suff_ch;
        }

        //! Get a list of the respective reverse children of node v.
        std::vector<node_type>
        prefix_children(const node_type &v) {
            std::vector<node_type> pref_ch;
            node_type c_node = v;
            char_type letter;
            size_type l_v = depth(v,1), aux_l, shift;
            if (v.l < l_v) {
                letter = edge(v, v.l + 1, 1);
                //c_node.node[0] = CST[0].wl(c_node.node[0], letter);
                c_node.node[0] = wl_all(c_node.node[0], letter, 0);
                c_node.l += 1; //update length
                pref_ch.push_back(c_node);
            }
            else{
                for (auto &n_child: CST[1].children(v.node[1])) {
                    c_node.node[1] = n_child;
                    //compute the reverse node
                    letter = edge(c_node, v.l + 1, 1);
                   // c_node.node[0] = CST[0].wl(c_node.node[0], letter);
                    c_node.node[0] = wl_all(c_node.node[0], letter, 0);
                    c_node.l += 1;
                    pref_ch.push_back(c_node);
                }
            }
            return pref_ch;
        }

        /*!
         * Search for pattern in backward direction, starting with a given node
		 * \param[in] pattern string of character to search for
		 * \param[in,out] s node (will be overwritten)
		 * \return size of the new interval
		 */
        template<class t_pat_iter>
        size_t
        backward_search(t_pat_iter begin, t_pat_iter end, node_type &s) {
            size_type aux_l = depth(s, 1) - s.l;
            size_type shift = 0;
            char_type letter;
            t_pat_iter it = begin;
            size_type char_pos = 0;
            while (it != end) {
                if (aux_l > 0) { //if aux_l>0 then there is only one possible child
                    shift = 0;
                    while (aux_l > 0 and it != end) { // only one possible child
                        shift++;
                        letter = edge(s, s.l + shift, 1);
                        if (letter != *it) {
                            s = root();
                            return 0;
                        }
                        aux_l--;
                        //s.node[0] = CST[0].wl(s.node[0], letter); //move the reverse
                        s.node[0] = wl_all(s.node[0], letter, 0);
                        it++;
                    }
                    s.l += shift;
                }
                else {  //explicit node
                    s.node[1] = CST[1].child(s.node[1], *it, char_pos);
                    if (s.node[1] == CST[1].root()) { //the child does not exist
                        s = root();
                        return 0;
                    }
                    //s.node[0] = CST[0].wl(s.node[0], *it); //move the reverse
                    s.node[0] = wl_all(s.node[0], *it, 0);
                    it++;
                    s.l++;
                    aux_l = depth(s, 1) - s.l;
                }
            }
            if (s == root() or s.l == 0) {
                s = root();
                return 0;
            }
            return rb(s, 0) - lb(s, 0) + 1;
        }

        /*!
        * Search for pattern in backward direction, starting with a given node (searchstate)
        * \param[in] character c to search for
        * \param[in,out] s searchstate (will be overwritten)
        * \return size of the new interval
        */
        size_type
        backward_search(char_type c, node_type &s) {
            s = wl(s, c, 0);
            if (s == root() or s.l == 0) {
                s = root();
                return 0;
            }
            return rb(s, 0) - lb(s, 0) + 1;
        }

        /*!
		 * Search for pattern in forward direction, starting with a given node
		 * \param[in] pattern string of character to search for
		 * \param[in,out] s node (will be overwritten)
		 * \return size of the new interval
		 */
        template<class t_pat_iter>
        size_t
        forward_search(t_pat_iter begin, t_pat_iter end, node_type &s) {
            size_type aux_l = depth(s, 0) - s.l;
            size_type shift = 0;
            char_type letter;
            t_pat_iter it = begin;
            size_type char_pos = 0;
            while (it != end) {
                if (aux_l > 0) {
                    shift = 0;
                    while (aux_l > 0 and it != end) { // only one possible child
                        shift++;
                        letter = edge(s, s.l + shift, 0);
                        if (letter != *it) {
                            s = root();
                            return 0;
                        }
                        aux_l--;
                        //s.node[1] = CST[1].wl(s.node[1], letter); //move the reverse
                        s.node[1] = wl_all(s.node[1], letter, 1);
                        it++;
                    }
                    s.l += shift;
                }
                else { // explicit node
                    s.node[0] = CST[0].child(s.node[0], *it, char_pos);
                    if (s.node[0] == CST[0].root()) { //the child does not exist
                        s = root();
                        return 0;
                    }
                    //s.node[1] = CST[1].wl(s.node[1], *it);//move the reverse
                    s.node[1] = wl_all(s.node[1], *it, 1);
                    it++;
                    s.l++;
                    aux_l = depth(s, 0) - s.l;
                }
            }
            if (s == root() or s.l == 0) {
                s = root();
                return 0;
            }
            return rb(s, 0) - lb(s, 0) + 1;
        }

        /*!
		 * Search for pattern in forward direction, starting with a given node (searchstate)
		 * \param[in] character c to search for
		 * \param[in,out] s searchstate (will be overwritten)
		 * \return size of the new interval
		 */
        size_type
        forward_search(char_type c, node_type &s) {
            s = wl(s, c, 1);
            if (s == root() or s.l == 0) {
                s = root();
                return 0;
            }
            return rb(s, 0) - lb(s, 0) + 1;
        }

        void
        print_node_path(const node_type &v, bool dir = 0) {
            size_type pos = CST[dir].csa[lb(v, dir)];
            std::string path = sdsl::extract(CST[dir].csa , pos, pos + v.l - 1);
            std::cout << path << std::endl;
        }

        size_type
        path_length(const node_type &v) {
            return v.l;
        }

    private:

        //! Weinier_link also for the leafs only for the given CST.
        /*!
         * \param v A valid node of a cst_cn.
         * \param c The character which should be prepended to the string of the current node.
         * \return  root() if the Weiner link of (v, c) does not exist,
         *          otherwise the Weiner link is returned.
         */
        s_node
        wl_all(const s_node& v, const char_type c, bool dir = 0) {
            size_type l, r;
            l = CST[dir].lb(v);
            r = CST[dir].rb(v);
            if (l == r) {  // is a leaf only one possible position
                size_type text_pos = CST[dir].csa[l];
                if (text_pos != 0) {
                    size_type sa_pos = CST[dir].csa.isa[text_pos - 1];
                    s_node w = CST[dir].node(sa_pos, sa_pos);
                    if (CST[dir].edge(w, 1) == c)
                        return w;
                }
                return CST[dir].root();
            }
            return CST[dir].wl(v, c);
        }


    }; //end class ACAT


    template<class t_cst, uint8_t t_width>
    SCAT<t_cst, t_width>::SCAT(sdsl::cache_config &config, sdsl::cache_config &rev_config) {
        {
            auto event = sdsl::memory_monitor::event("load cst");
            load_from_cache(CST[0], std::string(sdsl::conf::KEY_CST) + "_" + sdsl::util::class_to_hash(CST[0]), config);
            n_text = CST[0].size();
        }
        {
            auto event = sdsl::memory_monitor::event("load reverse cst");
            load_from_cache(CST[1], std::string(sdsl::conf::KEY_CST) + "_" + sdsl::util::class_to_hash(CST[1]),
                            rev_config);
        }
    }

    template<class node_type, class t_int>
    inline std::ostream &operator<<(std::ostream &os, const syn_node<node_type, t_int> &interval) {
        os << "node: " << interval.node[0] << "  rev_node: " << interval.node[1] << " l: " << interval.l;
        return os;
    }

} //end namespace

#endif //AFFIXTREE_AFFIX_TREE_CAT_SYN_HPP
