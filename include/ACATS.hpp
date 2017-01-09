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

/*! \file ACATS.hpp
    \brief ACATS.hpp contains an implementation of the Asynchronous Compressed Affix Array Sampled.
           The implementations of the Compressed Suffix Tree must use Suffix Array intervals as nodes.
    \author Rodrigo Canovas
*/
#ifndef AFFIXTREE_ACATS_HPP
#define AFFIXTREE_ACATS_HPP

#include <stack>
#include <vector>
#include <queue>
#include <sdsl/suffix_trees.hpp>
#include <sdsl/csa_alphabet_strategy.hpp>
#include <sdsl/bit_vectors.hpp>
#include <sdsl/suffix_tree_helper.hpp>
//#include "affix_tree.hpp"

namespace atds {

    template<class t_cst = sdsl::cst_sct3<>,
            uint8_t t_width = 0>
    class ACATS {

    public: //rename some variable types
        typedef typename t_cst::size_type                           size_type;
        typedef typename t_cst::alphabet_category                   alphabet_category;
        typedef t_cst                                               cst_type;
        typedef typename t_cst::node_type                           s_node; //single node
        typedef typename t_cst::char_type                           char_type;
        typedef async_node<typename t_cst::node_type, size_type>    node_type;
        typedef cat_tag                                             index_category;
        typedef sdsl::bit_vector::rank_1_type                       rank_1_t;

    private: //variables
        sdsl::bit_vector                afl_sample[2];
        rank_1_t                        afl_sample_rank[2];
        sdsl::int_vector<>              affixtab[2];     //Respective affixlinktab
        cst_type                        CST[2];
        size_type                       n_text;

    public: //methods and constructors
        const typename cst_type::csa_type::alphabet_type::char2comp_type& char2comp    = CST[1].csa.char2comp;
        const typename cst_type::csa_type::alphabet_type::comp2char_type& comp2char    = CST[1].csa.comp2char;
        const typename cst_type::csa_type::alphabet_type::C_type&         C            = CST[1].csa.C;
        const typename cst_type::csa_type::alphabet_type::sigma_type&     sigma        = CST[1].csa.sigma;

        ACATS() { }

        ACATS(sdsl::cache_config &config, sdsl::cache_config &rev_config);

        void
        swap(ACATS &aff) {
            if (this != &aff) {
                for (uint8_t i = 0; i < 2; ++i) {
                    CST[i].swap(aff.CST[i]);
                    affixtab[i].swap(aff.affixtab[i]);
                    afl_sample[i].swap(aff.afl_sample[i]);
                    sdsl::util::swap_support(afl_sample_rank[i], aff.afl_sample_rank[i], &afl_sample[i],
                                             &(aff.afl_sample[i]));
                }
                std::swap(n_text, aff.n_text);
            }
        }

        void
        load(std::istream &in) {
            for (uint8_t i = 0; i < 2; ++i) {
                CST[i].load(in);
                affixtab[i].load(in);
                afl_sample[i].load(in);
                afl_sample_rank[i].load(in, &(afl_sample[i]));
            }
            n_text = CST[0].size();
        }

        size_type
        serialize(std::ostream &out, sdsl::structure_tree_node *v = NULL, std::string name = "") const {
            sdsl::structure_tree_node *child = sdsl::structure_tree::add_child(v, name, sdsl::util::class_name(*this));
            size_type written_bytes = 0;
            written_bytes += CST[0].serialize(out, child, "cst");
            std::cout << "aff_size : " << affixtab[0].size() << std::endl;
            written_bytes += affixtab[0].serialize(out, child, "affixtab");
            written_bytes += afl_sample[0].serialize(out, child, "bitmap");
            written_bytes += afl_sample_rank[0].serialize(out, child, "rank bitmap");
            written_bytes += CST[1].serialize(out, child, "rev_cst");
            written_bytes += affixtab[1].serialize(out, child, "rev_affixtab");
            written_bytes += afl_sample[1].serialize(out, child, "rev_bitmap");
            written_bytes += afl_sample_rank[1].serialize(out, child, "rev_rank bitmap");
            sdsl::structure_tree::add_size(child, written_bytes);
            return written_bytes;
        }

        void
        copy(const ACATS &aff) {
            for (uint8_t i = 0; i < 2; ++i) {
                CST[i] = aff.CST[i];
                affixtab[i] = aff.affixtab[i];
                afl_sample[i] = aff.afl_sample[i];
                afl_sample_rank[i] = aff.afl_sample_rank[i];
                afl_sample_rank[i].set_vector(&afl_sample[i]);
            }
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
            return CST[v.dir].size(v.node);
        }

        node_type
        root(size_type index = 0) const {
            return node_type(CST[index].root(), 0, 0, index);
        }

        bool
        is_leaf(node_type v) const {
            return rb(v) == lb(v);
        }

        size_type
        degree(const node_type &v) const {
            size_type aux_l = depth(v) - v.l;
            if (aux_l > 0)
                return 1;
            else
                return CST[v.dir].degree(v.node);
        }


        //!Calculates the index of the leftmost leaf in the corresponding suffix array.
        size_type
        lb(const node_type &v) const {
            return CST[v.dir].lb(v.node);
        }

        //! Calculates the index of the rightmost leaf in the corresponding suffix array.
        size_type
        rb(const node_type &v) const {
            return CST[v.dir].rb(v.node);
        }

        node_type
        parent(const node_type &v) const {
            if (v == root(v.dir))
                return root(v.dir);
            node_type pnode(CST[v.dir].parent(v.node), v.l, v.k, v.dir);
            pnode.l = depth(pnode);
            if (v.k > pnode.l)
                return root(v.dir);
            return pnode;
        }

        node_type
        sl(const node_type &v, bool dir = 0) {
            node_type slink = v;
            if (v == root(v.dir))
                return root(v.dir);
            if (v.dir == dir) {
                slink.node = CST[v.dir].sl(v.node);
                slink.l -= 1;
                if (v.k > 0)
                    slink.k -= 1;
                return slink;
            }
            else {
                if (v.k == 0) {
                    slink = get_reverse(slink);
                    return sl(slink, slink.dir);
                }
                else {
                    size_type l_aux = v.l;
                    s_node p = CST[v.dir].parent(v.node);
                    l_aux -= CST[v.dir].depth(p);
                    if (l_aux == 1)
                        slink.node = p;
                    slink.l -= 1;
                    return slink;
                }
            }

        }

        //! Returns the maximum string depth of node v (LCP value).
        size_type
        depth(const node_type &v) const {
            return CST[v.dir].depth(v.node);
        }

        //! Returns the d-th character (1-based indexing) of the edge-label pointing to v (starting from the root).
        char_type
        edge(node_type v, size_type d) const {
            return CST[v.dir].edge(v.node, d);
        }


        //! Get the child w of node v which edge label (v,w) starts with character c.
        // \sa child(node_type v, const char_type c, size_type &char_pos)
        node_type
        child(const node_type &v, const char_type c) const {
            size_type char_pos;
            return child(v, c, char_pos);
        }

        node_type
        child(const node_type &v, const char_type c, size_type &char_pos) const {
            if (is_leaf(v))
                return root(v.dir);
            //two cases: l == string depth of the queried node or l is lower and there is only one possible child
            node_type c_node = v;
            size_type aux_l = depth(v) - v.l; //if aux_l > 0 then there is only one possible child for v
            char_type letter = c;
            if (aux_l == 0) { //we are in an explicit node
                c_node.node = CST[v.dir].child(v.node, c, char_pos);
                if (c_node == root(v.dir))
                    return root(v.dir);
                c_node.k = v.k;
                aux_l = depth(c_node) - v.l;
            }
            else { //only one possible child, we are in an implicit node
                letter = edge(c_node, v.l);
                if (letter != c)
                    return root(c_node.dir);
            }
            c_node.l += aux_l;
            return c_node;
        }

        //note that v.k must be equal to 0
        node_type
        get_reverse(const node_type &v) const {
            size_type n_home = home(v), length = rb(v) - lb(v);
            size_type apos = 0;
            node_type rev;
            if (is_leaf(v)) {
                size_type pos_sa = CST[v.dir].csa[rb(v)];
                size_type n = CST[v.dir].csa.size();
                if (pos_sa != n - 1)
                    pos_sa = n - pos_sa - 2 - (v.l - 1);
                pos_sa = CST[1 - v.dir].csa.isa[pos_sa];  //inverse SA of the reverse
                rev = node_type(CST[1 - v.dir].node(pos_sa, pos_sa), v.l, 0, 1 - v.dir);
            }
            else {
                size_type l = depth(v);
                size_type k = l - v.l;
                if (afl_sample[v.dir][n_home] == 1) {
                    apos = afl_sample_rank[v.dir](n_home); //rank - 1
                    apos = affixtab[v.dir][apos];
                    rev = node_type(CST[1 - v.dir].node(apos, apos + length), l, k, 1 - v.dir);
                }
                else {
                    char_type c = edge(v, 1);
                    node_type sl_node = node_type(CST[v.dir].sl(v.node), 0, 0, v.dir);
                    n_home = home(sl_node);
                    apos = afl_sample_rank[v.dir](n_home); //rank - 1
                    apos = affixtab[v.dir][apos];
                    length = rb(sl_node) - lb(sl_node);
                    rev = node_type(CST[1 - v.dir].node(apos, apos + length), l, k, 1 - v.dir);
                    if (depth(rev) == depth(sl_node))
                        rev = node_type(CST[1 - v.dir].child(rev.node, c), l, k, 1 - v.dir);
                }
            }
            return rev;
        }

        //! Compute the suffix number of a leaf node v (or the value of the left most sa).
        size_type
        sn(node_type v) const {
            node_type lf(CST[v.dir].leftmost_leaf(v), 0, 0, v.dir);
            return CST[v.dir].sn(lf.node);
        }


        //! Return all possible node extension of one character for v, returning the resulting nodes
        std::vector<node_type>
        suffix_children(const node_type &v) const {
            std::vector<node_type> suff_ch;
            node_type c_node = v;
            size_type l_v = depth(v);
            if (v.l < l_v) {
                c_node.l += 1;//=l_v;
                suff_ch.push_back(c_node);
            }
            else {// compute all the children of v
                for (auto &n_child: CST[v.dir].children(v.node)) {
                    c_node = node_type(n_child, 0, v.k, v.dir);
                    c_node.l = l_v + 1; //depth(c_node);
                    suff_ch.push_back(c_node);
                }
            }
            return suff_ch;
        }

        //! Get a list of the respective reverse children of node v.
        std::vector<node_type>
        prefix_children(const node_type &v) const {
            std::vector<node_type> pref_ch;
            node_type c_node = v;
            if (v.k > 0) {
                c_node.k -= 1; //TODO: wrong v.k could take negative values
                pref_ch.push_back(c_node);
            }
            else { //v is a leaf of an explicit node
                c_node = get_reverse(v); //move to the reverse node
                pref_ch = suffix_children(c_node);
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
            if (s.dir == 1) {
                return fwd_search(end, begin, s, 1);
            }
            else {
                return bwd_search(begin, end, s, 0);
            }
        }

        /*!
        * Search for pattern in backward direction, starting with a given node (searchstate)
        * \param[in] character c to search for
        * \param[in,out] s searchstate (will be overwritten)
        * \return size of the new interval
        */
        size_type
        backward_search(char_type c, node_type &s) {
            std::vector<char_type> v(1, c);
            return backward_search(v.begin(), v.end(), s);
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
            if (s.dir == 0) {
                return fwd_search(begin, end, s, 0);
            }
            else {
                return bwd_search(end, begin, s, 1);
            }
        }

        /*!
		 * Search for pattern in forward direction, starting with a given node (searchstate)
		 * \param[in] character c to search for
		 * \param[in,out] s searchstate (will be overwritten)
		 * \return size of the new interval
		 */
        size_type
        forward_search(char_type c, node_type &s) {
            std::vector<char_type> v(1, c);
            return forward_search(v.begin(), v.end(), s);
        }


        //! Return the local forward search (taking as reference the direction of the node s) of the
        // string from begin to end starting from the node s.
        template<class t_pat_iter>
        size_t
        fwd_search(t_pat_iter begin, t_pat_iter end, node_type &s, size_type it_dir) {
            size_type aux_l = depth(s) - s.l;
            size_type shift = 0;
            char_type letter;
            t_pat_iter it = begin;
            size_type char_pos = 0;
            while (it != end) {
                if (aux_l > 0) {
                    shift = 0;
                    while (aux_l > 0 and it != end) { // only one possible child
                        it -= it_dir;
                        shift++;
                        letter = edge(s, s.l + shift);
                        if (letter != *it) {
                            s = root(s.dir);
                            return 0;
                        }
                        it += (1 - it_dir);
                        aux_l--;
                    }
                    s.l += shift;
                }
                else { // explicit node
                    it -= it_dir;
                    s.node = CST[s.dir].child(s.node, *it, char_pos);
                    if (s.node == CST[s.dir].root()) { //the child does not exist
                        s = root(s.dir);
                        return 0;
                    }
                    it += (1 - it_dir);
                    s.l++;
                    aux_l = depth(s) - s.l;
                }
            }
            if (s == root(s.dir) or s.l == 0) {
                s = root(s.dir);
                return 0;
            }
            return rb(s) - lb(s) + 1;
        }

        //! Return the local backward search (taking as reference the direction of the node s) of the string
        // from end to begin starting from the node s.
        template<class t_pat_iter>
        size_t
        bwd_search(t_pat_iter begin, t_pat_iter end, node_type &s, size_type it_dir) {
            char_type c;
            size_type text_pos = 0, max_length = 0;
            t_pat_iter it = end;
            node_type r_node;
            if (s.k > 0) { //scan until arriving to the closer explicit node or the answer
                while (begin != it and s.k > 0) {
                    it -= (1 - it_dir); 
                    c = edge(s, s.k); 
                    if (c != *it) {
                        s = root(s.dir);
                        return 0;
                    }
                    it += it_dir;
                    s.k -= 1;
                }
                if (it == begin)
                    return size(s);
            }
            //s.k == 0   move to the reverse node and do fwd_search
            r_node = get_reverse(s);
            s = r_node;     //note that this can not be a leaf
            return fwd_search(it, begin, s, 1 - it_dir);
        }

        void
        print_node_path(const node_type &v) {
            size_type pos = CST[v.dir].csa[lb(v)];
            std::string path = sdsl::extract(CST[v.dir].csa, pos, pos + v.l - 1);
            std::cout << path << std::endl;
        }

        size_type
        path_length(const node_type &v) {
            return v.l - v.k;
        }

    private:

        size_type
        home(const node_type &v) const {
            if (v == root(v.dir))
                return rb(v);
            if (rb(v) == n_text - 1 or CST[v.dir].lcp[lb(v)] >= CST[v.dir].lcp[rb(v) + 1])
                return lb(v);
            else
                return rb(v);
        }

        void
        create_affixlinktab(sdsl::cache_config &config, size_type index){
            size_type min_width = sdsl::bits::hi(n_text) + 1;
            std::cout << "min_width: " << min_width << std::endl;
            sdsl::int_vector<> tmp_aff(n_text, n_text, min_width);
            afl_sample[index]  = sdsl::bit_vector(n_text, 0);
            tmp_aff[n_text - 1] = 0; //set the root link
            afl_sample[index][n_text - 1] = 1;
            std::cout << "generate the array" << std::endl;
            node_type v = root(index);
            traverse_alink_louds(v, tmp_aff);
            sdsl::util::init_support(afl_sample_rank[index], &afl_sample[index]);
            size_type sample_size = 0;
            affixtab[index] = sdsl::int_vector<>(afl_sample_rank[index](n_text), n_text, min_width);
            for (size_type i = 0; i < n_text; ++i) {
                if(afl_sample[index][i] == 1){
                    affixtab[index][sample_size] = tmp_aff[i];
                    ++sample_size;
                }
            }
            std::cout << "Sample size: " << sample_size << std::endl;
            std::cout << "Finished ALINK " << index << std::endl;
        }

        void
        traverse_alink_louds(node_type &v, sdsl::int_vector<> &tmp_aff) {
            if (is_leaf(v))  // if v is a leave, v has no child
                return;
            s_node tree_node;
            node_type r_node, c_node;
            size_type h = 0, h_sl = 0;
            size_type max_rev = 10, last_l = 0, aux_l = 0;
            std::queue<s_node> to_check;
            to_check.push(v.node);
            while (!to_check.empty()) {
                tree_node = to_check.front();
                aux_l = CST[v.dir].depth(tree_node);
                if (aux_l > last_l + max_rev) { //is too big put at the end
                    to_check.pop();
                    to_check.push(tree_node);
                }
                else {
                    if (aux_l > last_l)
                        last_l = aux_l;
                    c_node = node_type(tree_node, 0, 0, v.dir);
                    h = home(c_node);
                    if (tmp_aff[h] == n_text) {
                        r_node = reverse_interval(c_node, h, tmp_aff);
                        tmp_aff[h] = lb(r_node);
                    }
                    h_sl = home(node_type(CST[v.dir].sl(tree_node),0,0,v.dir));
                    if (afl_sample[v.dir][h_sl] == 0)
                        afl_sample[v.dir][h] = 1;
                    for (auto &child: CST[v.dir].children(tree_node)) {
                        if (!CST[v.dir].is_leaf(child)) {
                            to_check.push(child);
                        }
                    }
                    to_check.pop();
                }
            }
        }

        node_type
        reverse_interval(node_type &v, size_type n_home, sdsl::int_vector<> &tmp_aff) {
            node_type rev, sl_node;
            char_type c;
            size_type length = rb(v) - lb(v);
            if (tmp_aff[n_home] != n_text) {
                rev = node_type(CST[1 - v.dir].node(tmp_aff[n_home], tmp_aff[n_home] + length),
                                0, 0, 1 - v.dir);
            }
            else {
                c = edge(v, 1);
                sl_node = node_type(CST[v.dir].sl(v.node), 0, 0, v.dir);
                rev = reverse_interval(sl_node, home(sl_node), tmp_aff);
                if (depth(rev) == depth(sl_node))
                    rev = node_type(CST[1 - v.dir].child(rev.node, c), 0, 0, 1 - v.dir);
                tmp_aff[n_home] = lb(rev);
            }
            return rev;
        }

    }; //ends class ACATS


    template<class t_cst, uint8_t t_width>
    ACATS<t_cst, t_width>::ACATS(sdsl::cache_config &config, sdsl::cache_config &rev_config) {
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
        //create the affixlink
        create_affixlinktab(config, 0);
        create_affixlinktab(rev_config, 1);
    }

} //end namespace

#endif //AFFIXTREE_ACATS_HPP
