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
/*! \file AFA.hpp
    \brief AFA.hpp contains an implementation of Strothmann's Affix Array.
    \author Rodrigo Canovas
*/

#ifndef AFFIX_ARRAY_HPP
#define AFFIX_ARRAY_HPP

#include <stack>
#include <queue>
#include <vector>
#include <sdsl/csa_alphabet_strategy.hpp>
#include <sdsl/lcp_bitcompressed.hpp>
#include <sdsl/bit_vectors.hpp>
#include <sdsl/suffix_tree_helper.hpp>
//#include "affix_tree.hpp"

namespace atds {

    // Declaration of the AffixTree's node type
    template<class t_int = sdsl::int_vector<>::size_type>
    struct interval_sa {
        t_int i;        //!< The left border of the lcp-interval \f$\ell-[left..right]\f$.
        t_int j;        //!< The right border of the lcp-interval \f$\ell-[left..right]\f$.
        t_int l;        // value of the length of the path cover from the root to the interval [i+1,j]
        t_int k;        // value of the k-context of the node
        t_int n_l;      // negative length (i.e path that goes backward from from the root)
        bool dir;      // direction of the current tree (reverse = 1 or not = 0)

        //! Constructor
        interval_sa(t_int i=0, t_int j=0, t_int l=0, t_int k=0, bool dir=0, t_int n_l=0): i(i), j(j), l(l), k(k),
                                                                                             dir(dir), n_l(n_l){};

        //! Copy constructor
        interval_sa (const interval_sa  & iv) = default;
        //! Move copy constructor
        interval_sa(interval_sa && iv) = default;

        bool operator<(const interval_sa & interval)const {
            if (i!=interval.i)
                return i<interval.i;
            return j<interval.j;
        }

        //! Equality operator.
        /*! Two lcp-intervals are equal if and only if all their corresponding member variables have the same values.
         * */
        bool operator==(const interval_sa & interval)const {
            return i==interval.i and j==interval.j;
            //and l==interval.l
            //       and k==interval.k and dir==interval.dir and n_l==interval.n_l;
        }

        //! Inequality operator.
        /*! Two lcp-intervals are not equal if and only if not all their corresponding member variables have the same values.
         * */
        bool operator!=(const interval_sa & interval)const {
            return !(*this==interval);
        }

        //! Assignment operator.
        interval_sa & operator=(const interval_sa & interval) = default;
        //! Move assignment
        interval_sa & operator=(interval_sa && interval) = default;
    };


    template<class t_int>
    inline std::ostream& operator<<(std::ostream& os, const interval_sa<t_int>& interval);



    template<class t_alphabet_strat=sdsl::byte_alphabet, uint8_t t_width = 0>
    class AFA {

    public: //rename some variable types

        typedef t_alphabet_strat                            alphabet_type;
        typedef typename alphabet_type::alphabet_category   alphabet_category;
        typedef affixtree_tag                               index_category;
        typedef typename alphabet_type::char_type           char_type;
        typedef sdsl::int_vector<>::size_type               size_type;
        typedef interval_sa<size_type>                      node_type;

    private: //variables

        sdsl::int_vector<> Text;            //original text
        sdsl::int_vector<> SA[2];           //2 int vectors of size n*log(n) for SA and SA^r
        sdsl::lcp_byte<> LCP[2];            //2 int vectors of size n*maxlcp for LCP and LCP^r
        sdsl::int_vector<> childtab[2];     //Respective childtab info
        sdsl::int_vector<> affixtab[2];     //Respective affixlinktab
        alphabet_type      m_alphabet;      //alphabet policy of the text

    public: //methods and constructors
        const typename alphabet_type::char2comp_type& char2comp    = m_alphabet.char2comp;
        const typename alphabet_type::comp2char_type& comp2char    = m_alphabet.comp2char;
        const typename alphabet_type::C_type&         C            = m_alphabet.C;
        const typename alphabet_type::sigma_type&     sigma        = m_alphabet.sigma;


        /*default constructor*/
        AFA() { }

        AFA(sdsl::cache_config &config, sdsl::cache_config &rev_config);

        void
        swap(AFA &aff) {
            if (this != &aff) {
                Text.swap(aff.Text);
                for (uint8_t i = 0; i < 2; ++i) {
                    SA[i].swap(aff.SA[i]);
                    LCP[i].swap(aff.LCP[i]);
                    childtab[i].swap(aff.childtab[i]);
                    affixtab[i].swap(aff.affixtab[i]);
                }
                m_alphabet.swap(aff.m_alphabet);
            }
        }

        void
        load(std::istream &in) {
            Text.load(in);
            for (uint8_t i = 0; i < 2; ++i) {
                SA[i].load(in);
                LCP[i].load(in);
                childtab[i].load(in);
                affixtab[i].load(in);
            }
            m_alphabet.load(in);
        }

        //! Assignment Operator.
        /*!
         *    Required for the Assignable Concept of the STL.
         */
        AFA& operator=(const AFA& aff) {
            if (this != &aff) {
                copy(aff);
            }
            return *this;
        }

        //! Assignment Move Operator.
        /*!
         *    Required for the Assignable Concept of the STL.
         */
        AFA& operator=(AFA&& aff) {
            if (this != &aff) {
                Text = std::move(aff.Text);
                for (uint8_t i = 0; i < 2; ++i) {
                    SA[i] = std::move(aff.SA[i]);
                    LCP[i] = std::move(aff.LCP[i]);
                    childtab[i] = std::move(aff.childtab[i]);
                    affixtab[i] = std::move(aff.affixtab[i]);
                }
                m_alphabet     = std::move(aff.m_alphabet);
            }
            return *this;
        }

        size_type
        serialize(std::ostream &out, sdsl::structure_tree_node *v = NULL, std::string name = "") const {
            sdsl::structure_tree_node *child = sdsl::structure_tree::add_child(v, name, sdsl::util::class_name(*this));
            size_type written_bytes = 0;
            written_bytes += Text.serialize(out, child, "text");
            written_bytes += SA[0].serialize(out, child, "sa");
            written_bytes += LCP[0].serialize(out, child, "lcp");
            written_bytes += childtab[0].serialize(out, child, "childtab");
            written_bytes += affixtab[0].serialize(out, child, "affixtab");
            written_bytes += SA[1].serialize(out, child, "rev_sa");
            written_bytes += LCP[1].serialize(out, child, "rev_lcp");
            written_bytes += childtab[1].serialize(out, child, "rev_childtab");
            written_bytes += affixtab[1].serialize(out, child, "rev_affixtab");
            written_bytes += m_alphabet.serialize(out, child, "alphabet");
            sdsl::structure_tree::add_size(child, written_bytes);
            return written_bytes;
        }

        void
        copy(const AFA &aff) {
            Text = aff.Text;
            for (uint8_t i = 0; i < 2; ++i) {
                SA[i] = aff.SA[i];
                LCP[i] = aff.LCP[i];
                childtab[i] = aff.childtab[i];
                affixtab[i] = aff.affixtab[i];
            }
            m_alphabet     = aff.m_alphabet;
        }

        //Size of the text being represented
        size_type
        size() const {
            return Text.size();
        }

         //! Calculate the number of leaves in the subtree rooted at node v.
        /*! \param v A valid node.
         *  \return The number of leaves in the subtree rooted at node v.
         */
        size_type
        size(const node_type& v) const {
            return v.j - v.i + 1;
        }

        node_type
        root(bool index = 0) const {
            return node_type(0, size() - 1, 0, 0, index);
        }

        bool
        is_leaf(const node_type &v) const {
            return v.i == v.j;
        }

        //!Calculates the index of the leftmost leaf in the corresponding suffix array.
        size_type
        lb(const node_type& v)const {
            return v.i;
        }

        //! Calculates the index of the rightmost leaf in the corresponding suffix array.
        size_type
        rb(const node_type& v)const {
            return v.j;
        }

        size_type
        degree(const node_type &v) const {
            auto ch = suffix_children(v);
            return ch.size();
        }

        node_type
        parent(const node_type& v) {
            if (v == root(v.dir))
                return root(v.dir);
            else
                return root(v.dir); //not supported unless parent tables are stored;
        }


        node_type
        sl(const node_type& v)const {
            if (v == root(v.dir))
                return root(v.dir);
            return root(v.dir); //not supported unless slink tables are stored;
        }

        node_type
        child(const node_type &v, const char_type c, size_type &char_pos) const {
            if (is_leaf(v))
                return root(v.dir);
            size_type l = getlcp(v);
            size_type text_pos = SA[v.dir][v.i] + l;
            if (text_pos != Text.size() - 1) //always should be the first child if exist
                text_pos -= v.dir * (2 * text_pos - Text.size() + 2);
            node_type nc(v.i, 0, 0, v.k, v.dir);
            //get first child
            if (v.j < (Text.size() - 1) and (LCP[v.dir][v.j] > LCP[v.dir][v.j + 1])) {
                nc.j = childtab[v.dir][v.j] - 1;
                if (nc.j < nc.i)
                    nc.j = childtab[v.dir][nc.i] - 1;
            }
            else
                nc.j = childtab[v.dir][v.i] - 1;
            while (nc != root(v.dir) and (char_type) Text[text_pos] != c) {
                if (c < (char_type) Text[text_pos] or nc.j == (Text.size() - 1) or LCP[nc.dir][nc.i] > LCP[nc.dir][nc.j + 1])
                    nc = root(v.dir);  // v has not more children
                else { //get next sibling
                    nc.i = nc.j + 1;
                    if (childtab[nc.dir][nc.i] <= nc.i)
                        nc.j = v.j;
                    else {
                        if (LCP[nc.dir][nc.i] == LCP[nc.dir][childtab[nc.dir][nc.i]])
                            nc.j = childtab[nc.dir][nc.i] - 1; //using nextlindex
                        else
                            nc.j = v.j;
                    }
                }
                text_pos = SA[v.dir][nc.i] + l;
                text_pos -= v.dir * (2 * text_pos - Text.size() + 2);
            }
            char_pos = text_pos;  //note that this is the position in the original text (no reverse)
            nc.l = getlcp(nc);
            return nc;
        }

        //! Get the child w of node v which edge label (v,w) starts with character c.
        // \sa child(node_type v, const char_type c, size_type &char_pos)
        node_type
        child(const node_type &v, const char_type c) const {
            size_type char_pos;
            return child(v, c, char_pos);
        }


        /*!
         * Return the reverse of a node v assuming that v.k = 0 and
         * that the node is an explcit node. If  v is implicit move
         * to the closer explicit reverse node
         * */
        node_type
        get_reverse(const node_type &v) const {
            size_type n_home = home(v), length = v.j - v.i;
            node_type rev;
            if (is_leaf(v))
                rev = root(v.dir);
            else
                rev = node_type(affixtab[v.dir][n_home], affixtab[v.dir][n_home] + length,
                                    0, 0, !v.dir);
            rev.l = getlcp(rev);
            return rev;
        }

        //! Compute the suffix number of a leaf node v (or the value of the left most sa).
        size_type
        sn(node_type v) const {
            return SA[v.dir][v.i];
        }

        //!Return all possible node extension of one character for v, returning the resulting nodes
        std::vector<node_type>
        suffix_children(const node_type &v) const {
            std::vector<node_type> suff_ch;
            node_type c_node;
            size_type v_l = getlcp(v);
            if (v.l < getlcp(v)) {
                c_node = node_type(v.i, v.j, v_l + 1, v.k, v.dir);
                suff_ch.push_back(c_node);
            }
            else { // compute all the children of v
                c_node = node_type(v.i, 0, 0, 0, v.dir);
                if (v.j < (Text.size() - 1) and (LCP[v.dir][v.j] > LCP[v.dir][v.j + 1])) {
                    c_node.j = childtab[v.dir][v.j] - 1;
                    if (c_node.j < c_node.i)
                        c_node.j = childtab[v.dir][c_node.i] - 1;
                }
                else
                    c_node.j = childtab[v.dir][v.i] - 1;
                while (c_node != root(v.dir)) {
                    c_node.l += 1; // =getlcp(c_node);
                    c_node.k = v.k;
                    suff_ch.push_back(c_node);
                    //get next child
                    if (c_node.j == (Text.size() - 1) or
                            LCP[c_node.dir][c_node.i] > LCP[c_node.dir][c_node.j + 1])
                        c_node = root(v.dir);  // v has not more children
                    else { //get next sibling
                        c_node.i = c_node.j + 1;
                        if (childtab[c_node.dir][c_node.i] <= c_node.i)
                            c_node.j = v.j;
                        else {
                            if (LCP[c_node.dir][c_node.i] == LCP[c_node.dir][childtab[c_node.dir][c_node.i]])
                                c_node.j = childtab[c_node.dir][c_node.i] - 1; //using nextlindex
                            else
                                c_node.j = v.j;
                        }
                    }
                }
            }
            return suff_ch;
        }

        std::vector<node_type>
        prefix_children(const node_type &v) const {
            std::vector<node_type> pref_ch;
            node_type c_node;
            if (v.k > 0) {
                c_node = node_type(v.i, v.j, v.l, v.k - 1, v.dir);
                pref_ch.push_back(c_node);
            }
            else if(is_leaf(v)) {
                if (SA[v.dir][v.i] - v.n_l > 0)  // note that k == 0 if we got here
                    c_node = node_type(v.i, v.j, v.l, 0, v.dir, v.n_l + 1);
                else
                    return pref_ch; //not possible prefix child
                pref_ch.push_back(c_node);
            }
            else { //v is multiple with k == 0
                //move to the reverse
                c_node = get_reverse(v);
                c_node.l = getlcp(v);
                c_node.k = c_node.l - v.l;
                pref_ch = suffix_children(c_node);
            }
            return pref_ch;
        }

        /*!
		 * Search for pattern in backward direction (using the original text as reference),
         * starting with a given node
		 * \param[in] pattern string of character to search for
		 * \param[in,out] s node (will be overwritten)
		 * \return size of the new interval
		*/
        template<class t_pat_iter>
        size_t
        backward_search(t_pat_iter begin, t_pat_iter end, node_type &s) {
            if (s.dir == 1) {
              //  std::cout << "doing: fwd_search(end, begin, s, 1)" << std::endl;
                return fwd_search(end, begin, s, 1);
            }
            else {
             //   std::cout << "doing: bwd_search(end, begin, s)" << std::endl;
                return bwd_search(begin, end, s, 0);
            }
        }

        /*!
		 * Search for a char_type c in backward direction (using the original text as reference),
         * starting with a given node
		 * \param[in] character to search for
		 * \param[in,out] s node (will be overwritten)
		 * \return size of the new interval
		*/
        size_t
        backward_search(char_type c, node_type &s) {
            std::vector<char_type> v(1, c);
            return backward_search(v.begin(), v.end(), s);
        }

        /*!
		 * Search for pattern in forward direction (taking the original text as reference),
         * starting with a given node
		 * \param[in] pattern string of character to search for
		 * \param[in,out] s node (will be overwritten)
		 * \return size of the new interval
		 */
        template<class t_pat_iter>
        size_t
        forward_search(t_pat_iter begin, t_pat_iter end, node_type &s) {
            if(s.dir == 0) {
             //   std::cout << "doing: fwd_search(begin, end, s, 0)" << std::endl;
                return fwd_search(begin, end, s, 0);
            }
            else{
             //   std::cout << "doing: bwd_search(end, begin, s, 1)" << std::endl;
                return bwd_search(end, begin, s, 1);
            }
        }

        /*!
		 * Search for a char_type c in fordward direction (using the original text as reference),
         * starting with a given node
		 * \param[in] character to search for
		 * \param[in,out] s node (will be overwritten)
		 * \return size of the new interval
		*/
        size_t
        forward_search(char_type c, node_type &s) {
            std::vector<char_type> v(1, c);
            return forward_search(v.begin(), v.end(), s);
        }

        //! Return the local forward search (taking as reference the direction of the node s) of the string
        // from begin to end starting from the node s.
        template<class t_pat_iter>
        size_t
        fwd_search(t_pat_iter begin, t_pat_iter end, node_type &s, size_type it_dir) {
            char_type c;
            size_type text_pos = 0, max_length = 0;
            t_pat_iter it = begin;
            while (it != end) {
                text_pos = SA[s.dir][s.i] + s.l;
                if (text_pos != Text.size() - 1) //always should be the first child if exist
                    text_pos -= s.dir * (2 * text_pos - Text.size() + 2);
                if (is_leaf(s)) { //search the text from the leaf
                    max_length = Text.size() - SA[s.dir][s.i];
                    while (it != end and s.l < max_length) {
                        it -= it_dir;
                        c = Text[text_pos];
                        if (c != *it) {
                            s = root(s.dir);
                            return 0;
                        }
                        it += (1 - it_dir);
                        ++s.l;
                        text_pos += 1 - s.dir * 2;
                    }
                    if (s.l == max_length and it != end) {
                        s = root(s.dir);
                        return 0;
                    }
                }
                else { //s is a explicit or implicit node
                    max_length = getlcp(s);
                    while (it != end and max_length > s.l) { //we are inside of an implicit  node
                        it -= it_dir;
                        c = Text[text_pos];
                        if (c != *it) {
                            s = root(s.dir);
                            return 0;
                        }
                        it += (1 - it_dir);
                        ++s.l;
                        text_pos += 1 - s.dir * 2;
                    }
                    if(it != end) {
                        it -= it_dir;
                        c = *it;
                        s = child(s, c);
                        if(s == root(s.dir))
                            return 0;
                        s.l = max_length + 1;
                        it += (1 - it_dir);
                    }
                }
            }
            return size(s);
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
                text_pos = SA[s.dir][s.i] + s.k - 1;
                text_pos -= s.dir * (2 * text_pos - Text.size() + 2);
                while (begin != it and s.k > 0) {
                    //std::cout << "enter " << *it << std::endl;
                    it -= (1 - it_dir); //this may be wrong
                    c = Text[text_pos];
                    if (c != *it) {
                        s = root(s.dir);
                        return 0;
                    }
                    it += it_dir;
                    text_pos -= 1 + s.dir * 2;
                    --s.k;
                }
                if (it == begin)
                    return size(s);
            }
            //s.k == 0 in this point
            if (is_leaf(s)) {
                max_length = SA[s.dir][s.i] - s.n_l;
                text_pos = max_length - 1;
                text_pos -= s.dir * (2 * text_pos - Text.size() + 2);
                while (begin != it and max_length > 0) { //keep looking backward as far as possible
                    it -= (1 - it_dir);
                    c = Text[text_pos];
                    if (c != *it) {
                        s = root(s.dir);
                        return 0;
                    }
                    it += it_dir;
                    text_pos -= 1 + s.dir * 2;
                    ++s.n_l;
                    --max_length;
                }
                if (it == begin and max_length > 0)
                    return 1;
            }
            else { // move to the reverse node and do fwd_search
                r_node = get_reverse(s);
                r_node.l = getlcp(s);
                r_node.k = r_node.l - s.l;
                s = r_node;     //note that this can not be a leaf
                return fwd_search(it, begin, s, 1 - it_dir);
            }
            s = root(s.dir);
            return 0;
        }

        void
        print_node_path(const node_type &v) {
            std::string path = "";
            size_type text_pos = SA[v.dir][v.i] - v.n_l + v.k; //note that k and n_l are never at the same time >0
            size_type length_path = v.l - v.k + v.n_l;
            text_pos -= v.dir * (2 * text_pos - Text.size() + 2);
            for (size_type i = 0; i < length_path; ++i) {
                path += (char_type) Text[text_pos];
                text_pos+= 1 - v.dir * 2;
            }
            std::cout << path << std::endl;
        }

        size_type
        path_length(const node_type &v) {
            return v.l - v.k + v.n_l;
        }

    private:

        size_type
        home(const node_type &v) const {
            if (v == root(v.dir))
                return v.j;
            if (v.j == LCP[v.dir].size() - 1 or LCP[v.dir][v.i] >= LCP[v.dir][v.j + 1])
                return v.i;
            else
                return v.j;
        }

        size_type
        getlcp(const node_type &v) const {
            if (is_leaf(v))  //should not return
                return LCP[v.dir][v.i];
            if (v.i < childtab[v.dir][v.j] and childtab[v.dir][v.j] <= v.j)
                return LCP[v.dir][childtab[v.dir][v.j]];

            else
                return LCP[v.dir][childtab[v.dir][v.i]];
        }

        void
        create_structures(sdsl::cache_config &config, uint8_t index) {
            {
                sdsl::int_vector_buffer<> sa_buf(sdsl::cache_file_name(sdsl::conf::KEY_SA, config));
                SA[index] = sdsl::int_vector<t_width>(sa_buf.size(), 0, sa_buf.width());
                for (size_type i = 0; i < sa_buf.size(); ++i)
                    SA[index][i] = sa_buf[i];
            }
            {
                LCP[index] = sdsl::lcp_byte<>(config);
            }
            //create the extra structures
            compute_childtab(index);
            {
                auto event = sdsl::memory_monitor::event("SuffixLink tab");
                if (!cache_file_exists(KEY_SLi, config))
                    compute_suffixlinktab(config, index);
                register_cache_file(KEY_SLi, config);
                register_cache_file(KEY_SLj, config);
            }
        }

        void
        compute_childtab(uint8_t index) {
            childtab[index] = sdsl::int_vector<t_width>(SA[index].size(), 0, SA[index].width());
            uint64_t last_index = (uint64_t) -1;
            std::stack<uint64_t> prev_values;
            prev_values.push(0);
            for (uint64_t i = 1; i < LCP[index].size(); ++i) {
                while (LCP[index][i] < LCP[index][prev_values.top()]) {
                    last_index = prev_values.top();
                    prev_values.pop();
                    if (LCP[index][i] <= LCP[index][prev_values.top()] and
                        LCP[index][prev_values.top()] != LCP[index][last_index]) {
                        if (childtab[index][prev_values.top()] == 0) {
                            childtab[index][prev_values.top()] = last_index; //down[top]
                        }
                    }
                }
                if (last_index != (uint64_t) -1) {
                    childtab[index][i - 1] = last_index; //up[i]
                    last_index = (uint64_t) -1;
                }
                if (LCP[index][i] == LCP[index][prev_values.top()])
                    childtab[index][prev_values.top()] = i; //nextlindex[top]
                prev_values.push(i);
            }
            //special case for last child.down value
            while (!prev_values.empty()) {
                last_index = prev_values.top();
                prev_values.pop();
                if (!prev_values.empty()) {
                    if (0 <= LCP[index][prev_values.top()] and
                        LCP[index][prev_values.top()] != LCP[index][last_index]) {
                        if (childtab[index][prev_values.top()] == 0)
                            childtab[index][prev_values.top()] = last_index; //down[top]
                    }
                }
            }
        }


        void
        traverse_slink(node_type &v, sdsl::int_vector<> &sl_i, sdsl::int_vector<> &sl_j) {
            if (is_leaf(v))  // if v is a leave, v has no child
                return;
            node_type c_node(v.i, 0, 0, 0, v.dir);
            if (v.j < (Text.size() - 1) and (LCP[v.dir][v.j] > LCP[v.dir][v.j + 1])) {
                c_node.j = childtab[v.dir][v.j] - 1;
                if (c_node.j < c_node.i)
                    c_node.j = childtab[v.dir][c_node.i] - 1;
            }
            else
                c_node.j = childtab[v.dir][v.i] - 1;
            while (c_node != root(v.dir)) {
                c_node.l = getlcp(c_node);
                if (!is_leaf(c_node)) {
                    if (sl_i[home(c_node)] == Text.size())
                        compute_slink(c_node, v, sl_i, sl_j);
                    traverse_slink(c_node, sl_i, sl_j);
                }
                //get next child
                if (c_node.j == (Text.size() - 1) or
                    LCP[c_node.dir][c_node.i] > LCP[c_node.dir][c_node.j + 1])
                    c_node = root(v.dir);  // v has not more children
                else { //get next sibling
                    c_node.i = c_node.j + 1;
                    if (childtab[c_node.dir][c_node.i] <= c_node.i)
                        c_node.j = v.j;
                    else {
                        if (LCP[c_node.dir][c_node.i] == LCP[c_node.dir][childtab[c_node.dir][c_node.i]])
                            c_node.j = childtab[c_node.dir][c_node.i] - 1; //using nextlindex
                        else
                            c_node.j = v.j;
                    }
                }
            }//end while
        }

        void
        compute_slink(node_type &c_node, node_type &p_node, sdsl::int_vector<> &sl_i, sdsl::int_vector<> &sl_j) {
            size_type n_home = home(p_node), last_l = 0;
            size_type text_pos = SA[c_node.dir][c_node.i] + p_node.l;
            node_type slink(sl_i[n_home], sl_j[n_home], 0, 0, p_node.dir);
            slink.l = getlcp(slink);
            n_home = home(c_node);
            if (p_node == root(p_node.dir)) { //special case for the childs of the root
                if (c_node.l == 1)
                    slink = root(p_node.dir);
                else {
                    text_pos += 1;
                    text_pos -= p_node.dir * (2 * text_pos - Text.size() + 2);
                    slink = child(slink, (char_type) Text[text_pos]);
                }
            }
            else {
                text_pos -= p_node.dir * (2 * text_pos - Text.size() + 2);
                last_l = slink.l;
                slink = child(slink, (char_type) Text[text_pos]);
            }
            while (slink.l != c_node.l - 1) {
                text_pos += (slink.l - last_l) - p_node.dir * 2 * (slink.l - last_l);
                last_l = slink.l;
                slink = child(slink, (char_type) Text[text_pos]);
            }
            sl_i[n_home] = slink.i;
            sl_j[n_home] = slink.j;
        }

        void
        compute_suffixlinktab(sdsl::cache_config &config, size_type index) {
            size_type n = Text.size();
            sdsl::int_vector<> slink_i(n, n, SA[index].width());
            sdsl::int_vector<> slink_j(n, n, SA[index].width());
            slink_i[n - 1] = 0;
            slink_j[n - 1] = n - 1;
            node_type v = root(index);
            traverse_slink(v, slink_i, slink_j);
            std::cout << "Finished SLINK " << index << std::endl;
            store_to_cache(slink_i, KEY_SLi, config);
            store_to_cache(slink_j, KEY_SLj, config);
        }

        void
        create_affixlinktab(sdsl::cache_config &config, size_type index){
            size_type n= SA[index].size();
            {
                sdsl::int_vector_buffer<> sl_i(sdsl::cache_file_name(KEY_SLi, config));
                sdsl::int_vector_buffer<> sl_j(sdsl::cache_file_name(KEY_SLj, config));
                sdsl::int_vector<> sl_aux[2];   //suffixlink aux
                sl_aux[0] = sdsl::int_vector<t_width>(sl_i.size(), 0, sl_i.width());
                sl_aux[1] = sdsl::int_vector<t_width>(sl_j.size(), 0, sl_j.width());
                for (size_type i = 0; i < sl_i.size(); ++i) {
                    sl_aux[0][i] = sl_i[i];
                    sl_aux[1][i] = sl_j[i];
                }
                affixtab[index] = sdsl::int_vector<t_width>(n, n, SA[index].width());
                affixtab[index][n - 1] = 0; //set the root link
                node_type v = root(index);
                traverse_alink_louds(v, sl_aux[0], sl_aux[1]);
                //traverse_alink(v, sl_aux[0], sl_aux[1]);
            }
            std::cout << "Finished ALINK " << index << std::endl;
        }


         void
         traverse_alink(node_type &v, sdsl::int_vector<> &sl_i, sdsl::int_vector<> &sl_j) {
             if (is_leaf(v))  // if v is a leave, v has no child
                 return;
             node_type rev_node;
             node_type c_node(v.i, 0, 0, 0, v.dir);
             if (v.j < (Text.size() - 1) and (LCP[v.dir][v.j] > LCP[v.dir][v.j + 1])) {
                c_node.j = childtab[v.dir][v.j] - 1;
                if (c_node.j < c_node.i)
                    c_node.j = childtab[v.dir][c_node.i] - 1;
             }
             else
                 c_node.j = childtab[v.dir][v.i] - 1;
             while (c_node != root(v.dir)) {
                 c_node.l = getlcp(c_node);
                 if (!is_leaf(c_node)) {
                     if (affixtab[v.dir][home(c_node)] == Text.size()) {
                         rev_node = reverse_interval(c_node, sl_i, sl_j);
                         affixtab[v.dir][home(c_node)] = rev_node.i;
                     }
                     traverse_alink(c_node, sl_i, sl_j);
                 }
                 //get next child
                 if (c_node.j == (Text.size() - 1) or
                     LCP[c_node.dir][c_node.i] > LCP[c_node.dir][c_node.j + 1])
                     c_node = root(v.dir);  // v has not more children
                 else { //get next sibling
                     c_node.i = c_node.j + 1;
                     if (childtab[c_node.dir][c_node.i] <= c_node.i)
                         c_node.j = v.j;
                     else {
                         if (LCP[c_node.dir][c_node.i] == LCP[c_node.dir][childtab[c_node.dir][c_node.i]])
                             c_node.j = childtab[c_node.dir][c_node.i] - 1; //using nextlindex
                         else
                             c_node.j = v.j;
                     }
                 }
             }//end while
        }

        void
        traverse_alink_louds(node_type &v, sdsl::int_vector<> &sl_i, sdsl::int_vector<> &sl_j) {
            if (is_leaf(v))  // if v is a leave, v has no child
                return;
            std::pair<size_type, size_type> tree_node;
            std::queue<std::pair<size_type, size_type>> to_check;
            node_type rev_node, c_node;
            size_type max_rev = 10, last_l = 0;
            size_type n_text = Text.size(), i, j;
            to_check.push(std::make_pair(v.i, v.j));
            while (!to_check.empty()) {
                //check the node in the top of the queue
                tree_node = to_check.front();
                i = tree_node.first; j = tree_node.second;
                c_node = node_type(i, j, 0, 0, v.dir);
                c_node.l = getlcp(c_node);
                if(c_node.l > last_l + max_rev) { //is too big put at the again
                    to_check.pop();
                    to_check.push(tree_node);
                }
                else {
                    if (c_node.l > last_l)
                        last_l = c_node.l;
                   // std::cout << "checking: " << c_node << std::endl;
                    if (affixtab[v.dir][home(c_node)] == n_text) {
                        rev_node = reverse_interval(c_node, sl_i, sl_j);
                        affixtab[v.dir][home(c_node)] = rev_node.i;
                    }
                    //add its children that are not leaf to the queue
                    //get first child
                    c_node = node_type(i, 0, 0, 0, v.dir);
                    if (j < (n_text - 1) and (LCP[v.dir][j] > LCP[v.dir][j + 1])) {
                        c_node.j = childtab[v.dir][j] - 1;
                        if (c_node.j < c_node.i)
                            c_node.j = childtab[v.dir][c_node.i] - 1;
                    }
                    else
                        c_node.j = childtab[v.dir][i] - 1;
                    while (c_node != root(v.dir)) {
                        c_node.l = getlcp(c_node);
                        if (!is_leaf(c_node))
                            to_check.push(std::make_pair(c_node.i, c_node.j));
                        //get next child
                        if (c_node.j == (Text.size() - 1) or
                            LCP[c_node.dir][c_node.i] > LCP[c_node.dir][c_node.j + 1])
                            c_node = root(v.dir);  // v has not more children
                        else { //get next sibling
                            c_node.i = c_node.j + 1;
                            if (childtab[c_node.dir][c_node.i] <= c_node.i)
                                c_node.j = j;
                            else {
                                if (LCP[c_node.dir][c_node.i] == LCP[c_node.dir][childtab[c_node.dir][c_node.i]])
                                    c_node.j = childtab[c_node.dir][c_node.i] - 1; //using nextlindex
                                else
                                    c_node.j = j;
                            }
                        }
                    }//end while of children
                    to_check.pop();
                }
            }
        }

        node_type
        reverse_interval(node_type &v, sdsl::int_vector<> &sl_i, sdsl::int_vector<> &sl_j) {
            size_type n_home = home(v), length = v.j - v.i, text_pos;
            node_type rev, sl_node;
            if (affixtab[v.dir][n_home] != Text.size()) {
                rev = node_type(affixtab[v.dir][n_home], affixtab[v.dir][n_home] + length,
                                0,0, 1 - v.dir);
                rev.l = getlcp(rev);
            }
            else {
              //  std::cout << "aca" << std::endl;
                text_pos = SA[v.dir][v.i];
                if (text_pos != Text.size() - 1)
                    text_pos -= v.dir * (2 * text_pos - Text.size() + 2);
                sl_node = node_type(sl_i[n_home], sl_j[n_home], v.l -1, 0, v.dir);
                rev = reverse_interval(sl_node, sl_i, sl_j);
                if (getlcp(rev) <= getlcp(sl_node))
                    rev = child(rev, Text[text_pos]);
                affixtab[v.dir][n_home] = rev.i;
            }
            return rev;
        }

    }; //end class AFA

    template<class t_alphabet_strat, uint8_t t_width>
    AFA<t_alphabet_strat, t_width>::AFA(sdsl::cache_config &config, sdsl::cache_config &rev_config) {
        {
            // (1) get text
            std::string file = sdsl::cache_file_name(sdsl::key_trait<alphabet_type::int_width>::KEY_TEXT, config);
            sdsl::int_vector_buffer<alphabet_type::int_width> text_buf(file);
            Text = sdsl::int_vector<t_width>(text_buf.size(), 0, text_buf.width());
            std::cout << "Text length: " << text_buf.size() << std::endl;
            for (size_type i = 0; i < text_buf.size(); ++i)
                Text[i] = text_buf[i];
            alphabet_type tmp_alphabet(text_buf, text_buf.size());
            m_alphabet.swap(tmp_alphabet);
        }
        create_structures(config, 0);
        create_structures(rev_config, 1);
        create_affixlinktab(config, 0);
        create_affixlinktab(rev_config, 1);
    }


    template<class t_int>
    inline std::ostream& operator<<(std::ostream& os, const interval_sa<t_int>& interval) {
        os << "-[" << interval.i << "," << interval.j << "] l: " << interval.l <<
        " k: " << interval.k << " n_l: " << interval.n_l << " dir: " << interval.dir;
        return os;
    }


} //end namespace
#endif //AFFIX_ARRAY
