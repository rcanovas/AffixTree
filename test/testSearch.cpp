/* atds - affix tree data structures
 * Copyright (C)2015-2016 Rodrigo Canovas
 * This program is free software: you can redistribute it and/or modify
 * it under the terms of the GNU General Public License as published by
 * the Free Software Foundation, either version 3 of the License, or
 * (at your option) any later version.
 *
 * This program is distributed in the hope that it will be useful,
 * but WITHOUT ANY WARRANTY; without even the implied warranty of
 * MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
 * GNU General Public License for more details.
 *
 * You should have received a copy of the GNU General Public License
 * along with this program.  If not, see http://www.gnu.org/licenses/ .
 * */

// Created by Rodrigo Canovas on 2016.
// Some of the tests are based on the implementations of
// Thomas Schnattinger tests for the bidirectional
// wavelet tree.

#include <iostream>
#include "./../include/affix_tree.hpp"

using namespace std;

template<class idx_type>
void
search_index(string file, string sample_file) {
    typedef typename idx_type::node_type node;
    typedef typename idx_type::size_type size_type;
    typedef typename idx_type::char_type char_type;
    //load index
    idx_type idx;
    sdsl::load_from_file(idx, file);
    //load sample
    std::vector<string> queries;
    sdsl::int_vector<> text_sample;
    sdsl::load_vector_from_file(text_sample, sample_file, 1);
    size_type n = text_sample.size();
    size_type word_len = 0;
    size_type pos = 0;
    string word = "";
    while (pos < n) {
        if (text_sample[pos] == '\0') {
            queries.push_back(word);
            if (word_len != word.size()){
                cout << "word size: " << word.size() << endl;
                word_len = word.size();
            }
            word = "";
        }
        else
            word += (char)text_sample[pos];
        pos ++;
    }
    cout << "size of sample: " << queries.size() << endl;
    //run test
    using timer = std::chrono::high_resolution_clock;
    node s = idx.root();
    size_type fo, ba, len;
    double avg = 0;
    int i = 0;
    for (string path: queries){
        //cout << i << endl;
        //cout << path << endl;
        i++;
        word_len = path.size();
        if (word_len%2 != 0)
            continue;
        fo = word_len / 2;
        ba = fo - 1;
        len = 0;
        //start taking time  Take word_len/2 operations each dir
        auto start = timer::now();
        s= idx.root();
        while (fo < word_len) {
            len = idx.forward_search((char_type)path[fo], s);
            len = idx.backward_search((char_type)path[ba], s);
            ++fo;
            --ba;
        }
        auto stop = timer::now();
        auto elapsed = stop - start;
        avg += (double)((chrono::duration_cast<chrono::microseconds>(elapsed).count() * 1.0) / word_len);
    }
    cout << "Time avg per backward forward: " << avg / queries.size() << " us" << endl;
}

template<class idx_type>
void
testOperation(string file, string sample_file) {
    typedef typename idx_type::node_type node;
    typedef typename idx_type::size_type size_type;
    typedef typename idx_type::char_type char_type;
    //load index
    idx_type idx;
    sdsl::load_from_file(idx, file);
    //load sample
    std::vector<string> queries;
    sdsl::int_vector<> text_sample;
    sdsl::load_vector_from_file(text_sample, sample_file, 1);
    size_type n = text_sample.size();
    size_type word_len = 0;
    size_type pos = 0;
    string word = "";
    while (pos < n) {
        if (text_sample[pos] == '\0') {
            queries.push_back(word);
            if (word_len != word.size()){
                cout << "word size: " << word.size() << endl;
                word_len = word.size();
            }
            word = "";
        }
        else
            word += (char)text_sample[pos];
        pos ++;
    }
    cout << "size of sample: " << queries.size() << endl;
    //run test
    using timer = std::chrono::high_resolution_clock;
    node s = idx.root();
    size_type fo = 0, n_q = 0, len = 0;
    double avg = 0, avg2 = 0, avg3 = 0, avg_sl = 0;
    int i = 0;
    for (string path: queries){
        word_len = path.size();
        if (word_len%2 != 0)
            continue;
        fo = 0;
        //start taking time  Take word_len/2 operations each dir
        auto start = timer::now();
        s= idx.root();
        do {
            n_q++;
            auto start = timer::now();
            auto ch = idx.suffix_children(s); //computes s children
            auto stop = timer::now();
            auto elapsed = stop - start;
            avg += (double)((chrono::duration_cast<chrono::microseconds>(elapsed).count() * 1.0));
            //cout << "suff: " << ch.size();

            start = timer::now();
            auto d = idx.degree(s);  //computes degree
            stop = timer::now();
            elapsed = stop - start;
            avg3 += (double)((chrono::duration_cast<chrono::microseconds>(elapsed).count() * 1.0));
            //cout << "  degree: " << d;

            start = timer::now();
            auto node_sl = idx.sl(s);  //computes Suffix Link
            stop = timer::now();
            elapsed = stop - start;
            avg_sl += (double)((chrono::duration_cast<chrono::microseconds>(elapsed).count() * 1.0));
            //cout << "  : " << d;

            start = timer::now();
            ch = idx.prefix_children(s);  //computes p// children
            stop = timer::now();
            elapsed = stop - start;
            avg2 += (double)((chrono::duration_cast<chrono::microseconds>(elapsed).count() * 1.0));
            //cout << "  pref: " << ch .size()<< endl;

            len = idx.forward_search((char_type)path[fo], s);
            ++fo;
        } while (fo < word_len and s != idx.root() and !idx.is_leaf(s));
    }
    cout << "Time avg per suffix_children operation: " << avg / n_q << " us" << endl;
    cout << "Time avg per prefix children operation: " << avg2 / n_q << " us" << endl;
    cout << "Time avg per Degree operation: " << avg3 / n_q << " us" << endl;
    cout << "Time avg per Slink operation: " << avg_sl / n_q << " us" << endl;

}


int main(int argc, char* argv[]) {
    if(argc < 4) {
        cout << "Usage: " << argv[0] << " index_file_name index_id  sample_file" << endl;
        cout << "index_id|" << endl;
				cout << " 0 | BidWT" << endl;
        cout << " 1 | AFA" << endl;
        cout << " 2 | ACAT" << endl;
        cout << " 3 | ACATS" << endl;
				cout << " 4 | ACATN" << endl;
				cout << " 5 | RACATN" << endl;
				cout << " 6 | SCAT" << endl;
        return 0;
    }

    string file = argv[1];
    string sample = argv[3];
    int index_id = atoi(argv[2]);

    //load index
    switch (index_id) {
        case 0:
            testOperation<atds::BidWT<> >(file, sample);
            search_index<atds::BidWT<> >(file, sample);
            break;
        case 1:
            testOperation<atds::AFA<> >(file, sample);
            search_index<atds::AFA<> >(file, sample);
            break;
        case 2:
            testOperation<atds::ACAT<> >(file, sample);
            search_index<atds::ACAT<> >(file, sample);
            break;
        case 3:
            testOperation<atds::ACATS<> >(file, sample);
            search_index<atds::ACATS<> >(file, sample);
            break;
         case 4:
            testOperation<atds::ACATN<> >(file, sample);
            search_index<atds::ACATN<> >(file, sample);
            break;
         case 5:
            testOperation<atds::RACATN<> >(file, sample);
            search_index<atds::RACATN<> >(file, sample);
            break;
        case 6:
            testOperation<atds::SCAT<> >(file, sample);
            search_index<atds::SCAT<> >(file, sample);
            break;
        default:
            cout << "index_id must be a value in [0,6]" << endl;
            break;
    }

    return 0;
}
