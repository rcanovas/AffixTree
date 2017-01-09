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

#include <iostream>
#include "./../include/affix_tree.hpp"

using namespace std;

template<class idx_type>
void
create_index(string file, string tmp_dir, string out_file) {
	idx_type idx;
	string id = sdsl::util::basename(file);
	sdsl::cache_config config(true, tmp_dir, id); //true -->erase tmp
	construct(idx, file, config, 1);
	ofstream out(out_file);
	std::cout << "Size Text: " << idx.size() << std::endl;
	std::cout << "Size in bytes: " << sdsl::size_in_bytes(idx) << " bytes" << std::endl;
	std::cout << "Size in bytes: " << (sdsl::size_in_bytes(idx) * 1.0 / idx.size()) << "n bytes" << std::endl;
	idx.serialize(out);
	out.close();
}

int main(int argc, char* argv[]) {
	if (argc != 4) {
		cout << "Usage: " << argv[0] << " file_name  tmp_location index_type" << endl;
		cout << " # | Index_type" << endl;
		cout << "---+--------------------" << endl;
		cout << " 0 | BidWT" << endl;
		cout << " 1 | AFA" << endl;
		cout << " 2 | ACAT" << endl;
		cout << " 3 | ACATS" << endl; 
		cout << " 4 | ACATN" << endl;
		cout << " 5 | RACATN" << endl;
		cout << " 6 | SCAT" << endl;
		return 1;
	}

	string file = argv[1];
	string out_file = file;
	string tmp_dir = argv[2];
	uint32_t index_type = (uint32_t)atoi(argv[3]);

	using timer = std::chrono::high_resolution_clock;
	auto start = timer::now();
	
	//create index
	switch (index_type) {
		case 0:
			out_file += ".bidwt";
			std::cout << "index: bidwt" << std::endl;
			create_index<atds::BidWT<> >(file, tmp_dir, out_file);
			break;
		case 1:
			out_file += ".afa";
			std::cout << "index: afa" << std::endl;
			create_index<atds::AFA<> >(file, tmp_dir, out_file);
			break;
		case 2:
			out_file += ".acat";
			std::cout << "index: acat" << std::endl;
			create_index<atds::ACAT<> >(file, tmp_dir, out_file);
			break;
		case 3:
			out_file += ".acats";
			std::cout << "index: ACATS" << std::endl;
			create_index<atds::ACATS<> >(file, tmp_dir, out_file);
			break;
		case 4:
			out_file += ".acatn";
			std::cout << "index: ACATN" << std::endl;
			create_index<atds::ACATN<> >(file, tmp_dir, out_file);
			break;
		case 5:
			out_file += ".racatn";
			std::cout << "index: RACATN" << std::endl;
			create_index<atds::RACATN<> >(file, tmp_dir, out_file);
			break;
		case 6:
			out_file += ".scat";
			std::cout << "index: scat" << std::endl;
			create_index<atds::SCAT<> >(file, tmp_dir, out_file);
			break;
		default:
			cout << "index_type must be a value in [0,6]" << endl;
			break;
	}
	std::cout << "The index was created" << std::endl;
  auto stop = timer::now();
	auto elapsed = stop-start;
	std::cout << "Time = " << std::chrono::duration_cast<std::chrono::milliseconds>(elapsed).count() << " ms"<< std::endl;
	return 0;
}

