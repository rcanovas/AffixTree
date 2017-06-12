## atds

The Affix Tree Data Structure (or atds) implementation contains the C++11 
codes associated with the following reference: 
R. Canovas and E. Rivals. "Full Compressed Affix Tree Representations". To 
appear in Proc. DCC'17. 
This work includes the following data structures:
- BidWT: The bidirectional wavelet index. The code presented is simply a 
	 mask over the original code of Thomas Schnattinger supplied within 
	 the sdsl library. 
- AFA:  Our implementation of Strothmann's Affix Array.
- ACAT: The Asynchronous Compressed Affix Array. 
- ACATS: The Asynchronous Compressed Affix Array Sampled.
- ACATN: The Asynchronous Compressed Affix Array Non-sampled.
- RACATN: The Reduced Asynchronous Compressed Affix Array Non-sampled.
- SCAT: The Synchronous Compressed Affix Array.

Each of these data structures receive as template parameter, the suffix tree or 
wavelet tree used internally. The createIndex test shows an example of how to 
create the data structure using the default templates. 

## Compile

To be able to compile the atds codes: 
- Install sdsl-lite. Follow the installation guide here: (https://github.com/simongog/sdsl-lite)
- Modify the location of the sdsl library in the CMakeLists.txt if necessary.
- Go to the build folder and run: 
	- cmake ..
	- make


## Methods

-[createIndex]:

	Use: ./createIndex <file_name> <tmp_locationopt> index_id
      		<file_name>: Name of the file to be use to create the required data structure 
          	<tmp_location>: folder where temporal files will be generated during the construction
		index_id: 	0 | BidWT
 				1 | AFA
 				2 | ACAT
 				3 | ACATS
 				4 | ACATN
				5 | RACATN
 				6 | SCAT
          	output:  <file_name>.<index_type>

		Example: ./createIndex ./data/file ./tmp 2
		output:  file.acat
        

-[testSearch]:

	Use: ./testSearch <index_file_name> index_id  <sample_file>
		<file_name>: Name of the compressed index file to be tested 
		index_id: 	0 | BidWT
 				1 | AFA
 				2 | ACAT
 				3 | ACATS
 				4 | ACATN
				5 | RACATN
 				6 | SCAT
		<sample_file>: A file containing a set of sampled reads to be queried , separed by '\0' 
                              (see the GetSamplePhrases.py code as example of how to generate these files)
			  
		output:  Times per operation

		Example: ./testSearch ./data/file.acat 2 ./sample/file_sample  
		output:  
		        Time avg per suffix_children operation: 1000000 us
			Time avg per prefix children operation: 1000000 us
			Time avg per Degree operation: 1000000 us
			Time avg per Slink operation: 1000000 us
			Time avg per backward forward: 1000000 us
			
			
			
Note: These codes assume that the computer have enough RAM memory to read and store the complete input.
