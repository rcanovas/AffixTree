cmake_minimum_required(VERSION 3.3)
project(AffixTree)

#unset if you want to do debuging 
set(CMAKE_BUILD_TYPE "Release")

set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -std=c++11")

#modify these directories to include sdsl library location
include_directories(/home/your_user/include)  
link_directories(/home/your_user/lib)
link_libraries(sdsl divsufsort divsufsort64)

set(SOURCE_FILES test/createIndex.cpp include/affix_tree.hpp)
set(SOURCE_FILES2 test/testSearch.cpp include/affix_tree.hpp)

add_executable(createIndex ${SOURCE_FILES})
add_executable(testSearch ${SOURCE_FILES2})

