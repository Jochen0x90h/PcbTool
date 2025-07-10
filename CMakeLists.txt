cmake_minimum_required(VERSION 3.25)

# set project
project(pcb-tool)

# configure languages
enable_language(CXX)
set(CMAKE_CXX_STANDARD 20)

# enable address sanitizer
#set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -fsanitize=address")
#set(CMAKE_EXE_LINKER_FLAGS "${CMAKE_EXE_LINKER_FLAGS} -fsanitize=address")


# dependencies
#find_package(OpenXLSX CONFIG)
find_package(xlnt CONFIG)

add_subdirectory(src)
