cmake_minimum_required(VERSION 3.16)

project(MomentumCorrections LANGUAGES CXX)

set(CMAKE_CXX_STANDARD 17)
set(CMAKE_CXX_STANDARD_REQUIRED ON)

# Find ROOT from your environment (conda, cvmfs, etc.)
find_package(ROOT REQUIRED COMPONENTS RIO Tree Thread Core ROOTDataFrame)

# Make headers visible everywhere
include_directories(
    ${PROJECT_SOURCE_DIR}/include
    ${ROOT_INCLUDE_DIRS}
)

# Build subdirectories
add_subdirectory(lib)
add_subdirectory(CLAS12)
