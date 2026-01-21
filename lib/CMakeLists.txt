# Collect all .cc files in lib/
file(GLOB LIB_SOURCES "*.cpp")

add_library(momcorr ${LIB_SOURCES})

target_include_directories(momcorr PUBLIC
    ${PROJECT_SOURCE_DIR}/include
)

target_link_libraries(momcorr PUBLIC
    ROOT::Core
    ROOT::Tree
    ROOT::RDataFrame
)
