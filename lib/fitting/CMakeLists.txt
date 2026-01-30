add_library(fitting
    FitModel.cpp
    PeakFitter.cpp
)

target_include_directories(fitting PUBLIC ${CMAKE_CURRENT_SOURCE_DIR})
target_link_libraries(fitting PUBLIC ROOT::Core ROOT::Hist ROOT::MathCore)
