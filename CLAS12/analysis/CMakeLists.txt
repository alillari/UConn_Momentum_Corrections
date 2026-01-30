file(GLOB ANALYSIS_SOURCES "*.cpp")

foreach(src ${ANALYSIS_SOURCES})
    get_filename_component(exec_name ${src} NAME_WE)

    add_executable(${exec_name} ${src})

    # Link both your library and ROOT explicitly
    target_link_libraries(${exec_name} PRIVATE
        momcorr
	fitting
        ROOT::Core
        ROOT::RIO
        ROOT::Tree
        ROOT::Thread
    	ROOT::ROOTDataFrame
    )

    set_target_properties(${exec_name} PROPERTIES
        RUNTIME_OUTPUT_DIRECTORY
        ${PROJECT_SOURCE_DIR}/CLAS12/executable
    )
endforeach()
