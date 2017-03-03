if(ENABLE_MECHANICS)
    execute_process(COMMAND "${CMAKE_COMMAND}" "-E" "create_symlink"
        "${CMAKE_CURRENT_SOURCE_DIR}/meshes" "${CMAKE_CURRENT_BINARY_DIR}/meshes")

    function(add_custom_application source)
        set(parameters ${ARGN})

        add_executable(${source} ${source}.cpp)
        target_link_libraries(${source} NuToMechanics NuToMath NuToBase
            -lboost_program_options -lprofiler ${Boost_LIBRARIES}
            ${LAPACK_LIBRARIES} ${ANN_LIBRARIES} ${MUMPS_LIBRARIES})
        if(ENABLE_VISUALIZE)
            target_link_libraries(${source} NuToVisualize)
        endif()

        add_test(custom::${source} ${source} ${parameters})
    endfunction()

    add_subdirectory(damage)
    add_subdirectory(heatconduction)
    add_subdirectory(thermoelastic)
    add_subdirectory(thermoplastic)

endif()
