if(ENABLE_MECHANICS)
    create_symlink("meshes")

    function(add_custom_application source)
        set(parameters ${ARGN})

        add_executable(${source} ${source}.cpp)
        target_link_libraries(${source} NuToMechanics NuToMath NuToBase
            -lboost_program_options ${Boost_LIBRARIES}
            ${LAPACK_LIBRARIES} ${ANN_LIBRARIES} ${MUMPS_LIBRARIES})
        if(ENABLE_VISUALIZE)
            target_link_libraries(${source} NuToVisualize)
        endif()

        add_test(custom::${source} ${source} ${parameters})
    endfunction()

    add_subdirectory(damage)
    add_subdirectory(heatconduction)
    add_subdirectory(thermoplastic)
    #add_subdirectory(thermoelastic)
    #add_subdirectory(elastic)
    #add_subdirectory(loadtest)
    add_subdirectory(meshCreation)

endif()
