if(ENABLE_MECHANICS)
    execute_process(COMMAND "${CMAKE_COMMAND}" "-E" "create_symlink"
        "${CMAKE_CURRENT_SOURCE_DIR}/meshes" "${CMAKE_CURRENT_BINARY_DIR}/meshes")

    foreach(test 
            DamageShell
            DamageShellNonlinear
            DamageShellMeso
            khoury_cooling
            khoury_damage
            TrussHeatCoupledRestrained
            TemperatureXdTransient
            hancock1
            damage_bar
            )

        add_executable(${test} ${test}.cpp)
        target_link_libraries(${test} NuToMechanics NuToMath NuToBase -lboost_program_options -lprofiler ${Boost_LIBRARIES} ${LAPACK_LIBRARIES} ${ANN_LIBRARIES})
        if(ENABLE_VISUALIZE)
            target_link_libraries(${test} NuToVisualize)
        endif()
        if(MUMPS_FOUND)
            target_link_libraries(${test} ${MUMPS_LIBRARIES} ${CMAKE_THREAD_LIBS_INIT})
        endif()

    endforeach()
endif()
