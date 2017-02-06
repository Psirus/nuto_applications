if(ENABLE_MECHANICS)
    foreach(meshfile
            Temperature2DMeso.msh
            Temperature2DHomogeneous.msh
            OneStone.msh
            TwoElements.msh
            cylinder.msh
            cylinder.geo
            Coarse2D.msh
            )

        configure_file(${CMAKE_CURRENT_SOURCE_DIR}/${meshfile} ${CMAKE_CURRENT_BINARY_DIR}/${meshfile} COPYONLY)
    endforeach()

    foreach(test 
            SlabTransient
            Coupled2DMesoTransient
            Coupled2DDamageTransient
            DamageShell
            DamageShellNonlinear
            DamageShellMeso
            khoury_cooling
            khoury_damage
            TrussHeatCoupled
            TrussHeatCoupledRestrained
            TemperatureXdTransient
            ThermoDamage
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
