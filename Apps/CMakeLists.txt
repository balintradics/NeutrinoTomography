set(APPS
         test_neutrinoosc
	 Nu_EarthCell
	 Run_NuFlux_Simple
	 Run_NuFlux_SimpleGlobe
	 Run_NuFlux_SimpleLongitude
	 Run_NuFlux_Osc
	 Run_NuFlux_OscLongitude
	 Run_NuProb_OscLongitude
	 RotateEarth
	 Print_SurfaceCells_Long
	 Run_CalcOscProb_Cell)


foreach(a ${APPS})
          add_executable(${a} ${a}.cpp)

          add_dependencies(${a} NeutrinoOsc EarthModel)
          target_include_directories(${a} PRIVATE ${CMAKE_SOURCE_DIR})
          target_link_libraries(${a} -L${CMAKE_BINARY_DIR}/lib NeutrinoOsc EarthModel )

          install(TARGETS ${a} DESTINATION bin)
endforeach()
