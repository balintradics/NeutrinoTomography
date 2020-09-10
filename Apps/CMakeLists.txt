set(APPS
         test_neutrinoosc
	 Nu_EarthCell
	 Run_NuFlux_Simple
	 Run_NuFlux_Osc
	 Run_NuFlux_OscLongitude)


foreach(a ${APPS})
          add_executable(${a} ${a}.cpp)

          add_dependencies(${a} NeutrinoOsc EarthModel)
          target_include_directories(${a} PRIVATE ${CMAKE_SOURCE_DIR})
          target_link_libraries(${a} -L${CMAKE_BINARY_DIR}/lib NeutrinoOsc EarthModel )

          install(TARGETS ${a} DESTINATION bin)
endforeach()
