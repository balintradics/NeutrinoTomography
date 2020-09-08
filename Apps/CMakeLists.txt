set(APPS
         testme
         test_dearth
	 Nu_EarthCell)


foreach(a ${APPS})
          add_executable(${a} ${a}.cpp)

          add_dependencies(${a} NeutrinoOsc EarthModel)
          target_include_directories(${a} PRIVATE ${CMAKE_SOURCE_DIR})
          target_link_libraries(${a} -L${CMAKE_BINARY_DIR}/lib NeutrinoOsc EarthModel )

          install(TARGETS ${a} DESTINATION bin)
endforeach()
