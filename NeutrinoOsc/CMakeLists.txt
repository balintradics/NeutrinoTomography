set(NO_implementation_files
	complex.cpp
	neutrino_osc.cpp  
	nrutil.cpp
	odeint.cpp
	thematrix.cpp )
        

set(NO_header_files
	complex.h  
	neutrino_osc.h  
	nrutil.h  
	odeint.h  
	thematrix.h )


add_library(NeutrinoOsc SHARED ${NO_implementation_files})
target_include_directories(NeutrinoOsc PRIVATE ${CMAKE_SOURCE_DIR})
target_link_libraries(NeutrinoOsc)
add_dependencies(NeutrinoOsc EarthModel)

install(TARGETS NeutrinoOsc DESTINATION lib)