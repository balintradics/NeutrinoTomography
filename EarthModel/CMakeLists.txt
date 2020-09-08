set(DE_implementation_files
        DiscreteEarth.cpp)
        

set(DE_header_files
        DiscreteEarth.h)

add_library(EarthModel SHARED ${DE_implementation_files})
target_include_directories(EarthModel PRIVATE ${CMAKE_SOURCE_DIR})
target_link_libraries(EarthModel)

install(TARGETS EarthModel DESTINATION lib)
