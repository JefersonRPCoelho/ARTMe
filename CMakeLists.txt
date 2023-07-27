cmake_minimum_required(VERSION 2.8)
project (ARTMe)

#Bring the headers into the project.
include_directories(include)

set(EXECUTABLE_OUTPUT_PATH "/bin/")
set(CMAKE_INSTALL_PREFIX "./")

#Adding sources.
file(GLOB SOURCES "src/*.cpp")


#Generate the static library from the sources
add_library(ARTMeStatic STATIC ${SOURCES})

add_executable (ARTMe ${SOURCES})

#Generate the shared library from the sources
add_library(ARTMeShared SHARED ${SOURCES})

set_target_properties(ARTMeStatic PROPERTIES OUTPUT_NAME ARTMe)
#set_target_properties(ARTMeShared PROPERTIES OUTPUT_NAME ARTMe)

#install(TARGETS ARTMeStatic DESTINATION lib)
install(TARGETS ARTMeShared DESTINATION lib)
install(TARGETS ARTMe DESTINATION bin)

# Compile features: requires c++11
target_compile_features(ARTMe PUBLIC cxx_std_11)
# Avoid extensions such as "gnu++11" instead of "c++11", which could include non-standard features
set_target_properties(ARTMe PROPERTIES CXX_EXTENSIONS OFF)

