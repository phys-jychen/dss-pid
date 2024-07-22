#----------------------------------------
# Project: dss-pid
#----------------------------------------

# Minimum version of CMake
cmake_minimum_required(VERSION 3.16)
#set(CMAKE_CXX_STANDARD 17)

# Name and version of the project
project(dss_pid VERSION 1.1.0 LANGUAGES C CXX)

# Type of CMake
set(CMAKE_BUILD_TYPE Debug)

# External Packages
find_package(ROOT REQUIRED COMPONENTS ROOTDataFrame TMVA ROOTVecOps)

# Set runtime output directory as bin
set(CMAKE_RUNTIME_OUTPUT_DIRECTORY ${CMAKE_CURRENT_SOURCE_DIR}/bin)

# Add the directories
include_directories(./include)
include_directories(./src)
include(${ROOT_USE_FILE})

# Add executables
add_executable(iPID pid.cxx src/BDT.cxx
                            src/Fit3D.cxx
                            src/Hits.cxx
                            src/Variables.cxx
)

# Link libraries
target_link_libraries(iPID ${ROOT_LIBRARIES})

# Add commands to set up the environment with the help of setup.sh...
execute_process(COMMAND cp ${CMAKE_CURRENT_SOURCE_DIR}/config/setup.sh ${PROJECT_BINARY_DIR})
execute_process(COMMAND sed -i "s:<PROJECT>:${CMAKE_CURRENT_SOURCE_DIR}:g" ${PROJECT_BINARY_DIR}/setup.sh)
execute_process(COMMAND chmod 774 ${PROJECT_BINARY_DIR}/setup.sh)
