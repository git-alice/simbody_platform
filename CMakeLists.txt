cmake_minimum_required(VERSION 3.19)
project(simbody)

set(CMAKE_CXX_STANDARD 14)

#add_executable(simbody long_pendulum.cpp two_wheels.cpp motor_torque_limit.cpp)
#add_executable(simbody two_wheels_2.cpp)
add_executable(simbody one_wheel_2.cpp)
#add_executable(simbody motor_torque_limit.cpp)

include_directories(simbody "/usr/local/Cellar/simbody/3.7/include/simbody")

target_link_libraries(simbody "/usr/local/Cellar/simbody/3.7/lib/libSimTKcommon.dylib")
target_link_libraries(simbody "/usr/local/Cellar/simbody/3.7/lib/libSimTKmath.dylib")
target_link_libraries(simbody "/usr/local/Cellar/simbody/3.7/lib/libSimTKsimbody.dylib")