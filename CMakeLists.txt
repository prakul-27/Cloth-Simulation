cmake_minimum_required(VERSION 3.0)

project(minor)
set(CMAKE_CXX_STANDARD 11)

find_package(GLEW REQUIRED)
find_package(glm REQUIRED)
find_package(OpenGL REQUIRED)

add_library(minor Cloth.cpp Particle.cpp Spring.cpp)
target_link_libraries(minor GLEW::GLEW glm::glm GL GLU glut freeimage )

add_executable(Main Main.cpp)
target_link_libraries(Main minor)

