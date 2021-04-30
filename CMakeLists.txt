cmake_minimum_required(VERSION 3.17)
project(SIFT)

set(CMAKE_CXX_STANDARD 11)

# 指定lib目录
link_directories(${PROJECT_SOURCE_DIR}/lib)
# 指定头文件搜索路径
include_directories(${PROJECT_SOURCE_DIR}/include)

set(CMAKE_EXE_LINKER_FLAGS "-static-libgcc -static-libstdc++")
add_executable(SIFT main.cpp SIFT.cpp include/SIFT.h)
# 将库链接到项目中
target_link_libraries(${PROJECT_NAME} libspot.a)

