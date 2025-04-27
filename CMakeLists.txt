project(LOAD_HIGHS LANGUAGES CXX)

set(HIGHS_DIR path_to_highs_install/lib/cmake/highs)

find_package(HIGHS REQUIRED)
find_package(Threads REQUIRED)

add_executable(main main.cpp)
target_link_libraries(main highs::highs)

