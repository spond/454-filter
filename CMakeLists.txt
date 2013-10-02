cmake_minimum_required(VERSION 2.8)
project(qfilt)

# set(CMAKE_CXX_COMPILER "g++-4.7")
# set(CMAKE_EXE_LINKER_FLAGS "-g -pg -Wall")
set(INSTALL_PREFIX /usr/local CACHE PATH "Installation prefix")
set(CMAKE_INSTALL_PREFIX ${INSTALL_PREFIX} CACHE INTERNAL "Installation prefix" FORCE)

set(CMAKE_CONFIGURATION_TYPES Release)

include_directories(
    src/
)

add_executable(
    qfilt
    src/argparse.cpp
    src/ifile.cpp
    src/main.cpp
    src/seq.cpp
    src/strtok.cpp
)

target_link_libraries(qfilt m)

# do not remove -Wall and -Werror: please fix the errors instead of being lazy
set_target_properties(
    qfilt
    PROPERTIES
    COMPILE_FLAGS "-O3 -std=c++98 -pedantic-errors -Wall -Werror"
    # COMPILE_FLAGS "-O3 -g -pg -Wall -Werror"
)

install(
    TARGETS qfilt
	RUNTIME DESTINATION bin
	OPTIONAL
)

find_program(
    ASTYLE
    astyle
    PATHS ENV PATH
)

if(ASTYLE)
    add_custom_target(
        style
        COMMAND ${ASTYLE} --options=astyle.txt --suffix=none --recursive 'src/*.cpp' 'src/*.hpp' 'src/*.tpp'
        DEPENDS astyle.txt
    )
endif(ASTYLE)
