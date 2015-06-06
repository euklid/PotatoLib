project(PotatoLib)
cmake_minimum_required(VERSION 2.8)
#aux_source_directory(. SRC_LIST)
set(SOURCES
            cell.cpp
            element.cpp
            kernel.cpp
            main.cpp
            tree.cpp
            point.cpp
            fmm.cpp
)

set(HEADERS
            cell.h
            element.h
            kernel.h
            tree.h
            point.h
            fmm.h
)

#add_executable(${PROJECT_NAME} ${SRC_LIST})
add_executable(${PROJECT_NAME} ${SOURCES} ${HEADERS})
