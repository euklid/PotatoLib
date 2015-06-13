project(PotatoLib)
cmake_minimum_required(VERSION 2.8)
#aux_source_directory(. SRC_LIST)
set(SOURCES
            cell.cpp
            kernel.cpp
            main.cpp
            tree.cpp
            point.cpp
            fmm.cpp
            fmm2d.cpp
)

set(HEADERS
            cell.h
            element.h
            kernel.h
            tree.h
            point.h
            fmm.h
            fmm2d.h
)

#add_executable(${PROJECT_NAME} ${SRC_LIST})
add_executable(${PROJECT_NAME} ${SOURCES} ${HEADERS})
