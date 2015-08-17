project(PotatoLib)
cmake_minimum_required(VERSION 2.8)
#aux_source_directory(. SRC_LIST)

include_directories(/usr/local/include/)

set(SOURCES
            cell.cpp
            main.cpp
            tree.cpp
            point.cpp
            fmm.cpp
            fmm2d.cpp
            tree2d.cpp
            cell2d.cpp
            kernel_laplace_point_2d.cpp
            fmm_gmres_solver.cpp
)

set(HEADERS
            cell.h
            element.h
            point_element.h
            kernel.h
            tree.h
            point.h
            fmm.h
            fmm2d.h
            tree2d.h
            cell2d.h
            complex_t.h
            kernel_laplace_point_2d.h
            gmres.h
            fmm_gmres_solver.h
)

#add_executable(${PROJECT_NAME} ${SRC_LIST})
add_executable(${PROJECT_NAME} ${SOURCES} ${HEADERS})
