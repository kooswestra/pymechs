target_sources(cmechanism
    PRIVATE
        solver.h
        solver.cpp
        rungekutta4.h
        rungekutta4.cpp
        trapezoid.h
        trapezoid.cpp
    )

target_include_directories(cmechanism PUBLIC ${CMAKE_CURRENT_LIST_DIR})
