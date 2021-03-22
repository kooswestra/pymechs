target_sources(cmechanism
    PUBLIC
        controller.h
    PRIVATE
        pid.h
        pid.cpp
        controller.cpp
        feedforward.h
        feedforward.cpp
    )

target_include_directories(cmechanism PUBLIC ${CMAKE_CURRENT_LIST_DIR})
