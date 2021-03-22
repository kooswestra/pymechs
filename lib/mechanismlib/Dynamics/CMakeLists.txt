target_sources(cmechanism
    PRIVATE
        connection.h
        connection.cpp
        constraintconnection.h
        constraintconnection.cpp
        hinge.h
        hinge.cpp
        motorelement.h
        motorelement.cpp
        motor.h
        motor.cpp
        dynamicelement.h
        dynamicelement.cpp
        spring.h
        spring.cpp
        torsionspring.h
        torsionspring.cpp
        stateequation.h
        stateequation.cpp
        link.h
        link.cpp
        endeffector.h
        endeffector.cpp
    )

target_include_directories(cmechanism PUBLIC ${CMAKE_CURRENT_LIST_DIR})
