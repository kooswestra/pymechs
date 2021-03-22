target_sources(cmechanism
    PUBLIC
        dna.h
        mechanism.h
        objective.h
        mutator.h
        pickandplace.h
    PRIVATE
        mechanism.cpp
        dna.cpp
        objective.cpp
        mutator.cpp
        pickandplace.cpp
    )

target_include_directories(cmechanism PUBLIC ${CMAKE_CURRENT_LIST_DIR})
