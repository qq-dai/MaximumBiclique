include_directories(../tools)

add_library(MBCS MBCS.h MBCS.cpp)
target_link_libraries(MBCS biGraph)