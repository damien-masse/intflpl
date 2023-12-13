# ==================================================================
#  codac / tests - cmake configuration file
# ==================================================================

set(TESTS_NAME flpl-tests)

list(APPEND SRC_TESTS 

  main.cpp
  tests_flpl_permut.cpp
  tests_flpl_LU.cpp
  tests_flpl_simplex.cpp
  tests_flpl_cstrvect.cpp
  tests_flpl_poly.cpp

)

add_executable(${TESTS_NAME} ${SRC_TESTS})
set(FLPL_HEADERS_DIR ../src/core catch/)
target_include_directories(${TESTS_NAME} SYSTEM PUBLIC ibex ${FLPL_HEADERS_DIR})
target_link_libraries(${TESTS_NAME} PUBLIC Ibex::ibex intflpl)
add_dependencies(check ${TESTS_NAME})
add_test(NAME ${TESTS_NAME} COMMAND ${TESTS_NAME})
