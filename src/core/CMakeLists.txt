# ==================================================================
#  invariant-lib - cmake configuration file
# ==================================================================


# Sources

set(SRC
  permut.cpp
  intLU.cpp
  intsimplex.cpp
  intcstrvect.cpp
  intpoly.cpp
)

set(HEADERS
  flpl_def.h
  permut.h
  intLU.h
  intsimplex.h
  intcstrvect.h
  intpoly.h
  intpoly_inline.hpp
  intflpl.h
)


add_library(intflpl #SHARED
    ${SRC}
    ${HEADERS}
)

target_link_libraries(intflpl
    ${LIBS} Ibex::ibex
)

target_include_directories(intflpl
    PUBLIC ibex
)

install(FILES ${HEADERS} DESTINATION include/intflpl)
install(TARGETS intflpl DESTINATION lib/intflpl)
