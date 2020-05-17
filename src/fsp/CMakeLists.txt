set (MOD_NAME KRYLOVFSP)
set (${MOD_NAME}_COMPONENTS
        KrylovSolver
        )
set (EXTENSION f90)
set (DEPS STATE_SPACE_OBJS EXPOKIT_OBJS PARSER_OBJS)

foreach(MOD ${${MOD_NAME}_COMPONENTS})
    list(APPEND SRCS ${MOD}.${EXTENSION})
endforeach()

add_fortran_object(${MOD_NAME}_OBJS ${SRCS})
target_link_libraries(${MOD_NAME}_OBJS ${DEPS})

list(APPEND KFSSA_OBJS "$<TARGET_OBJECTS:${MOD_NAME}_OBJS>")
set(KFSSA_OBJS ${KFSSA_OBJS} PARENT_SCOPE)