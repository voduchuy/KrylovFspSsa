set (MOD_NAME STATE_SPACE)
set (${MOD_NAME}_COMPONENTS
        StateSpace
        )
set (EXTENSION f90)
set (DEPS MODEL_OBJS BRENT_HASH_OBJS)

foreach(MOD ${${MOD_NAME}_COMPONENTS})
    list(APPEND SRCS ${MOD}.${EXTENSION})
endforeach()

message(${KFSSA_OBJS})

add_fortran_object(${MOD_NAME}_OBJS ${SRCS} )
target_link_libraries(${MOD_NAME}_OBJS ${DEPS})

list(APPEND KFSSA_OBJS "$<TARGET_OBJECTS:${MOD_NAME}_OBJS>")
set(KFSSA_OBJS ${KFSSA_OBJS} PARENT_SCOPE)