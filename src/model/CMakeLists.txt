set (MOD_NAME MODEL)
set (${MOD_NAME}_COMPONENTS
        ModelModule
        )
set (EXTENSION f90)
set (DEPS PARSER_OBJS)

foreach(MOD ${${MOD_NAME}_COMPONENTS})
    list(APPEND SRCS ${MOD}.${EXTENSION})
endforeach()

add_fortran_object(${MOD_NAME}_OBJS ${SRCS} ${KFSSA_OBJS})
target_link_libraries(${MOD_NAME}_OBJS ${DEPS})

list(APPEND KFSSA_OBJS "$<TARGET_OBJECTS:${MOD_NAME}_OBJS>")
set(KFSSA_OBJS ${KFSSA_OBJS} PARENT_SCOPE)