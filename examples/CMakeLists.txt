set(PROGS
        toggle
        transcr6d
        repressilator
        )
set(DEPS PARSER_OBJS MODEL_OBJS STATE_SPACE_OBJS KRYLOVFSP_OBJS EXPOKIT_OBJS -lblas -llapack)

foreach ( PROG ${PROGS})
    add_executable(${PROG} ${PROG}.f90)
    target_link_libraries(${PROG} KFSSA ${DEPS})
endforeach()