set (MOD_NAME EXPOKIT)
set (${MOD_NAME}_COMPONENTS
        dgchbv
        dgpadm
        clock
        )
foreach(MOD ${${MOD_NAME}_COMPONENTS})
    list(APPEND SRCS ${MOD}.f)
endforeach()

add_library(${MOD_NAME}_OBJS OBJECT ${SRCS})

list(APPEND KFSSA_OBJS "$<TARGET_OBJECTS:${MOD_NAME}_OBJS>")
set(KFSSA_OBJS ${KFSSA_OBJS} PARENT_SCOPE)

