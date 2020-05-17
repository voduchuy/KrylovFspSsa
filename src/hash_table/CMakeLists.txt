set (MOD_NAME BRENT_HASH)
set (${MOD_NAME}_COMPONENTS
        big_integer_module
        HashTable
        )
set (EXTENSION f90)

foreach(MOD ${${MOD_NAME}_COMPONENTS})
    list(APPEND SRCS ${MOD}.${EXTENSION})
endforeach()

#add_library(${MOD_NAME}_OBJS OBJECT ${SRCS})
add_fortran_object(${MOD_NAME}_OBJS ${SRCS})

list(APPEND KFSSA_OBJS "$<TARGET_OBJECTS:${MOD_NAME}_OBJS>")
set(KFSSA_OBJS ${KFSSA_OBJS} PARENT_SCOPE)