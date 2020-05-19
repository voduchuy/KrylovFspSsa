function(add_fortran_object LIB)
    add_library(${LIB} OBJECT ${ARGN})

    # set module path to LIB_DIR/mod
    get_target_property(LIB_DIR ${LIB} BINARY_DIR)
    set_target_properties(${LIB} PROPERTIES Fortran_MODULE_DIRECTORY ${LIB_DIR}/mod)

    # making LIB_DIR/mod available for libraries linking LIB
    target_include_directories(${LIB} INTERFACE ${LIB_DIR}/mod)
endfunction(add_fortran_object)

add_subdirectory(expokit)
add_subdirectory(hash_table)
add_subdirectory(parser)
add_subdirectory(model)
add_subdirectory(state_space)
add_subdirectory(fsp)

add_library(KFSSA ${KFSSA_OBJS})
