set(dir "${CMAKE_CURRENT_SOURCE_DIR}")

list(
    APPEND srcs
    "${dir}/gauss.cpp"
    "${dir}/main.cpp"
    "${dir}/Structure.cpp"
)

set(srcs "${srcs}" PARENT_SCOPE)