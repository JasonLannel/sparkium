add_subdirectory(tex)

make_directory(custom_meshes)
file(GLOB mesh_files *.obj *.ply *.mtl)
file(COPY ${mesh_files} DESTINATION ${CMAKE_BINARY_DIR}/custom_meshes)