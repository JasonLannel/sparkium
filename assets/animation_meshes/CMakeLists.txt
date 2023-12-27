make_directory(animation_meshes)
file(GLOB mesh_files *.obj)
file(COPY ${mesh_files} DESTINATION ${CMAKE_BINARY_DIR}/animation_meshes)