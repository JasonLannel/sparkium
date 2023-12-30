make_directory(models)
file(GLOB obj_files *.obj)
file(COPY ${obj_files} DESTINATION ${CMAKE_BINARY_DIR}/custom_scene_veach/models)