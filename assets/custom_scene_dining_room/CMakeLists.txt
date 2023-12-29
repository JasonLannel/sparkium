add_subdirectory(textures)
add_subdirectory(models)

make_directory(custom_scene_dining_room)
file(GLOB scene_files *.xml *.hdr)
file(COPY ${scene_files} DESTINATION ${CMAKE_BINARY_DIR}/custom_scene_dining_room)