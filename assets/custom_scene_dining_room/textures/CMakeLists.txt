make_directory(textures)
file(GLOB texture_files *.png *.jpg *.hdr)
file(COPY ${texture_files} DESTINATION ${CMAKE_BINARY_DIR}/custom_scene_dining_room/texures)