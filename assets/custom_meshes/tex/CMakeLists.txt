make_directory(tex)
file(GLOB texture_files *.png *.jpg *.bmp)
file(COPY ${texture_files} DESTINATION ${CMAKE_BINARY_DIR}/custom_meshes/tex)