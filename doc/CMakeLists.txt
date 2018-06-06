# doxygen target
find_package(Doxygen)

if (DOXYGEN_FOUND)
  # subdirs to be processed by doxygen
  set(DOXYGEN_SRC_DIRS 
    ${PROJECT_SOURCE_DIR}/src
    ${PROJECT_SOURCE_DIR}/doc/concepts
    ${PROJECT_BINARY_DIR}/doc/extra_dox)

  # create space-separated list of documentation input directories
  string(REGEX REPLACE ";" " " DOXYGEN_INPUT "${DOXYGEN_SRC_DIRS}")

  # set path to images
  set(DOXYGEN_BASE ${PROJECT_SOURCE_DIR}/doc/doxygen)
  set(DOXYGEN_IMAGE_PATH ${PROJECT_BINARY_DIR}/doc/doxygen/images)

  configure_file(${PROJECT_SOURCE_DIR}/doc/Doxyfile.in ${PROJECT_BINARY_DIR}/doc/Doxyfile)
  configure_file(${PROJECT_SOURCE_DIR}/doc/concepts/main_page.dox.in ${PROJECT_BINARY_DIR}/doc/extra_dox/main_page.dox)

  # Doxygen target:
  # 1) copy over doc/doxygen directory from source tree
  # 2) remove doc/doxygen/html if it exists
  # 3) run doxygen Doxyfile to regenerate the documentation
  add_custom_target(doxygen 
    COMMAND ${CMAKE_COMMAND} -E copy_directory ${PROJECT_SOURCE_DIR}/doc/doxygen ${PROJECT_BINARY_DIR}/doc/doxygen
    COMMAND ${CMAKE_COMMAND} -E remove_directory ${PROJECT_BINARY_DIR}/doc/doxygen/html
    COMMAND ${DOXYGEN_EXECUTABLE} ${PROJECT_BINARY_DIR}/doc/Doxyfile)

  # Create empty build directory, to avoid error if make doxygen is not run.
  execute_process(COMMAND ${CMAKE_COMMAND} -E make_directory ${PROJECT_BINARY_DIR}/doc/doxygen)

  # Install doxygen documentation if it has been generated.
  install(DIRECTORY ${PROJECT_BINARY_DIR}/doc/doxygen
    DESTINATION share/hiflow3/doc)

endif (DOXYGEN_FOUND)


