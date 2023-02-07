
set(bin_PROGRAMS rotcurve)

set(common_LINKLIB ${EIGEN3_LIBRARIES})

set(common_INCLUDE $<INSTALL_INTERFACE:include>
  $<BUILD_INTERFACE:${PROJECT_SOURCE_DIR}/include/>
  $<BUILD_INTERFACE:${PROJECT_SOURCE_DIR}/src/>
  $<BUILD_INTERFACE:${PROJECT_SOURCE_DIR}/> ${CMAKE_BINARY_DIR}
  ${DEP_INC} ${CMAKE_CURRENT_SOURCE_DIR}
  ${CMAKE_CURRENT_SOURCE_DIR}/.. ${EIGEN3_INCLUDE_DIR})


add_executable(rotcurve       rotcurve.cc     )

foreach(program ${bin_PROGRAMS})
  target_link_libraries(${program} ${common_LINKLIB})
  target_include_directories(${program} PUBLIC ${common_INCLUDE})
  install(TARGETS ${program} DESTINATION bin)
endforeach()
