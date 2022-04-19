
add_compile_definitions(HAVEEIGEN=${HAVEEIGEN})

set(bin_PROGRAMS rotcurve densitymap fullintegrate comparediscs)

set(common_LINKLIB yaml-cpp ${EIEN3_LIBRARIES})

set(common_INCLUDE $<INSTALL_INTERFACE:include>
  $<BUILD_INTERFACE:${PROJECT_SOURCE_DIR}/include/>
  $<BUILD_INTERFACE:${PROJECT_SOURCE_DIR}/> ${CMAKE_BINARY_DIR}
  ${DEP_INC} ${CMAKE_CURRENT_SOURCE_DIR}
  ${CMAKE_CURRENT_SOURCE_DIR}/.. ${EIGEN3_INCLUDE_DIR})

add_executable(fullintegrate  fullintegrate.cc)
add_executable(densitymap     density_map.cc  )
add_executable(rotcurve       rotcurve.cc     )
add_executable(comparediscs   compare_discs.cc)

if(HAVEPYBIND11)
  pybind11_add_module(mwlmc mwlmc.cc)
endif()

foreach(program ${bin_PROGRAMS})
  target_link_libraries(${program} ${common_LINKLIB})
  target_include_directories(${program} PUBLIC ${common_INCLUDE})
  install(TARGETS ${program} DESTINATION bin)
endforeach()
