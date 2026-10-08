add_library(inflator SHARED Inflator.cpp)
add_library(inflator::inflator ALIAS inflator)
target_compile_definitions(inflator PRIVATE INFLATOR_BUILDING_LIBRARY)
target_compile_features(inflator PUBLIC cxx_std_17)
target_include_directories(inflator PUBLIC
    $<BUILD_INTERFACE:${PROJECT_SOURCE_DIR}/include>
    $<INSTALL_INTERFACE:include>)
target_link_libraries(inflator PRIVATE micro::isosurface_inflator micro::mesh_fem)
set_target_properties(inflator PROPERTIES VERSION 1.0.0 SOVERSION 1)
if(UNIX AND NOT APPLE)
    target_link_options(inflator PRIVATE
        "LINKER:--version-script=${CMAKE_CURRENT_SOURCE_DIR}/exports.map"
        "LINKER:--exclude-libs,ALL")
endif()

include(GNUInstallDirs)
install(TARGETS inflator EXPORT inflatorTargets
    LIBRARY DESTINATION ${CMAKE_INSTALL_LIBDIR}
    RUNTIME DESTINATION ${CMAKE_INSTALL_BINDIR}
    ARCHIVE DESTINATION ${CMAKE_INSTALL_LIBDIR})
install(FILES ${PROJECT_SOURCE_DIR}/include/inflator/Inflator.hpp DESTINATION ${CMAKE_INSTALL_INCLUDEDIR}/inflator)
install(EXPORT inflatorTargets NAMESPACE inflator:: DESTINATION ${CMAKE_INSTALL_LIBDIR}/cmake/inflator)
file(WRITE ${CMAKE_CURRENT_BINARY_DIR}/inflatorConfig.cmake "include(\"\${CMAKE_CURRENT_LIST_DIR}/inflatorTargets.cmake\")\n")
install(FILES ${CMAKE_CURRENT_BINARY_DIR}/inflatorConfig.cmake DESTINATION ${CMAKE_INSTALL_LIBDIR}/cmake/inflator)

if(MICRO_BUILD_LIBRARY_TESTS)
    enable_testing()
    add_executable(inflator_api_test test_api.cpp)
    target_link_libraries(inflator_api_test PRIVATE inflator::inflator)
    add_test(NAME inflator_api_2d COMMAND inflator_api_test 2)
    add_test(NAME inflator_api_3d COMMAND inflator_api_test 3)
    set_tests_properties(inflator_api_2d inflator_api_3d PROPERTIES TIMEOUT 120)
endif()
