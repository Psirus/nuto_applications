add_custom_application(PlateCooling)
add_custom_application(DynamicXd)

add_test(custom::MesoShell ${CMAKE_CURRENT_SOURCE_DIR}/MesoShell.py)
set_property(TEST custom::MesoShell PROPERTY ENVIRONMENT "PYTHONPATH=${CMAKE_BINARY_DIR}/python")
