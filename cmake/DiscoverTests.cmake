include_guard(GLOBAL)

# Register every visible test once. Only OOC-enabled builds split discovery to label
# cases tagged [OOC] (Catch2 tags ignore case) for `ctest -L "^ooc$"`.
# PROPERTIES: optional ctest properties that are set on every discovered test.
function(simplnx_discover_tests UNIT_TEST_TARGET)
  cmake_parse_arguments(ARGS "" "" "PROPERTIES" ${ARGN})
  if(DREAM3D_ENABLE_OOC)
    catch_discover_tests(${UNIT_TEST_TARGET} TEST_SPEC "~[OOC]~[.]" PROPERTIES ${ARGS_PROPERTIES})
    catch_discover_tests(${UNIT_TEST_TARGET} TEST_SPEC "[OOC]~[.]" TEST_LIST ${UNIT_TEST_TARGET}_OOC_TESTS PROPERTIES LABELS "ooc" ${ARGS_PROPERTIES})
  else()
    catch_discover_tests(${UNIT_TEST_TARGET} TEST_SPEC "~[.]" PROPERTIES ${ARGS_PROPERTIES})
  endif()
endfunction()
