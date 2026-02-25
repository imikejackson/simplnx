# ------------------------------------------------------------------------------
# Required McPlotty and H5Support
# ------------------------------------------------------------------------------

#if(NOT SIMPLNX_BUILD_McPlotty)
#  find_package(H5Support REQUIRED)
#  find_package(McPlotty REQUIRED)
#else()
    if(NOT TARGET McPlotty::McPlotty)
        if(EXISTS "${simplnx_SOURCE_DIR}/../McPlotty")
            set(McPlotty_SOURCE_DIR "${simplnx_SOURCE_DIR}/../McPlotty")
        else()
            message(FATAL_ERROR "McPlotty_SOURCE_DIR was not set. Where is the McPlotty project directory. Please set the McPlotty_SOURCE_DIR variable to the McPlotty directory.")
        endif()
        message(STATUS "McPlotty_SOURCE_DIR: ${McPlotty_SOURCE_DIR}")


        set(McPlotty_ENABLE_Qt6 OFF)
        add_subdirectory( ${McPlotty_SOURCE_DIR} ${PROJECT_BINARY_DIR}/McPlotty)
    endif()
#endif()

