# QuartetScores - Code for computing various support scores for internodes.
# Copyright (C) 2016-2017 Sarah Lutteropp and Lucas Czech
#
# This program is free software: you can redistribute it and/or modify
# it under the terms of the GNU General Public License as published by
# the Free Software Foundation, either version 3 of the License, or
# (at your option) any later version.
#
# This program is distributed in the hope that it will be useful,
# but WITHOUT ANY WARRANTY; without even the implied warranty of
# MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
# GNU General Public License for more details.
#
# You should have received a copy of the GNU General Public License
# along with this program.  If not, see <http://www.gnu.org/licenses/>.
#
# Contact:
# Sarah Lutteropp <sarah.lutteropp@h-its.org> or
# Lucas Czech <lucas.czech@h-its.org>
# Exelixis Lab, Heidelberg Institute for Theoretical Studies
# Schloss-Wolfsbrunnenweg 35, D-69118 Heidelberg, Germany

# --------------------------------------------------------------------------------------------------
#   Download Genesis
# --------------------------------------------------------------------------------------------------

IF( NOT EXISTS ${PROJECT_SOURCE_DIR}/genesis/CMakeLists.txt )
    message (STATUS "Genesis not found")
    message (STATUS "Will now download Genesis")

    # If Genesis was not found, we download and unpack it (at configure time). This roughly follows
    # https://github.com/google/googletest/tree/master/googletest#incorporating-into-an-existing-cmake-project

    configure_file(
        ${PROJECT_SOURCE_DIR}/cmake/GenesisDownload.cmake
        ${CMAKE_BINARY_DIR}/genesis-download/CMakeLists.txt
    )

    execute_process( COMMAND ${CMAKE_COMMAND} -G "${CMAKE_GENERATOR}" .
        RESULT_VARIABLE result
        WORKING_DIRECTORY ${CMAKE_BINARY_DIR}/genesis-download
    )

    if(result)
        message (STATUS "Cannot configure Genesis: ${result}")
        return()
    endif()

    execute_process( COMMAND ${CMAKE_COMMAND} --build .
        RESULT_VARIABLE result
        WORKING_DIRECTORY ${CMAKE_BINARY_DIR}/genesis-download
    )

    if(result)
        message (STATUS "Cannot build Genesis: ${result}")
        return()
    endif()

    # If the header still does not exists, something went wrong.
    IF( NOT EXISTS ${PROJECT_SOURCE_DIR}/genesis/CMakeLists.txt )
        message (STATUS "Genesis not found.")
        return()
    ENDIF()

    message (STATUS "Finished downloading Genesis")
ENDIF()

# --------------------------------------------------------------------------------------------------
#   Download Tclap
# --------------------------------------------------------------------------------------------------

IF( NOT EXISTS ${PROJECT_SOURCE_DIR}/tclap/include )
    message (STATUS "Tclap not found")
    message (STATUS "Will now download Tclap")

    # If Tclap was not found, we download and unpack it (at configure time). This roughly follows
    # https://github.com/google/googletest/tree/master/googletest#incorporating-into-an-existing-cmake-project

    configure_file(
        ${PROJECT_SOURCE_DIR}/cmake/TclapDownload.cmake
        ${CMAKE_BINARY_DIR}/tclap-download/CMakeLists.txt
    )

    execute_process( COMMAND ${CMAKE_COMMAND} -G "${CMAKE_GENERATOR}" .
        RESULT_VARIABLE result
        WORKING_DIRECTORY ${CMAKE_BINARY_DIR}/tclap-download
    )

    if(result)
        message (STATUS "Cannot configure Tclap: ${result}")
        return()
    endif()

    execute_process( COMMAND ${CMAKE_COMMAND} --build .
        RESULT_VARIABLE result
        WORKING_DIRECTORY ${CMAKE_BINARY_DIR}/tclap-download
    )

    if(result)
        message (STATUS "Cannot build Tclap: ${result}")
        return()
    endif()

    # If the header still does not exists, something went wrong.
    IF( NOT EXISTS ${PROJECT_SOURCE_DIR}/tclap/include )
        message (STATUS "Tclap not found.")
        return()
    ENDIF()

    message (STATUS "Finished downloading Tclap")
ENDIF()
