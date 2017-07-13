################################################################################
### Copyrights (c) 2015
### Marwan Abdellah <abdellah.marwan@gmail.com>
################################################################################

include( FindZLIB )
if( ZLIB_FOUND )
    set( LIBRARY_PATHS
            /usr/lib
            /usr/local/lib
            /sw/lib
            /opt/local/lib
            $ENV{PROGRAM_FILES}/OpenEXR/lib/static )

    find_path( OPENEXR_INCLUDE_PATH ImfRgbaFile.h
        PATH_SUFFIXES OpenEXR
            /usr/include
            /usr/local/include
            /sw/include
            /opt/local/include )

    include_directories()
    include_directories( ${OPENEXR_INCLUDE_PATH} )
    include_directories( ${OPENEXR_INCLUDE_PATH}/OpenEXR )
    include_directories( /usr/include )
    include_directories( /usr/local/include )
    include_directories( /opt/local/include )


    find_library( OPENEXR_HALF_LIBRARY
        NAMES Half
        PATHS ${LIBRARY_PATHS} )
    link_libraries( ${OPENEXR_HALF_LIBRARY} )

    find_library( OPENEXR_IEX_LIBRARY
        NAMES Iex
        PATHS ${LIBRARY_PATHS} )
    link_libraries( ${OPENEXR_IEX_LIBRARY} )

    find_library( OPENEXR_IMATH_LIBRARY
        NAMES Imath
        PATHS ${LIBRARY_PATHS} )
    link_libraries( ${OPENEXR_IMATH_LIBRARY} )

    find_library( OPENEXR_ILMIMF_LIBRARY
        NAMES IlmImf
        PATHS ${LIBRARY_PATHS} )
    link_libraries( ${OPENEXR_ILMIMF_LIBRARY} )

    find_library( OPENEXR_ILMTHREAD_LIBRARY
        NAMES IlmThread
        PATHS ${LIBRARY_PATHS} )
    link_libraries( ${OPENEXR_ILMTHREAD_LIBRARY} )
ENDIF( ZLIB_FOUND )

if( OPENEXR_INCLUDE_PATH AND
    OPENEXR_IMATH_LIBRARY AND
    OPENEXR_ILMIMF_LIBRARY AND
    OPENEXR_IEX_LIBRARY AND
    OPENEXR_HALF_LIBRARY
)
    set( OPENEXR_FOUND TRUE )
    set( OPENEXR_INCLUDE_PATHS ${OPENEXR_INCLUDE_PATH} CACHE STRING
        "The include paths needed to use OpenEXR" )
    include_directories( ${OPENEXR_INCLUDE_PATHS} )

    set( OPENEXR_LIBRARIES
            ${OPENEXR_IMATH_LIBRARY}
            ${OPENEXR_ILMIMF_LIBRARY}
            ${OPENEXR_IEX_LIBRARY}
            ${OPENEXR_HALF_LIBRARY}
            ${OPENEXR_ILMTHREAD_LIBRARY}
            ${ZLIB_LIBRARY} CACHE STRING "The libraries needed to use OpenEXR" )
    link_libraries( ${OPENEXR_LIBRARIES} )

endif( OPENEXR_INCLUDE_PATH AND
    OPENEXR_IMATH_LIBRARY AND
    OPENEXR_ILMIMF_LIBRARY AND
    OPENEXR_IEX_LIBRARY AND
    OPENEXR_HALF_LIBRARY )

if( OPENEXR_FOUND )
    if( NOT OPENEXR_FIND_QUIETLY )
        message( STATUS "Found OpenEXR: ${OPENEXR_ILMIMF_LIBRARY}" )
    endif( NOT OPENEXR_FIND_QUIETLY )
else( OPENEXR_FOUND )
    if( OPENEXR_FIND_REQUIRED )
        message( FATAL_ERROR "Could not find OpenEXR library" )
    endif( OPENEXR_FIND_REQUIRED )
endif( OPENEXR_FOUND )

mark_as_advanced(
    OPENEXR_INCLUDE_PATHS
    OPENEXR_LIBRARIES
    OPENEXR_ILMIMF_LIBRARY
    OPENEXR_IMATH_LIBRARY
    OPENEXR_IEX_LIBRARY
    OPENEXR_HALF_LIBRARY
)
