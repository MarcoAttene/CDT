cmake_minimum_required(VERSION 3.10)

# Use CMake -DLGPL=ON if you opt for the LGPL license
option(LGPL "LGPL" OFF)

# Only one of the following USE options can be ON

# Indirect Predicates (recommended)
option(USE_INDIRECT_PREDS "USE_INDIRECT_PREDS" ON) 

# These options are only useful for comparisons
option(USE_DOUBLE "USE_DOUBLE" OFF)
option(USE_LAZY_CORE "USE_LAZY_CORE" OFF)
option(USE_PLAIN_CORE "USE_PLAIN_CORE" OFF)
option(USE_LAZY_GMPQ "USE_LAZY_GMPQ" OFF)
option(USE_PLAIN_GMPQ "USE_PLAIN_GMPQ" OFF)

# set the project name
project(cdt)

# add path to cmake module path
list(PREPEND CMAKE_MODULE_PATH
    ${CMAKE_CURRENT_LIST_DIR}/cmake/
)

# specify the C++ standard
set(CMAKE_CXX_STANDARD 17)
set(CMAKE_CXX_STANDARD_REQUIRED True)

# add the executable
add_executable(${PROJECT_NAME} 
    src/giftWrap.cpp
    src/PLC.cpp
    src/delaunay.cpp
    src/main.cpp
)

if(CMAKE_CXX_COMPILER_ID STREQUAL "MSVC")
# optimize for speed
	target_compile_options(${PROJECT_NAME} PUBLIC /Ot)
# grant IEEE 754 compliance
	target_compile_options(${PROJECT_NAME} PUBLIC "/fp:strict")
# use intrinsic functions
	target_compile_options(${PROJECT_NAME} PUBLIC "/Oi")
# turn off annoying warnings
	target_compile_options(${PROJECT_NAME} PUBLIC "/D _CRT_SECURE_NO_WARNINGS")
# reserve enough stack size
	target_link_options(${PROJECT_NAME} PUBLIC "/STACK:8421376")
	if(NOT LGPL)
		target_compile_options(${PROJECT_NAME} PUBLIC "/D USE_MAROTS_METHOD")
	endif()
else()
# set standard optimization level
	target_compile_options(${PROJECT_NAME} PUBLIC -O2)
# grant IEEE 754 compliance
	target_compile_options(${PROJECT_NAME} PUBLIC -frounding-math)
# reserve enough stack size	
	target_compile_options(${PROJECT_NAME} PUBLIC -Wl,-z,stacksize=8421376)
	if(NOT LGPL)
		target_compile_options(${PROJECT_NAME} PUBLIC "-DUSE_MAROTS_METHOD")
	endif()
endif()

# Public include directory
target_include_directories(${PROJECT_NAME} PUBLIC
	src
	include
)

# Add define for indirect predicates
if(USE_INDIRECT_PREDS)
	target_compile_definitions(${PROJECT_NAME} PRIVATE USE_INDIRECT_PREDS)
endif()

# Add define for double
if(USE_DOUBLE)
target_compile_definitions(${PROJECT_NAME} PRIVATE USE_DOUBLE)
endif()

# The other number types require to link against CGAL/GMP/MPFR
if(USE_LAZY_CORE OR USE_PLAIN_CORE OR USE_LAZY_GMPQ OR USE_PLAIN_GMPQ)
	include(cgal)
	target_link_libraries(${PROJECT_NAME} CGAL::CGAL CGAL::CGAL_Core)

	if(USE_LAZY_CORE)
		target_compile_definitions(${PROJECT_NAME} PRIVATE USE_LAZY_CORE)
	endif()

	if(USE_PLAIN_CORE)
		target_compile_definitions(${PROJECT_NAME} PRIVATE USE_PLAIN_CORE)
	endif()

	if(USE_LAZY_GMPQ)
		target_compile_definitions(${PROJECT_NAME} PRIVATE USE_LAZY_GMPQ)
	endif()

	if(USE_PLAIN_GMPQ)
		target_compile_definitions(${PROJECT_NAME} PRIVATE USE_PLAIN_GMPQ)
	endif()
endif()
