cmake_minimum_required(VERSION 3.10)

# Use CMake -DLGPL=ON if you opt for the LGPL license
option(LGPL "LGPL" OFF)

# set the project name
project(cdt)

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

# The following enable compilation and linking against CGAL on my own MSVC
# CGAL is used only for comparison purposes. It is safe to comment
# all the lines below this one for use with indirect predicates.

#if(CMAKE_CXX_COMPILER_ID STREQUAL "MSVC")
#	set(GMP_HOME_DIR "D:/SYNC_DATA/Sviluppo_Software/Third_Party/mpir-master/msvc/vs19")
#	set(GMPC_HOME_DIR ${GMP_HOME_DIR}/lib_mpir_gc/x64/Release)
#	set(GMPXX_HOME_DIR ${GMP_HOME_DIR}/lib_mpir_cxx/x64/Release)
#	set(GMP_INCLUDE_DIR ${GMPXX_HOME_DIR})
#	set(CGAL_DIR "D:/SYNC_DATA/vcpkg/packages/cgal_x64-windows/share/cgal")
#	set(GMP_LIBRARIES "D:/SYNC_DATA/vcpkg/installed/x64-windows/lib/gmp.lib")
#	set(MPFR_LIBRARIES "D:/SYNC_DATA/vcpkg/installed/x64-windows/lib/mpfr.lib")
#	set(MPFR_INCLUDE_DIR "D:/SYNC_DATA/vcpkg/installed/x64-windows/include")
#	set(Boost_INCLUDE_DIR "D:/SYNC_DATA/vcpkg/installed/x64-windows/include")
#	find_package(CGAL REQUIRED COMPONENTS Core)
#	target_link_libraries(${PROJECT_NAME} CGAL::CGAL CGAL::CGAL_Core)
#endif()
