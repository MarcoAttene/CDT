cmake_minimum_required(VERSION 3.10)
if(POLICY CMP0135)
  cmake_policy(SET CMP0135 NEW)
endif()

# Set the target architecture.
# All modern x86/x64 processors support AVX2.
# Older x86/x64 processors may support SSE2 but not AVX2.
# Very old x86/x64 processors, or non x86/x64 
# processors, do not support any of the two.
# Modern ARM processors support NEON.

# --- Uncomment max one of the following lines ---

#							#Generic build
#set(ENABLE_SSE2 True)		#Build for x86-64 with SSE2
set(ENABLE_AVX2 True)		#Build for x86-64 with AVX2
#set(ENABLE_ARMNEON True)	#Build for ARM with NEON

# --------- DON NOT EDIT BELOW THIS LINE! --------

set(CMAKE_CXX_STANDARD 20)
set(CMAKE_CXX_STANDARD_REQUIRED True)

# disable SSE and AVX on MacOS
if(APPLE)
set(ENABLE_SSE2 False)
set(ENABLE_AVX2 False)
endif()

# AVX includes SSE -> make that explicit
if(ENABLE_AVX2)
set(ENABLE_SSE2 True)
endif()

# Check that no incompatible architectures are slelected
if(ENABLE_ARMNEON AND ENABLE_SSE2)
message(FATAL_ERROR "ARMNEON and SSE2/AVX2 cannot be both set!")
endif()

# Use CMake -DLGPL=ON if you opt for the LGPL license
option(LGPL "LGPL" ON)

# set the project name
project(cdt)

# add the executable
add_executable(${PROJECT_NAME} 
    main.cpp
)


# Dependencies
include(FetchContent)

if(EXISTS ../.local.cmake)
	set(nfg_SOURCE_DIR "../NFG/")
	set(indirect_predicates_SOURCE_DIR "../Indirect_Predicates/")
	set(delaunay3d_SOURCE_DIR "../Delaunay3D/")
else()
	FetchContent_Declare(nfg
		GIT_REPOSITORY https://github.com/MarcoAttene/nfg.git
		GIT_TAG main
		SOURCE_SUBDIR do-not-configure
	)

	FetchContent_Declare(indirect_predicates
		GIT_REPOSITORY https://github.com/MarcoAttene/Indirect_Predicates.git
		SOURCE_SUBDIR do-not-configure
	)

	FetchContent_Declare(Delaunay3D
		GIT_REPOSITORY https://github.com/MarcoAttene/Delaunay3D.git
		GIT_TAG main
		SOURCE_SUBDIR do-not-configure
	)

	FetchContent_MakeAvailable(nfg indirect_predicates delaunay3d)
	message(STATUS "NFG is available in " ${nfg_SOURCE_DIR})
	message(STATUS "Indirect_Predicates is available in " ${indirect_predicates_SOURCE_DIR})
	message(STATUS "Delaunay3D is available in " ${delaunay3d_SOURCE_DIR})
endif()

if(ENABLE_ARMNEON)
	FetchContent_Declare(simde
		GIT_REPOSITORY https://github.com/simd-everywhere/simde.git
		GIT_TAG v0.8.2
	)
	FetchContent_MakeAvailable(simde)
	message(STATUS "simde is available in " ${simde_SOURCE_DIR})

	target_include_directories(${PROJECT_NAME} PUBLIC ${simde_SOURCE_DIR}/simde)
	target_compile_options(${PROJECT_NAME} PUBLIC -DSIMDE_ENABLE_NATIVE_ALIASES)
endif()

# Compiler-specific options
if(CMAKE_CXX_COMPILER_ID STREQUAL "MSVC")
	# grant IEEE 754 compliance; remove security controls; turn off annoying warnings
	target_compile_options(${PROJECT_NAME} PUBLIC "/fp:strict;/GS-;/D _CRT_SECURE_NO_WARNINGS")
	# set target architecture
	if(ENABLE_AVX2)
		target_compile_options(${PROJECT_NAME} PUBLIC "/arch:AVX2")
	elseif(ENABLE_SSE2)
		target_compile_options(${PROJECT_NAME} PUBLIC "/arch:SSE2")
	elseif(ENABLE_ARMNEON)
		target_compile_options(${PROJECT_NAME} PUBLIC "/arch:ARMV8.0")
	endif()
else()
	# set max optimization
	target_compile_options(${PROJECT_NAME} PUBLIC -O3)
	# grant IEEE 754 compliance
	target_compile_options(${PROJECT_NAME} PUBLIC -frounding-math)
	# the exact predicates rely on error-free transformations (e.g. Dekker
	# two-product): forbid the compiler from contracting a*b+c into an FMA, which
	# rounds once instead of twice and silently breaks those identities. This
	# mirrors the /fp:strict used on the MSVC path.
	target_compile_options(${PROJECT_NAME} PUBLIC -ffp-contract=off)
	# set target architecture
	if(ENABLE_AVX2)
		target_compile_options(${PROJECT_NAME} PUBLIC "-mavx2")
		target_compile_options(${PROJECT_NAME} PUBLIC "-mfma")
	elseif(ENABLE_SSE2)
		target_compile_options(${PROJECT_NAME} PUBLIC "-msse2")
	elseif(ENABLE_ARMNEON)
		target_compile_options(${PROJECT_NAME} PUBLIC "-march=armv8-a+simd")
	endif()
endif()

# Public include directory
target_include_directories(${PROJECT_NAME} PUBLIC
	include
	${nfg_SOURCE_DIR}/include
	${indirect_predicates_SOURCE_DIR}/include
	${delaunay3d_SOURCE_DIR}/include
)

# Tests. test_cdt tetrahedrizes one OFF file and audits the result -- no ghosts or dangling
# adjacency in the solid region, a closed and consistently oriented boundary, positive tets, and the
# volume the PLC encloses. Run them with `ctest --test-dir build`.
#
# The tests judge by the CDT3 TEST PASSED line rather than the exit status: a failure raised inside
# the library ends the process through ip_error(), which upstream implements as exit(0), so a
# status-only test would read a mid-run abort as a pass.
enable_testing()

# The test is built the same way as the main executable -- same include path, same floating-point
# and architecture flags. The exact predicates depend on those flags, and the test reads coordinates
# those predicates produced, so it must not be compiled under different rules.
add_executable(test_cdt test_cdt.cpp)
target_include_directories(test_cdt PUBLIC
	include
	${nfg_SOURCE_DIR}/include
	${indirect_predicates_SOURCE_DIR}/include
	${delaunay3d_SOURCE_DIR}/include
)
if(CMAKE_CXX_COMPILER_ID STREQUAL "MSVC")
	target_compile_options(test_cdt PUBLIC "/fp:strict;/GS-;/D _CRT_SECURE_NO_WARNINGS")
	if(ENABLE_AVX2)
		target_compile_options(test_cdt PUBLIC "/arch:AVX2")
	elseif(ENABLE_SSE2)
		target_compile_options(test_cdt PUBLIC "/arch:SSE2")
	elseif(ENABLE_ARMNEON)
		target_compile_options(test_cdt PUBLIC "/arch:ARMV8.0")
	endif()
else()
	target_compile_options(test_cdt PUBLIC -O3 -frounding-math -ffp-contract=off)
	if(ENABLE_AVX2)
		target_compile_options(test_cdt PUBLIC "-mavx2" "-mfma")
	elseif(ENABLE_SSE2)
		target_compile_options(test_cdt PUBLIC "-msse2")
	elseif(ENABLE_ARMNEON)
		target_compile_options(test_cdt PUBLIC "-march=armv8-a+simd")
	endif()
endif()
if(ENABLE_ARMNEON)
	target_include_directories(test_cdt PUBLIC ${simde_SOURCE_DIR}/simde)
	target_compile_options(test_cdt PUBLIC -DSIMDE_ENABLE_NATIVE_ALIASES)
endif()

function(cdt_add_off_test name file)
	add_test(NAME ${name} COMMAND test_cdt ${CMAKE_CURRENT_SOURCE_DIR}/Input_file/${file} ${ARGN})
	set_tests_properties(${name} PROPERTIES PASS_REGULAR_EXPRESSION "CDT3 TEST PASSED")
endfunction()

# A single cell of a background grid, as an immersed-boundary mesher hands it to CDT3: a closed,
# manifold, consistently oriented 69-vertex shell that CDT3 used to crash on.
#
# Its PLC face 28 has five vertices lying exactly in the plane of the face itself.
# getTetsIntersectingFace() gathered the tet stars of those flat vertices and de-duplicated them
# with std::unique WITHOUT sorting first, so a tet shared by two stars survived twice.
# recoverFaceHSi() then moved the cavity tets to the mesh tail by index, rewriting only the entry it
# held a reference to, and each duplicate dragged an unrelated tet into the cavity. A face of one
# such intruder lay in the plane of face 28, which makes it a hull facet of the cavity's Delaunay --
# one side a ghost that markInnerTets() pins to DT_OUT -- so neither side could be classified DT_IN,
# and meshCavity() reconnected the mesh to a tet it was about to delete. The segfault landed two PLC
# faces later, in an unrelated walk, which is why the audit checks the whole mesh.
#
# The expected volume is what the input shell encloses (the cell is 1x1x1 and the solid fills 0.646
# of it), computed by the divergence theorem over the input triangles; a CDT of a closed PLC must
# reproduce it. The tolerance is 1e-7 rather than machine epsilon because the sum runs over the
# DOUBLE approximations of the Steiner points, which are implicit points held exactly inside the
# library -- 1e-7 is still seven significant digits here, while any topological defect (a lost or
# doubled tet, an unrecovered face) moves the volume by a thousand times more.
cdt_add_off_test(cdt_flat_cavity_face flat_cavity_face.off 0.64642867950838023 1e-7)

# The sample model shipped with CDT3, as a check that the audit above is not tuned to one input.
cdt_add_off_test(cdt_boeing_part boeing_part.off)
