SET(SRC_TRIMESHTRACER
	main.cpp
	maths.h
	maths.cpp
	scene.h
	scene.cpp
	
	./external/objparser.h
	./external/objparser.cpp
	
	./external/sokol_time.h
	./external/stb_image_write.h
)

add_executable(trimeshtracer ${SRC_TRIMESHTRACER})
target_include_directories(trimeshtracer
	PUBLIC
		${CMAKE_CURRENT_SOURCE_DIR}
		${CMAKE_CURRENT_SOURCE_DIR}/external
)

add_executable(trimeshtracer_omp ${SRC_TRIMESHTRACER})
target_include_directories(trimeshtracer_omp
	PUBLIC
		${CMAKE_CURRENT_SOURCE_DIR}
		${CMAKE_CURRENT_SOURCE_DIR}/external
)

IF(WIN32)
	set_target_properties(trimeshtracer_omp PROPERTIES COMPILE_FLAGS "/openmp")
ENDIF(WIN32)
IF(UNIX AND NOT APPLE)
	set_target_properties(trimeshtracer_omp PROPERTIES COMPILE_FLAGS "-fopenmp")
	set_target_properties(trimeshtracer_omp PROPERTIES LINK_FLAGS "-fopenmp") 
ENDIF(UNIX AND NOT APPLE)
