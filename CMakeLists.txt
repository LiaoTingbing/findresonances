# CMakeList.txt: findresonances 的 CMake 项目，在此处包括源代码并定义
# 项目特定的逻辑。
#
cmake_minimum_required(VERSION 3.10)


set(CMAKE_BUILD_PARALLEL_LEVEL 4)
project ("findresonances")


include_directories(
	E:/VisualStudioFiles/vcpkg/ku/armadillo-14.4.1/include

)

# 将源代码添加到此项目的可执行文件。
add_executable(findresonances "src/main.cpp" "include/main.h" "include/findresonances.h" "src/findresonances.cpp"
)


# TODO: 如有需要，请添加测试并安装目标。

target_link_libraries(
	${PROJECT_NAME}
	E:/VisualStudioFiles/vcpkg/ku/armadillo-14.4.1/lib_win64/libopenblas.lib
)