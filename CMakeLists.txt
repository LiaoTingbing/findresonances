# CMakeList.txt: findresonances 的 CMake 项目，在此处包括源代码并定义
# 项目特定的逻辑。
#
cmake_minimum_required(VERSION 3.10)

 
project ("findresonances")

# 将源代码添加到此项目的可执行文件。
add_executable(findresonances "src/example1.cpp"  "include/findresonances.h" "src/findresonances.cpp"
 )


# TODO: 如有需要，请添加测试并安装目标。

find_package(Armadillo CONFIG REQUIRED)
target_link_libraries(findresonances PRIVATE armadillo)
