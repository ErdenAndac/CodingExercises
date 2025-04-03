add_library(eigen3 INTERFACE)
target_include_directories(eigen3
    INTERFACE
        $<BUILD_INTERFACE:${CMAKE_CURRENT_SOURCE_DIR}/eigen>
        $<INSTALL_INTERFACE:include>
)

# Create an alias target with the namespace
add_library(Eigen3::Eigen ALIAS eigen3)

# Add matplotlibcpp (Header-only)
add_library(matplotlibcpp INTERFACE)
target_include_directories(matplotlibcpp INTERFACE ${CMAKE_CURRENT_SOURCE_DIR}/matplotlib)

# Find Python3 for Matplotlib (Explicitly Set Paths for Python 3.10 and NumPy)
set(Python3_EXECUTABLE "C:/Users/erden/AppData/Local/Programs/Python/Python310/python.exe")
set(Python3_INCLUDE_DIR "C:/Users/erden/AppData/Local/Programs/Python/Python310/include")
set(Python3_LIBRARIES "C:/Users/erden/AppData/Local/Programs/Python/Python310/libs/python310.lib")

set(Python3_NumPy_INCLUDE_DIRS "C:/Users/erden/AppData/Local/Programs/Python/Python310/lib/site-packages/numpy/_core/include")

find_package(Python3 COMPONENTS Interpreter Development NumPy REQUIRED)

# Link Python3::Python to matplotlibcpp
target_include_directories(matplotlibcpp INTERFACE ${Python3_INCLUDE_DIR} ${Python3_NumPy_INCLUDE_DIRS})
target_link_libraries(matplotlibcpp INTERFACE Python3::Python)

# Create an alias target for better usage
add_library(MatplotlibCpp::Matplotlib ALIAS matplotlibcpp)

