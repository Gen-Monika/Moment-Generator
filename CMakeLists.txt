cmake_minimum_required(VERSION 3.21)

project(MomentGenerator
    VERSION 0.1.0
    DESCRIPTION "Qt + KaTeX central moment formula generator"
    LANGUAGES CXX
)

set(CMAKE_CXX_STANDARD 17)
set(CMAKE_CXX_STANDARD_REQUIRED ON)
set(CMAKE_AUTOMOC ON)
set(CMAKE_AUTORCC ON)
set(CMAKE_AUTOUIC ON)

find_package(Qt6 6.5 REQUIRED COMPONENTS Widgets WebEngineWidgets)

qt_add_executable(MomentGenerator WIN32
    src/main.cpp
    src/MainWindow.cpp
    src/MainWindow.h
    src/MomentEngine.cpp
    src/MomentEngine.h
)

target_link_libraries(MomentGenerator
    PRIVATE
        Qt6::Widgets
        Qt6::WebEngineWidgets
)

target_compile_definitions(MomentGenerator
    PRIVATE
        QT_NO_KEYWORDS
)
