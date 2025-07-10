add_executable(${PROJECT_NAME}
    main.cpp
    kicad.cpp
    kicad.hpp
)
target_link_libraries(${PROJECT_NAME}
   #OpenXLSX::OpenXLSX
   xlnt::xlnt
)

# install
install(TARGETS ${PROJECT_NAME})
