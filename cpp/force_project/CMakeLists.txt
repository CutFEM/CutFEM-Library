set (FORCE_PROJECT_SOURCE
        ./input_handler.cpp
        )


add_library(force_project STATIC "")
add_library(force_project::force_project ALIAS force_project)
target_sources(force_project PRIVATE ${FORCE_PROJECT_SOURCE})
target_link_libraries(force_project PUBLIC ryml::ryml )

