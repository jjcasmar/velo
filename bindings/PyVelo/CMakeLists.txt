find_package(pybind11)

set(SOURCES
    PyVelo.cpp
)

set(HEADERS
)

pybind11_add_module(PyVelo ${SOURCES} ${HEADERS})
target_link_libraries(PyVelo PUBLIC VeloPhysics pybind11::module pybind11::headers)