foreach(dir libefp
            pybind11
            qcelemental
)
   add_subdirectory(${dir})
endforeach()

