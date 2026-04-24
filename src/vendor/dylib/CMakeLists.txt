target_sources(kynema_fmb_library PRIVATE)

install(FILES
        dylib.hpp
        LICENSE.txt
        DESTINATION include/KynemaFMB/vendor/dylib/
)
