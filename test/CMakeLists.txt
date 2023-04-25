# Parent project does not export its library target, so this CML implicitly
# depends on being added from it, i.e. the testing is done only from the build
# tree and is not feasible from an install location

project(MrLavaLobaTests LANGUAGES CXX)

# ---- Dependencies ----

find_package(Catch2 REQUIRED)
include(Catch)

# ---- Tests ----

add_executable(MrLavaLoba_test source/MrLavaLoba_test.cpp)
target_link_libraries(
    MrLavaLoba_test PRIVATE
    MrLavaLoba_lib
    Catch2::Catch2WithMain
)
target_compile_features(MrLavaLoba_test PRIVATE cxx_std_20)

catch_discover_tests(MrLavaLoba_test)

# ---- End-of-file commands ----

add_folders(Test)
