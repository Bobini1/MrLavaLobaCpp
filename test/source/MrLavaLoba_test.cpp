#include <catch2/catch_test_macros.hpp>

#include "lib.hpp"

TEST_CASE("Name is MrLavaLoba", "[library]")
{
  auto const lib = library {};
  REQUIRE(lib.name == "MrLavaLoba");
}
