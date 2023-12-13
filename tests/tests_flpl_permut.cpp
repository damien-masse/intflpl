#include <cstdio>
#include "catch.hpp"
#include <ibex.h>
#include <flpl.h>

using namespace Catch;
using namespace Detail;
using namespace std;
using namespace ibex;
using namespace intflpl;

TEST_CASE("Permut")
{
  SECTION("Test Permut basic")
  {
    Permut P(3,6);
    P >> 2 >> 4;
    CHECK(!P.is_full());
    P.addNewImg(5);
    CHECK(P.is_full());
    CHECK(P.img(0)==2);
    CHECK(P.rev(4)==1);
    CHECK(P.rev(3)==-1);
    P.exchangeIdx(1,2);
    CHECK(P.rev(4)==2);
    CHECK(P.img(1)==5);
    P.removeAndAddBack(1,3);
    CHECK(P.img(1)==4);
    CHECK(P.img(2)==3);
    CHECK(P.rev(4)==1);
    CHECK(P.rev(3)==2);
    P.reset();
    CHECK(!P.is_defined(1));
    CHECK(!P.is_image(3));
  }
}

