#include <iostream>
#include <vector>

// see https://codereview.stackexchange.com/a/136373
#pragma once

#if __has_include(<optional>)

#   include <optional>
namespace stdx {
  using namespace ::std;
}
#elif __has_include(<experimental/optional>)
#   include <experimental/optional>
namespace stdx {
  using namespace ::std;
  using namespace ::std::experimental;
}

#else
#   error <experimental/optional> and <optional> not found
#endif


//GiNaC includes
#include <ginac/matrix.h>
#include <ginac/print.h>
#include <ginac/basic.h>
#include <ginac/function.h>
#include <ginac/symbol.h>
#include <ginac/mul.h>
#include <ginac/add.h>
#include <ginac/power.h>
#include <ginac/operators.h>

namespace g = GiNaC; 
typedef std::vector< g::matrix > mat_vec;