#include <iostream>
#include <clmagic/calculation/number/floating.h>
#include <clmagic/calculation/number/function.h>
#include <clmagic/calculation/formulas/legendre.h>

int main() {
  using namespace::calculation;

  std::cout << "P(x) = 0, x = {" << std::endl;
  auto polynomial = LegendrePolynomial(10);
  Float128 roots[10];
  solve_zeros(polynomial, roots);
  for (size_t i = 0; i != 10; ++i) {
    std::cout << roots[i] << "," << std::endl;
  }
  std::cout << "}" << std::endl;

  /*
  P(x) = 0, x = {
  0.97390652851717172007796401208445210,
  0.86506336668898451073209668842349307,
  0.67940956829902440623432736511487355,
  0.43339539412924719079926594316578418,
  0.14887433898163121088482600112971998,
  -0.14887433898163121088482600112971998,
  -0.43339539412924719079926594316578418,
  -0.67940956829902440623432736511487355,
  -0.86506336668898451073209668842349307,
  -0.97390652851717172007796401208445210,
  }
  */

  return 0;
}