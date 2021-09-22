#include <clmagic/calculation/number/floating.h>
#include <iostream>

int main() {
  //using Number = calculation::FloatX<52,11>;
  using Number = calculation::Float128;//modern CPUs use xmm0 calculate double and float. we can view by expanded bits

  std::cout.precision(16);
  for (size_t i = 0; i < 10000; i++) {
    double a = rand() + rand() / double(RAND_MAX);
    double b = (rand() + rand() / double(RAND_MAX));
    double c = a + b;
    Number c1 = Number(a) + Number(b);
    std::cout << c << "\t" << c1 << std::endl;

    // not( relative error <= epsilon()/2 )
    if ( !(abs(c1 - c) <= (std::numeric_limits<double>::epsilon()/2)*c) ) {
      std::cout << "Error" << std::endl;
      std::cout << "error = " << abs(c1 - c) << std::endl;
      std::cout << "eps   = " << c * Number(std::numeric_limits<double>::epsilon()/2) << std::endl;
      std::cin.get();
    }
  }
  
  return 0;
}
