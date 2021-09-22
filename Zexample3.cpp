#include <iostream>
#include <clmagic/calculation/number/floating.h>

int main() {
  {
    using Number = calculation::FloatX<7, 3>;
    Number a = 1;
    std::cout << "a = " << a << std::endl;
    std::cout << "a + 7.5 = " << a + 7.5 << std::endl;
    std::cout << "a - 7.5 = " << a - 7.5 << std::endl;
    std::cout << "a * 7.5= " << a * 7.5 << std::endl;
    std::cout << "a / 7.5= " << a / 7.5 << std::endl;
  }

  {
    using Number = calculation::FloatX<2, 2>;
    /* Max = 0 10 11 = 1.11*exp2(10-01) = dex 1+0.5+0.25 * 2 = 3.5 */
    Number a = 3.5;
    std::cout << "a = " << a << std::endl;
    std::cout << "a + 0.5 = " << a + 0.5 << std::endl;
    std::cout << "a - 0.5 = " << a - 0.5 << std::endl;
    std::cout << "a * 0.5= " << a * 0.5 << std::endl;
    std::cout << "a / 0.5= " << a / 0.5 << std::endl;/* 0.5 cvt to 0.0, X/0.0 = inf */
    std::cout << "a + 1 = " << a + 1 << std::endl;
    std::cout << "a - 1 = " << a - 1 << std::endl;
    std::cout << "a * 1.5 = " << a * 1.5 << std::endl;
    std::cout << "a / 1.5= " << a / 1.5 << std::endl;

    /*std::cout << int(a) << std::endl;
    std::cout << int(a + 0.5) << std::endl;
    std::cout << int(a + 1.0) << std::endl;
    std::cout << int(a + 1.5) << std::endl;*/
    /*std::cout << (long long)(a) << std::endl;
    std::cout << (long long)(a + 0.5) << std::endl;
    std::cout << (long long)(a + 1.0) << std::endl;
    std::cout << (long long)(a + 1.5) << std::endl;*/
  }

  {
    using Number = calculation::FloatX<23, 8>;
    std::cout << Number(4294967000u) <<std::endl;
    std::cout << (unsigned int)Number(4294967000u) <<std::endl;
  }

  return 0;
}