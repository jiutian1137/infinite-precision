# float_infinite_precision
float_&lt;23,8>, float_&lt;52,11>, float_&lt;112,15>,  . . . the library target is  infinite precision

clmagic/calculation/fundamental/float.h:{
	Author:"LongJiangnan",
	Date:"2019-2021",
	Owner:"LongJiangnan",
	License:"PIO(Please identify the Owner)",
	Ref:["https://www.rfwireless-world.com/Tutorials/floating-point-tutorial.html",
		 "https://github.com/gcc-mirror/gcc/blob/master/gcc/real.c"
		 "https://www.exploringbinary.com/binary-division/"],
	Example:{
		#include "clmagic/calculation/fundamental/float.h"
		#include <iostream>

		int main(int, char**) {
			using namespace::calculation;

			std::cout << to_string(float256(1.5)) << std::endl;
			std::cout << to_string(float256(1.5) * float256(2.0)) << std::endl;
			std::cout << to_string(float256(5.5) * float256(129.0)) << std::endl << std::endl;

			std::cout << to_string(float128(1.0) / float128(3.0)) << std::endl;
			std::cout << to_string(float128(1.0) / float128(9.0)) << std::endl;
			std::cout.precision(16);
			std::cout << double(float128(1.0) / (float128(9) + (float128(9) / float128(10)))) << std::endl;
			std::cout << to_string(float128(1.0) / (float128(9) + (float128(9) / float128(10)))) << std::endl << std::endl;

			std::cout << to_string(float256(1.0) / float256(3.0)) << std::endl;
			std::cout << to_string(float256(1.0) / float256(9.0)) << std::endl;
			std::cout << to_string(float256(1.0) / (float256(9) + (float256(9) / float256(10)))) << std::endl << std::endl;

			std::cout << to_string(float512(1.0) / float512(3.0)) << std::endl;
			std::cout << to_string(float512(1.0) / float512(9.0)) << std::endl;
			std::cout << to_string(float512(1.0) / (float512(9) + (float512(9) / float512(10)))) << std::endl;

			return 0;
		}
	}
}
