/* clmagic/calculation/fundamental/float.h:{
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
  },
  Example2:{
    #include "clmagic/calculation/fundamental/float.h"
    #include <iostream>
    using namespace::calculation;

    template<typename __traits, size_t __size>
    struct Spectrum {
      using Real = __traits;
      static constexpr double wavelength_lower = static_cast<double>(380);
      static constexpr double wavelength_upper = static_cast<double>(740);
      static constexpr size_t sampling_segments = __size;
      static constexpr double sampling_steplength = (wavelength_upper - wavelength_lower) / static_cast<double>(sampling_segments);
      
      Real operator()(Real lambda) const {
        assert(wavelength_lower <= lambda && lambda < wavelength_upper);
        return spectral_power_array[
          static_cast<size_t>(floor((lambda - wavelength_lower) / sampling_steplength))
        ];
      }

      Real spectral_power_array[sampling_segments];
    };

    int main(int, char**) {
      Spectrum<float128, 10> rod = {
        float128(1.0), float128(9)/10, float128(6)/10, float128(4)/10, float128(9)/10,
        float128(9)/10, float128(4)/10, float128(4)/10, float128(2)/10, float128(1)/10
      };
      
      for (float i = 385.f; i < 740.0f; i += 10.0f) {
        std::cout << i << " = " << to_string(rod(float128(i))) << std::endl;
      }
	
      return 0;
    }
  }
} */
#pragma once

#include <cassert>
#include <deque>
#include <algorithm>
#include <xstring>
namespace calculation {
	// { used for float_<X,X> to decimal to decimal_string }
	struct decimal {
		std::deque<uint8_t> numbers;
		// index
		size_t decimal_point = 0;
		// flag
		bool negative = false;
		// flags
		int mode = 0;

		decimal() = default;

		int compare_without_sign(const decimal& right) const {
			size_t this_integer_length = this->numbers.size() - this->decimal_point;
			size_t right_integer_length = right.numbers.size() - right.decimal_point;
			if (this_integer_length != right_integer_length) {
				return this_integer_length > right_integer_length ? 1 : -1;
			}

			// setup first1,last1,first2
			auto first1 = numbers.rbegin();
			auto last1  = numbers.rend();
			auto first2 = right.numbers.rbegin();
			auto last2  = right.numbers.rend();
			size_t this_decimal_length = this->numbers.size() - this_integer_length;
			size_t right_decimal_length = right.numbers.size() - right_integer_length;
			if (this_decimal_length < right_decimal_length) {
				last2 -= (right_decimal_length - this_decimal_length);
			} else {
				last1 -= (this_decimal_length - right_decimal_length);
			}

			// compare
			for ( ; first1 != last1 && first2 != last2; ++first1, ++first2) {
				if (*first1 != *first2) {
					return *first1 > *first2 ? 1 : -1;
				}
			}
			if (first1 != last1 && first2 == last2) {
				return -1;
			}
			if (first1 == last1 && first2 != last2) {
				return 1;
			}
			return 0;
		}

		void addev_without_sign(const decimal& right) {
			// push Zero-precision to lowest-bit(decimal)
			if ( this->decimal_point < right.decimal_point ) {
				size_t diff = right.decimal_point - this->decimal_point;
				for (size_t i = 0; i != diff; ++i) {
					this->numbers.push_front(0);
				}
				this->decimal_point += diff;
			}

			// push Zero-precision to highest-bit(integer)
			size_t this_integer_length  = this->numbers.size() - this->decimal_point;
			size_t right_integer_length = right.numbers.size() - right.decimal_point;
			if (this_integer_length < right_integer_length) {
				size_t diff = right_integer_length - this_integer_length;
				for (size_t i = 0; i != diff; ++i) {
					this->numbers.push_back(0);
				}
			}
			this->numbers.push_back(0);

			// addss this[i], right[i]
			auto this_number      = std::next(this->numbers.begin(), this->decimal_point - right.decimal_point);
			auto this_back_number = std::prev(this->numbers.end());
			auto right_number      = right.numbers.begin();
			auto right_back_number = std::prev(right.numbers.end());
			while (true) {
				// addss this, right
				*this_number += *right_number;
				
				// carry
				if ( *this_number > 9 ) {
					uint8_t carry = *this_number / 10;
					*this_number  = *this_number % 10;

					auto this_next_number = std::next(this_number);
					while (carry != 0) {
						*this_next_number += carry;
						carry             = *this_next_number / 10;
						*this_next_number %= 10;
						++this_next_number;
					}
				}
				
				// break if back_number calculated
				if (this_number == this_back_number || right_number == right_back_number) {
					break;
				}

				++this_number;
				++right_number;
			}

			// pop Zero-precisions, {234321.0}
			while (numbers.size() > decimal_point+1 && numbers.back() == 0) {
				numbers.pop_back();
			}
			while (decimal_point != 0 && numbers.front() == 0) {
				numbers.pop_front();
				--decimal_point;
			}
		}

		void subev_without_sign(const decimal& right) {
			assert( compare_without_sign(right) == 1 );
			
			// push Zero-precision to lowest-bit(decimal)
			if ( this->decimal_point < right.decimal_point ) {
				size_t diff = right.decimal_point - this->decimal_point;
				
				for (size_t i = 0; i != diff; ++i) {
					this->numbers.push_front(0);
				}
				this->decimal_point += diff;
			}

			// subss this[i], right[i]
			auto this_back_number = std::next(this->numbers.begin(), this->decimal_point - right.decimal_point);
			auto this_number = std::next(this_back_number, right.numbers.size() - 1);
			auto right_number = std::prev(right.numbers.end());
			auto right_back_number = right.numbers.begin();
			while (true) {
				// borrow if ...
				if (*this_number < *right_number) {
					auto this_borrow = std::next(this_number);
					while (*this_borrow == 0) {
						*this_borrow++ = 9;
					}
					*this_borrow -= 1;
					*this_number += 10;
				}

				// subss this, right
				*this_number -= *right_number;
				
				// break if back_number calculated
				if (this_number == this_back_number || right_number == right_back_number) {
					break;
				}

				--this_number;
				--right_number;
			}

			// pop Zero-precisions, {234321.0}
			while (numbers.size() > decimal_point+1 && numbers.back() == 0) {
				numbers.pop_back();
			}
			while (decimal_point != 0 && numbers.front() == 0) {
				numbers.pop_front();
				--decimal_point;
			}
		}

		void mulev2() {
			// mulss this[i], 2
			std::for_each(numbers.begin(), numbers.end(), 
				[](uint8_t& number) { number *= 2; });
			
			// push a Zero-precision to highest-bit(integer)
			numbers.push_back(0);
			
			// carry
			auto this_number = numbers.begin();
			auto this_back_number = std::prev(numbers.end());
			while (true) {
				if ( *this_number > 9 ) {
					uint8_t carry = *this_number / 10;
					*this_number  = *this_number % 10;

					auto this_next_number = std::next(this_number);
					while (carry != 0) {
						*this_next_number += carry;
						carry             = *this_next_number / 10;
						*this_next_number %= 10;
						++this_next_number;
					}
				}

				if (this_number == this_back_number) {
					break;
				}

				++this_number;
			}

			// pop Zero-precisions, {234321.0}
			while (numbers.size() > decimal_point+1 && numbers.back() == 0) {
				numbers.pop_back();
			}
			while (decimal_point != 0 && numbers.front() == 0) {
				numbers.pop_front();
				--decimal_point;
			}
		}

		void divev2() {
			// push 5 Zero-precision to lowest-bit(decimal)
			size_t try_precision = 5;
			for (size_t i = 0; i != try_precision; ++i) {
				numbers.push_front(0);
			}
			decimal_point += 5;

			// divss this[i], 2
			auto this_number = std::prev(numbers.end());
			auto this_back_number = numbers.begin();
			uint8_t remainder = 0;
			while (true) {
				uint8_t divisor = remainder*10  + *this_number;
				remainder = divisor % 2;
				divisor = divisor / 2;
				*this_number = divisor;

				if (this_number == this_back_number) {
					break;
				}

				--this_number;
			}

			// pop Zero-precisions, {234321.0}
			while (numbers.size() > decimal_point+1 && numbers.back() == 0) {
				numbers.pop_back();
			}
			while (decimal_point != 0 && numbers.front() == 0) {
				numbers.pop_front();
				--decimal_point;
			}
		}

		void reset_precision(size_t precision = 15) {
			auto significant_last = std::find_if(numbers.rbegin(), numbers.rend(), 
				[](uint8_t number) { return number != uint8_t(0); }
				);
			if ( significant_last != numbers.rend() ) {
				size_t _Myprecision = 
					significant_last.base() < std::next(numbers.begin(), decimal_point) 
						? decimal_point
						: std::distance(numbers.begin(), significant_last.base());
				while (decimal_point != 0 && _Myprecision > precision) {
					if (_Myprecision - 1 == precision) {
						// round if ...
						if (numbers.front() >= 5) {

							numbers.push_back(0);
							auto next_number = std::next(numbers.begin());
							*next_number += 1;
							while (*next_number > 9) {
								*next_number = 0;
								++next_number;
								*next_number += 1;
							}
							if (numbers.back() == 0) {
								numbers.pop_back();
							}

						}
					}
					numbers.pop_front();
					--decimal_point;
					--_Myprecision;
				}
			}
		}

		decimal& operator+=(const decimal& right) {
			if ( negative == right.negative ) {
				this->addev_without_sign(right);
			} else {
				int cmp = this->compare_without_sign(right);
				if (cmp == 0) {
					*this = decimal();
				} else if (cmp == 1) {
					this->subev_without_sign(right);
				} else {
					auto temp = right;
					temp += *this;
					*this = std::move(temp);
				}
			}

			return *this;
		}

		decimal operator+(const decimal& right) const {
			return decimal(*this) += right;
		}

		decimal& operator-=(const decimal& right) {
			if ( negative != right.negative ) {// -5 - +3 = -5 - 3 = -(5 + 3)
				this->addev_without_sign(right);
			} else {
				int cmp = this->compare_without_sign(right);
				if (cmp == 0) {
					*this = decimal();
				} else if (cmp == 1) {// 5 - 3 = 2
					this->subev_without_sign(right);
				} else {// 3 - 5 = -(5-3)
					auto temp = right;
					temp -= *this;
					temp.negative ^= 1;
					*this = std::move(temp);
				}
			}

			return *this;
		}

		decimal operator-(const decimal& right) const {
			return decimal(*this) -= right;
		}

		// { small to larger, meaning 0,1,2,3,4,5,6... or back to front or forward, this is we style }
		std::string to_string(bool string_is_right_to_left = true) const {
			auto str = std::string(numbers.size() + 1, ' ');
			if (negative) {
				str.push_back('-');
			}

			auto str_dest = 
			std::transform(numbers.begin(), std::next(numbers.begin(),decimal_point), str.begin(), 
				[](uint8_t number) { return number + 48; });
			
			*str_dest++ = '.';
			
			std::transform(std::next(numbers.begin(),decimal_point), numbers.end(), str_dest,
				[](uint8_t number) { return number + 48; });
			
			if (string_is_right_to_left) {
				std::reverse(str.begin(), str.end());
			}

			return std::move(str);
		}
		
		explicit decimal(std::string str, bool string_is_right_to_left = true) {
			assert(!str.empty());
			if (string_is_right_to_left) {
				std::reverse(str.begin(), str.end());
			}
			
			bool is_negative = str.back() == '-';
			bool is_integer = str.find('.') == std::string::npos;

			size_t sign_length = size_t(is_negative) + size_t(!is_integer);
			size_t number_length = str.size() - sign_length;
			this->numbers.resize(number_length, 0);
			
			this->negative = is_negative;
			if ( !is_integer ) {
				this->decimal_point = str.find('.');
				auto number_dest = 
				std::transform(
					str.begin(), 
					std::next(str.begin(), decimal_point), 
					this->numbers.begin(), 
					[](char ch) { return uint8_t(ch - 48); }
					);
				std::transform(
					std::next(str.begin(), decimal_point + 1), 
					std::prev(str.end(), size_t(is_negative)),
					number_dest,
					[](char ch) { return uint8_t(ch - 48); }
					);
			} else {
				this->decimal_point = 0;
				std::transform(
					str.begin(), 
					std::prev(str.end(), size_t(is_negative)),
					this->numbers.begin(),
					[](char ch) { return uint8_t(ch - 48); }
					);
			}
		}
	};
}



#include <bitset>
#include <stdexcept>
namespace std {
	// { slower than <add> for 5 times }
	template<size_t _Bits>
	std::bitset<_Bits>& operator+=(std::bitset<_Bits>& _Left, std::bitset<_Bits> _Right) {
		auto _Carry = _Left & _Right;
		for ( ; _Carry.any(); _Carry = _Left & _Right) {
			_Left    = _Left ^ _Right;
			_Right   = _Carry << 1;
		}

		return (_Left ^= _Right);
		/*
			00110101                      
		+ 00000001                    
		= 00110110                     
		forward_bit = A & B
		forward_bit <<= 1
		result = A | forward_bit

			10101111 A
		+ 00110110 B
		= 11100101 Result
		 
			10011001 A ^ B        = A2
		+ 01001100 (A & B) << 1 = B2
		= 11100101
		  
			11010101 A2 ^ B2        = A3
		+ 00010000 (A2 & B2) << 1 = B3
		= 11100101

			11000101 A3 ^ B3        = A4
		+ 00100000 (A3 & B3) << 1 = B4
		= 11100101

			00000000 A4 & B4
			11100101 A4 ^ B4 = Result
		*/
	}

	// { slower than <sub> for 3.1 times }
	template<size_t _Bits>
	std::bitset<_Bits>& operator-=(std::bitset<_Bits>& _Left, std::bitset<_Bits> _Right) {
		for ( ; ; ) {
			_Left  = _Left ^ _Right;
			_Right = _Left & _Right;// (_Left^_Right)&_Right
			if ( _Right.none() ) {
				break;
			}
		#if _DEBUG
			if ( _Right.test(_Bits-1) ) {// _Right > _Left
				throw std::underflow_error("std::operator-=(std::bitset<...>&, std::bitset<...>)");
			}
		#endif
			_Right <<= 1;
		}

		return _Left;
		/*
			10100111 A
		  - 10001111 B
		  = 00011000 Result
		
			00101000 A ^ B         = A1
		  - 00010000 (A1 & B) << 1 = B1
		  = 24(10digit)

			00111000 A1 ^ B1        = A2
		  - 00100000 (A2 & B1) << 1 = B2
		  = 24(10digit)

			00011000 A2 ^ B2 = A3 = Result
			00000000 A3 & B2
		*/
		/*
			00000000 A
		  - 01001000 B
		  = -72
	   
			01001000 A ^ B         = A1
		  - 10010000 (A1 & B) << 1 = B1
		  = 72-144(10digit)

			11011000 A1 ^ B1        = A2
		  - 00100000 (A2 & B1) << 1 = B2
		  = 216-32(10digit) error
		*/
		/*
			A^B = A1, A1 must less A if A >= B
			A^B must eliminate leftmost-bit if A.leftmost-bit == B.leftmost-bit == 1, so B1.test(farleft) == 1 must B > A, 
		*/
	}

	template<size_t _Bits>
	std::bitset<_Bits>& operator*=(std::bitset<_Bits>& _Left, std::bitset<_Bits> _Right) {
		/*{
			      10101001 = A
		        * 10101111 = B
		------------------
		         +10101001 = A<<0
		 	    +10101001  = A<<1
			   +10101001   = A<<2
			  +10101001    = A<<3
		     +00000000
		    +10101001      = A<<5
		   +00000000
          +10101001        = A<<7
		 -----------------
		= 0111001110000111
		}*/
		std::bitset<_Bits> _Temp = _Left;
		_Left.reset();
		
		while (_Right.any()) {
			if (_Right.test(0)) {
				_Left += _Temp;
			}
		#if _DEBUG
			if ( _Temp.test(_Bits - 1) ) {
				throw std::overflow_error("std::operator*=(std::bitset<...>&, std::bitset<...>)");
			}
		#endif
			_Temp <<= 1;
			_Right >>= 1;
		}

		return _Left;
	}

	template<size_t _Bits> inline
	std::bitset<_Bits>& operator+=(std::bitset<_Bits>& _Left, const int _Right) noexcept {
		return _Left += std::bitset<_Bits>(_Right);
	}

	template<size_t _Bits> inline
	std::bitset<_Bits>& operator-=(std::bitset<_Bits>& _Left, const int _Right) noexcept {
		return _Left -= std::bitset<_Bits>(_Right);
	}

	template<size_t _Bits> inline
	std::bitset<_Bits> operator+(const std::bitset<_Bits>& _Left, const std::bitset<_Bits>& _Right) noexcept {
		std::bitset<_Bits> _Ans = _Left;
		return _Ans += _Right;
	}

	template<size_t _Bits> inline
	std::bitset<_Bits> operator-(const std::bitset<_Bits>& _Left, const std::bitset<_Bits>& _Right) noexcept {
		std::bitset<_Bits> _Ans = _Left;
		return _Ans -= _Right;
	}

	template<size_t _Bits> inline
	std::bitset<_Bits> operator*(const std::bitset<_Bits>& _Left, const std::bitset<_Bits>& _Right) noexcept {
		std::bitset<_Bits> _Ans = _Left;
		return _Ans *= _Right;
	}

	template<size_t _Bits> inline
	std::bitset<_Bits> operator+(const std::bitset<_Bits>& _Left, const int _Right) noexcept {
		std::bitset<_Bits> _Ans = _Left;
		return _Ans += std::bitset<_Bits>(_Right);
	}

	template<size_t _Bits> inline
	std::bitset<_Bits> operator-(const std::bitset<_Bits>& _Left, const int _Right) noexcept {
		std::bitset<_Bits> _Ans = _Left;
		return _Ans -= std::bitset<_Bits>(_Right);
	}

	template<size_t _Bits> inline
	bool operator<(const std::bitset<_Bits>& _Left, const std::bitset<_Bits>& _Right) {
		using _Word = typename std::bitset<_Bits>::_Ty;
		size_t _Words = sizeof(std::bitset<_Bits>) / sizeof(_Word);
		for (size_t _Wpos = _Words - 1; _Wpos != size_t(-1); --_Wpos) {
			if (_Left._Getword(_Wpos) != _Right._Getword(_Wpos)) {
				return _Left._Getword(_Wpos) < _Right._Getword(_Wpos);
			}
		}

		return false;
		//return _Right != _Left && (_Left - _Right == std::bitset<_Bits>(0));
	}

	template<size_t _Bits> inline
	bool operator>(const std::bitset<_Bits>& _Left, const std::bitset<_Bits>& _Right) {
		return _Right < _Left;
	}

	template<size_t _Bits> inline
	bool operator<=(const std::bitset<_Bits>& _Left, const std::bitset<_Bits>& _Right) {
		return !(_Left > _Right);
	}

	template<size_t _Bits> inline
	bool operator>=(const std::bitset<_Bits>& _Left, const std::bitset<_Bits>& _Right) {
		return !(_Left < _Right);
	}

	template<size_t _OutBits, size_t _InBits>
	std::bitset<_OutBits> bitset_cast(const std::bitset<_InBits>& _Source) {
	if constexpr (_OutBits < _InBits) {
		return reinterpret_cast<const std::bitset<_OutBits>&>(_Source);
	} else if constexpr(_InBits < _OutBits) {
		std::bitset<_OutBits> _Destination;
		std::copy(reinterpret_cast<const char*>(&_Source), 
			reinterpret_cast<const char*>(&_Source) + sizeof(_Source), 
			reinterpret_cast<char*>(&_Destination));
		return _Destination;
	} else {
		return _Source;
	}
}
}

#include <math.h>
namespace calculation {
	// { template float_<MantissaBits, ExponentBits>, (1 + 0.Mantissa) * 2^Exponent * (-1)^Sign }
	template<size_t Mn, size_t En, bool Isbase = false>
	class float_ {
	public:
		static constexpr size_t mantissa_bits = Mn;
		static constexpr size_t exponent_bits = En;
		static constexpr size_t sign_bits = 1;
		static constexpr size_t bits = En + Mn + 1;
		static_assert(exponent_bits != 0, "assert(exponent_bits != 0)");
		static_assert(mantissa_bits != 0, "assert(mantissa_bits != 0)");

		static constexpr size_t mantissa_offset_bits = 0;
		static constexpr size_t exponent_offset_bits = mantissa_offset_bits + mantissa_bits;
		static constexpr size_t sign_offset_bits = exponent_offset_bits + exponent_bits;
		
		static std::bitset<bits> sign_mask() {
			static auto _Mask = std::bitset<bits>(1) << (bits - 1);
			return _Mask;
		}
		
		static std::bitset<bits> exponent_mask() {
			static auto _Mask = ((std::bitset<bits>(1) << exponent_bits) - 1) << exponent_offset_bits;
			return _Mask;
		}
		
		static std::bitset<bits> mantissa_mask() {
			static auto _Mask = (std::bitset<bits>(1) << mantissa_bits) - 1 /* << 0 */;
			return _Mask;
		}

		static std::bitset<bits> hidden_significant() {
			static auto _Mask = std::bitset<bits>(1) << mantissa_bits;
			return _Mask;
		}

		static std::bitset<bits> exponent_bias() {
			static auto _Mask = ((std::bitset<bits>(1) << (exponent_bits - 1)) - 1) << exponent_offset_bits;
			return _Mask;
		}

		static std::bitset<bits> zero() {
			return std::bitset<bits>(0);
		}

		static std::bitset<bits> infinite() {
			return exponent_mask();
		}

	public:

		std::bitset<bits> _Mybitset;

		const std::bitset<bits>& bitset() const {
			return _Mybitset;
		}

		std::bitset<bits>& bitset() {
			return _Mybitset;
		}

		float_() = default;

		explicit 
		float_(std::bitset<bits> __bitset) : _Mybitset(__bitset) {}

		template<size_t En2, size_t Mn2> explicit 
		float_(const float_<En2, Mn2>& other) {
			using other_float = float_<En2, Mn2>;
			constexpr 
			size_t other_bits = other_float::bits;
			auto other_bitset = other.bitset();

			// _My_mantissa = shifted(other_mantissa)
			std::bitset<other_bits> other_mantissa = other_bitset & other_float::mantissa_mask();
			std::bitset<bits> _My_mantissa;
			if constexpr ( mantissa_bits < other_float::mantissa_bits ) {
				std::bitset<other_bits> shifted_mantissa = other_mantissa >> (other_float::mantissa_bits - mantissa_bits);
				std::bitset<bits> casted_shifted_mantissa = std::bitset_cast<bits>(shifted_mantissa);
				_My_mantissa = casted_shifted_mantissa;
			} else if constexpr ( other_float::mantissa_bits < mantissa_bits ) {
				std::bitset<bits> casted_mantissa = std::bitset_cast<bits>( other_mantissa );
				std::bitset<bits> shifted_casted_mantissa = casted_mantissa << (mantissa_bits - other_float::mantissa_bits);
				_My_mantissa = shifted_casted_mantissa;
			} else {
				_My_mantissa = std::bitset_cast<bits>( other_mantissa );
			}
			
			// decompose other_exponent
			std::bitset<other_bits> other_abs_exponent = (other_bitset & other_float::exponent_mask());
			bool other_exponent_sign;
			if ( other_abs_exponent < other_float::exponent_bias() ) {
				other_exponent_sign = true;//negative
				other_abs_exponent = other_float::exponent_bias() - other_abs_exponent;
			} else {
				other_exponent_sign = false;//negative
				other_abs_exponent -= other_float::exponent_bias();
			}
			
			// other_exponent is infinite ?, abandon an exponent-upper [-125,128] to [-125,127]
			std::bitset<bits> _My_exponent;
			std::bitset<other_bits> _My_abs_infinite;
			if constexpr ( exponent_offset_bits < other_float::exponent_offset_bits ) {
				std::bitset<other_bits> casted_infinite = std::bitset_cast<other_bits>( infinite() - exponent_bias() - 1 );
				std::bitset<other_bits> shifted_casted_infinite = casted_infinite << (other_float::exponent_offset_bits - exponent_offset_bits);
				_My_abs_infinite = shifted_casted_infinite;
			} else if constexpr ( other_float::exponent_offset_bits < exponent_offset_bits ) {
				std::bitset<bits> shifted_infinite = (infinite() - exponent_bias() - 1) >> (exponent_offset_bits - other_float::exponent_offset_bits);
				std::bitset<other_bits> casted_shifted_infinite = std::bitset_cast<other_bits>(shifted_infinite);
				_My_abs_infinite = casted_shifted_infinite;
			} else {
				_My_abs_infinite = reinterpret_cast<const std::bitset<other_bits>&>( infinite() - exponent_bias() - 1 );
			}

			if (_My_abs_infinite < other_abs_exponent) {
				_My_exponent = infinite();
			} else {

				// _My_exponent = shifted(other_exponent)
				if constexpr ( exponent_offset_bits < other_float::exponent_offset_bits ) {
					std::bitset<other_bits> shifted_exponent = other_abs_exponent >> (other_float::exponent_offset_bits - exponent_offset_bits);
					std::bitset<bits> casted_shifted_exponent = std::bitset_cast<bits>(shifted_exponent);
					_My_exponent = casted_shifted_exponent;
				} else if constexpr ( other_float::exponent_offset_bits < exponent_offset_bits ) {
					std::bitset<bits> casted_exponent = std::bitset_cast<bits>( other_abs_exponent );
					std::bitset<bits> shifted_casted_exponent = casted_exponent << (exponent_offset_bits - other_float::exponent_offset_bits);
					_My_exponent = shifted_casted_exponent;
				} else {
					_My_exponent = std::bitset_cast<bits>( other_abs_exponent );
				}

				if (other_exponent_sign) {
					_My_exponent = exponent_bias() - _My_exponent;
				} else {
					_My_exponent += exponent_bias();
				}

			}
			
			bool is_negative = (other_bitset & other_float::sign_mask()) 
				== (decltype(other_bitset)(1) << other_float::sign_offset_bits);
			auto _My_sign = std::bitset<bits>(is_negative) << sign_offset_bits;

			_Mybitset = _My_sign | _My_exponent | _My_mantissa;
			if (other_float::iszero(other)) {
				_Mybitset = zero();
			}
		}
		
		float_(float ohter) : float_(reinterpret_cast<const float_<23,8>&>(ohter)) {}
		
		float_(double ohter) : float_(reinterpret_cast<const float_<52,11>&>(ohter)) {}

		float_(unsigned int other) {
			if (other == 0U) {
				_Mybitset = zero();
				return;
			}

			// 000000000.00000000000010100100010 * pow(2,mantissa_bits)
			auto this_exponent = exponent_bias() + (std::bitset<bits>(mantissa_bits) << exponent_offset_bits);
			auto this_significant = std::bitset<bits>(other);

			//  normalize significant-bits and exponent-bits
			auto hidden_significant_mask = ~mantissa_mask();
			auto exponent_one = std::bitset<bits>(1) << exponent_offset_bits;
			if ( this_significant < hidden_significant() ) {
				while ( (this_significant & hidden_significant_mask) != hidden_significant() ) {
					this_significant <<= 1;
					this_exponent -= exponent_one;
				}
			} else {
				while ( (this_significant & hidden_significant_mask) != hidden_significant() ) {
					this_significant >>= 1;
					this_exponent += exponent_one;
				}
			}

			_Mybitset = this_exponent | (this_significant & mantissa_mask());
		}

		float_(unsigned long long other) {
			if (other == 0ULL) {
				_Mybitset = zero();
				return;
			}

			// 1expexpexp.mantissamantissa * pow(2,mantissa_bits)
			auto this_exponent = exponent_bias() + (std::bitset<bits>(mantissa_bits) << exponent_offset_bits);
			auto this_significant = std::bitset<bits>(other);

			//  normalize significant-bits and exponent-bits
			auto hidden_significant_mask = ~mantissa_mask();
			auto exponent_one = std::bitset<bits>(1) << exponent_offset_bits;
			if ( this_significant < hidden_significant() ) {
				while ( (this_significant & hidden_significant_mask) != hidden_significant() ) {
					this_significant <<= 1;
					this_exponent -= exponent_one;
				}
			} else {
				while ( (this_significant & hidden_significant_mask) != hidden_significant() ) {
					this_significant >>= 1;
					this_exponent += exponent_one;
				}
			}

			_Mybitset = this_exponent | (this_significant & mantissa_mask());
		}

		float_(int other) : float_(other < 0 ? static_cast<unsigned int>(-other) : static_cast<unsigned int>(other)) {
			_Mybitset |= std::bitset<bits>(other < 0) << sign_offset_bits;
		}

		float_(long long other) : float_(other < 0 ? static_cast<unsigned long long>(-other) : static_cast<unsigned long long>(other)) {
			_Mybitset |= std::bitset<bits>(other < 0) << sign_offset_bits;
		}
		
		explicit operator float() const {
			auto the_float = float_<23,8>(*this);
			return reinterpret_cast<const float&>(the_float);
		}

		explicit operator double() const {
			auto the_double = float_<52,11>(*this);
			return reinterpret_cast<const double&>(the_double);
		}

		explicit operator unsigned int() const {
			const auto zero_exponent = exponent_bias();

			// only fraction
			auto this_exponent = (_Mybitset & exponent_mask());
			if ( this_exponent < zero_exponent ) {
				return static_cast<unsigned int>(0);
			}

			// compute saved_bits, @floor(float_)
			auto exp2_mantissa_bits = std::bitset<bits>(mantissa_bits) << exponent_offset_bits;
			auto truncated_exponent = this_exponent - exp2_mantissa_bits;
			auto trunc_exponent = zero_exponent - truncated_exponent;
			size_t trunc_bits = (trunc_exponent >> exponent_offset_bits).to_ulong();
			size_t saved_bits = (mantissa_bits+1) - trunc_bits;
			if (saved_bits > sizeof(unsigned int) * 8) {
				throw std::overflow_error("float_<...>::operator unsigned int() const");
			}

			// smaller-shift to correct integer-bits
			auto this_significant = _Mybitset & mantissa_mask() | hidden_significant();
			this_significant >>= ((mantissa_bits+1) - saved_bits);

			auto dest = std::bitset_cast<sizeof(unsigned int)*8>(this_significant);
			return reinterpret_cast<const unsigned int&>(dest);
		}

		explicit operator unsigned long long() const {
			const auto zero_exponent = exponent_bias();

			// only fraction
			auto this_exponent = (_Mybitset & exponent_mask());
			if ( this_exponent < zero_exponent ) {
				return static_cast<unsigned long long>(0);
			}

			// compute saved_bits, @floor(float_)
			auto exp2_mantissa_bits = std::bitset<bits>(mantissa_bits) << exponent_offset_bits;
			auto truncated_exponent = this_exponent - exp2_mantissa_bits;
			auto trunc_exponent = zero_exponent - truncated_exponent;
			size_t trunc_bits = (trunc_exponent >> exponent_offset_bits).to_ulong();
			size_t saved_bits = (mantissa_bits+1) - trunc_bits;
			if (saved_bits > sizeof(unsigned int) * 8) {
				throw std::overflow_error("float_<...>::operator unsigned int() const");
			}

			// smaller-shift to correct integer-bits
			auto this_significant = _Mybitset & mantissa_mask() | hidden_significant();
			this_significant >>= ((mantissa_bits+1) - saved_bits);

			auto dest = std::bitset_cast<sizeof(unsigned long long)*8>(this_significant);
			return reinterpret_cast<const unsigned long long&>(dest);
		}
		
		explicit operator int() const {
			abort();
		}

		explicit operator long long() const {
			abort();
		}

		template<typename Ty> 
		const Ty& as() const {
			return *reinterpret_cast<const Ty*>(&_Mybitset);
		}

		template<typename Ty> 
		Ty& as() { 
			return *reinterpret_cast<Ty*>(&_Mybitset); 
		}

		static bool iszero(const float_& x) {
			return (x.bitset() & (exponent_mask() | mantissa_mask())) == zero();
		}

		static bool isinf(const float_& x) {
			return (x.bitset() & exponent_mask()) == exponent_mask();
		}

		bool operator==(const float_& right) const {
			return _Mybitset == right._Mybitset 
				|| (iszero(*this) && iszero(right));
		}

		bool operator!=(const float_& right) const {
			return !(*this == right);
		}

		bool operator<(const float_& right) const {
			if (*this == right) {
				return false;
			}

			if ((_Mybitset & sign_mask()) > (right.bitset() & sign_mask())) {
				return true;
			}

			if ((_Mybitset & (~sign_mask())) < (right.bitset() & (~sign_mask()))) {
				return true;
			}

			return false;
		}

		bool operator>(const float_& right) const {
			return right < *this;
		}
		
		bool operator<=(const float_& right) const {
			return !(*this > right);
		}

		bool operator>=(const float_& right) const {
			return !(*this < right);
		}

		float_ operator-() const {
			auto neg_sign = std::bitset<bits>(1) << sign_offset_bits;
			return float_{ _Mybitset | neg_sign };
		}

		float_ operator+(const float_& right) const {
		// { slower than <addss> 40 times }
			if (iszero(*this)) {
				return right;
			} 
			if (iszero(right)) {
				return *this;
			}

			// get sign-bits, exponent-bits, significant-bits
			auto this_sign = _Mybitset & sign_mask();
			auto right_sign = right.bitset() & sign_mask();
			auto this_exponent = (_Mybitset & exponent_mask());
			auto right_exponent = (right.bitset() & exponent_mask());
			auto this_significant = (_Mybitset & mantissa_mask()) | hidden_significant();
			auto right_significant = (right.bitset() & mantissa_mask()) | hidden_significant();
			const auto exp2_one = std::bitset<bits>(1) << exponent_offset_bits;

			// sync exponent-bits, to greater
			if ( this_exponent < right_exponent ) {
				do {
					this_significant >>= 1;
					this_exponent += exp2_one;
				} while (this_exponent != right_exponent);
			} else if (this_exponent != right_exponent) {
				do {
					right_significant >>= 1;
					right_exponent += exp2_one;
				} while (right_exponent != this_exponent);
			}

			// add significant-bits, update sign-bits
			if (this_sign == right_sign) {
				this_significant += right_significant;
			} else if (this_significant > right_significant) {
				this_significant -= right_significant;
			} else if (this_significant != right_significant) {
				this_significant = right_significant - this_significant;
				this_sign = right_sign;
			} else {// this + right = 0
				return float_{ zero() };
			}

			//  normalize significant-bits and exponent-bits
			if ( this_significant < hidden_significant() ) {
				while ( (this_significant & hidden_significant()).none() ) {
					this_significant <<= 1;
					this_exponent -= exp2_one;
				}
			} else {
				// hidden-significant && (exponent|sign-bits) == 0
				const auto highbit_mask = ~mantissa_mask();
				while ( (this_significant & highbit_mask) != hidden_significant() ) {
					this_significant >>= 1;
					this_exponent += exp2_one;
				}
			}

			return float_{ this_sign | this_exponent | (this_significant & mantissa_mask()) };
		}

		float_ operator-(const float_& right) const {
			return *this + (-right);
		}

		float_ operator*(const float_& right) const {
			if (iszero(*this) || iszero(right)) {
				return float_{ zero() };
			}

			// this_sign = this_sign * right_sign
			auto this_sign = _Mybitset & sign_mask();
			auto right_sign = right.bitset() & sign_mask();
			this_sign ^= right_sign;
			
			// this_exponent = this_exponent + right_exponent
			auto this_exponent = (_Mybitset & exponent_mask());
			auto right_exponent = (right.bitset() & exponent_mask());
			this_exponent += right_exponent;// sign bit, so no overflow
			this_exponent -= exponent_bias();

			// this_significant * exp2(mantissa_bits) = (this_significant * right_significant)
			auto this_significant = std::bitset_cast<bits * 2>( (_Mybitset & mantissa_mask()) | hidden_significant() );
			auto right_significant = std::bitset_cast<bits * 2>( (right.bitset() & mantissa_mask()) | hidden_significant() );
			this_significant *= right_significant;

			// normalize...
			const auto exp2_one = std::bitset<bits>(1) << exponent_offset_bits;
			auto hidden_significant_ex = std::bitset<bits * 2>( 1 ) << (mantissa_bits * 2);
			if ( this_significant < hidden_significant_ex ) {
				while ( (this_significant & hidden_significant_ex).none() ) {
					this_significant <<= 1;
					this_exponent -= exp2_one;
				}
			} else {
				// hidden-significant && (exponent|sign-bits) == 0
				const auto highbit_mask = ~(hidden_significant_ex - 1);
				while ( (this_significant & highbit_mask) != hidden_significant_ex ) {
					this_significant >>= 1;
					this_exponent += exp2_one;
				}
			}
			this_significant >>= mantissa_bits;
			
			return float_{ this_sign | this_exponent | (std::bitset_cast<bits>(this_significant) & mantissa_mask()) };
		}
		
		float_ operator/(const float_& right) const {
			auto this_sign = _Mybitset & sign_mask();
			auto right_sign = right.bitset() & sign_mask();
			this_sign ^= right_sign;
			if (iszero(*this)) {
				return float_{ zero() };
			}
			if ( iszero(right) ) {
				return float_{ this_sign | infinite() };
			}

			auto this_exponent = (_Mybitset & exponent_mask());
			auto right_exponent = (right.bitset() & exponent_mask());
			this_exponent += exponent_bias();// avoid underflow
			this_exponent -= right_exponent;

			auto this_significant = std::bitset<bits>(0);
			auto dividend = (_Mybitset & mantissa_mask()) | hidden_significant();
			auto divisor = (right.bitset() & mantissa_mask()) | hidden_significant();
			size_t offset = mantissa_bits;// contains hidden-significant
			while (true) {
				if (dividend >= divisor) {
					this_significant |= (std::bitset<bits>(1) << offset);
					dividend -= divisor;
				}

				if (--offset == size_t(-1)) {
					break;
				}

				dividend <<= 1;
			}

			//  normalize significant-bits and exponent-bits
			const auto exp2_one = std::bitset<bits>(1) << exponent_offset_bits;
			if ( this_significant < hidden_significant() ) {
				while ( (this_significant & hidden_significant()).none() ) {
					this_significant <<= 1;
					this_exponent -= exp2_one;
				}
			} else {
				// hidden-significant && (exponent|sign-bits) == 0
				const auto highbit_mask = ~mantissa_mask();
				while ( (this_significant & highbit_mask) != hidden_significant() ) {
					this_significant >>= 1;
					this_exponent += exp2_one;
				}
			}

			return float_{ this_sign | this_exponent | (this_significant & mantissa_mask()) };
		}
		
		float_ operator+=(const float_& right) {
			*this = *this + right;
			return *this;
		}

		float_ operator-=(const float_& right) {
			*this = *this - right;
			return *this;
		}

		float_ operator*=(const float_& right) {
			*this = *this * right;
			return *this;
		}

		float_ operator/=(const float_& right) {
			*this = *this / right;
			return *this;
		}

		// numeric functions

		static float_ abs(const float_& left) {
			const auto abs_mask = ~(sign_mask());
			return float_{ left.bitset() & abs_mask };
		}

		static float_ floor(const float_& left) {
			const auto zero_exponent = exponent_bias();

			// only fraction, 1.010101010111... * exp2(0)
			auto left_exponent = left.bitset() & exponent_mask();
			if ( left_exponent < zero_exponent ) {
				return float_();
			}

			// check only integer, 1010101010111... * exp2(exponent - mantissa_bits), 
			auto exp2_mantissa_bits = std::bitset<bits>(mantissa_bits) << exponent_offset_bits;
			auto truncated_exponent = left_exponent - exp2_mantissa_bits;
			if ( truncated_exponent >= zero_exponent ) {
				return left;
			}

			auto trunc_exponent = zero_exponent - truncated_exponent;
			size_t trunc_bits = (trunc_exponent >> exponent_offset_bits).to_ulong();
			auto trunc_mask = ~( ( std::bitset<bits>(1) << trunc_bits ) - 1 );
			return float_{ left.bitset() & trunc_mask };
		}

		static float_ fract(const float_& left) {
			const auto zero_exponent = exponent_bias();

			// 1.010101010111... * exp2(0), check only fraction
			auto left_exponent = left.bitset() & exponent_mask();
			if ( left_exponent < zero_exponent ) {
				return left;
			}

			// 1010101010111... * exp2(exponent - mantissa_bits), check only integer
			auto exp2_mantissa_bits = std::bitset<bits>(mantissa_bits) << exponent_offset_bits;
			auto truncated_exponent = left_exponent - exp2_mantissa_bits;
			if ( truncated_exponent >= zero_exponent ) {
				return float_();
			}

			auto trunc_exponent = zero_exponent - truncated_exponent;
			size_t trunc_bits = (trunc_exponent >> exponent_offset_bits).to_ulong();
			auto fract_mask = (( std::bitset<bits>(1) << trunc_bits ) - 1);
			
			//  normalize significant-bits and exponent-bits, assert(larger_shift)
			auto left_significant = left.bitset() & fract_mask;
			auto exp2_one = std::bitset<bits>(1) << exponent_offset_bits;
			assert ( left_significant < hidden_significant() );
			while ( (left_significant & hidden_significant()) != hidden_significant() ) {
				left_significant <<= 1;
				left_exponent -= exp2_one;
			}
			return float_{ (left.bitset() & sign_mask()) | left_exponent | (left_significant & mantissa_mask()) };
		}
		
		static float_ ceil(const float_& left) {
			abort();
		}

		static float_ round(const float_& left) {
			abort();
		}
		
		static float_ exp(const float_& left) {
			abort();
		}

		static float_ log(const float_& left) {
			abort();
		}

		static float_ sqrt(const float_& left) {
			abort();
		}

		static float_ cbrt(const float_& left) {
			abort();
		}

		static float_ pow(const float_& left, const float_& right) {
			abort();
		}

		static float_ mod(const float_& left, const float_& right) {
			abort();
		}

		// period functoins

		static float_ sin(const float_& left) {
			abort();
		}

		static float_ cos(const float_& left) {
			abort();
		}

		static float_ tan(const float_& left) {
			abort();
		}

		static float_ asin(const float_& left) {
			abort();
		}
		
		static float_ acos(const float_& left) {
			abort();
		}

		static float_ atan(const float_& left) {
			abort();
		}
	};

	// { IEEE754 single-precision }
	template<>
	class float_<23,8> {
		using _Mybase = float_<23, 8, true>;
	public:
		static constexpr size_t mantissa_bits = 23;
		static constexpr size_t exponent_bits = 8;
		static constexpr size_t sign_bits = 1;
		static constexpr size_t bits = 32;

		static constexpr size_t mantissa_offset_bits = 0;
		static constexpr size_t exponent_offset_bits = mantissa_offset_bits + mantissa_bits;
		static constexpr size_t sign_offset_bits = exponent_offset_bits + exponent_bits;

		static constexpr std::bitset<bits> sign_mask() {
			return std::bitset<bits>(0b10000000000000000000000000000000);
		}
		static constexpr std::bitset<bits> exponent_mask() {
			return std::bitset<bits>(0b01111111100000000000000000000000);
		}
		static constexpr std::bitset<bits> mantissa_mask() {
			return std::bitset<bits>(0b00000000011111111111111111111111);
		}
		
		static constexpr std::bitset<bits> hidden_significant() {
			return std::bitset<bits>(0b00000000100000000000000000000000);
		}
		static constexpr std::bitset<bits> exponent_bias() {
			return std::bitset<bits>(0b00111111100000000000000000000000);
		}
		
		static constexpr std::bitset<bits> zero() {
			return std::bitset<bits>(0b00000000000000000000000000000000);
		}
		static constexpr std::bitset<bits> infinite() {
			return exponent_mask();
		}
	public:
		union {
			std::bitset<bits> _Mybitset;
			float _Myfp;
		};
		std::bitset<bits>& bitset() { return _Mybitset; }
		const std::bitset<bits>& bitset() const { return _Mybitset; }
		
		float_() : _Myfp() {}
		float_(const float_&) = default;
		float_(float_&&) = default;
		template<size_t Mn,size_t En>
		float_(const float_<Mn,En>& other) : _Mybitset(reinterpret_cast<const std::bitset<bits>&>( _Mybase(other) )) {}
		float_& operator=(const float_&) = default;

		constexpr float_(float other) : _Myfp(other) {}
		constexpr float_(double other) : float_(static_cast<float>(other)) {}
		constexpr float_(int other) : float_(static_cast<float>(other)) {}
		constexpr float_(long long other) : float_(static_cast<float>(other)) {}
		constexpr float_(unsigned int other) : float_(static_cast<float>(other)) {}
		constexpr float_(unsigned long long other) : float_(static_cast<float>(other)) {}
		constexpr operator float() const { return _Myfp; }
		constexpr explicit operator double() const { return static_cast<double>(this->operator float()); }
		constexpr explicit operator int() const { return static_cast<int>(this->operator float()); }
		constexpr explicit operator long long() const { return static_cast<long long>(this->operator float()); }
		constexpr explicit operator unsigned int() const { return static_cast<unsigned int>(this->operator float()); }
		constexpr explicit operator unsigned long long() const { return static_cast<unsigned long long>(this->operator float()); }

		inline constexpr bool operator==(float_ right) const {
			return _Myfp == right._Myfp;
		}
		inline constexpr bool operator!=(float_ right) const {
			return !(*this == right);
		}
		inline constexpr bool operator<(float_ right) const {
			return _Myfp < right._Myfp;
		}
		inline constexpr bool operator>(float_ right) const {
			return right < *this;
		}
		inline constexpr bool operator<=(float_ right) const {
			return !(*this > right);
		}
		inline constexpr bool operator>=(float_ right) const {
			return !(*this < right);
		}

		inline constexpr float_ operator-() const {
			return float_( -_Myfp );
		}
		inline constexpr float_ operator+(float_ right) const {
			return float_( _Myfp + right._Myfp );
		}
		inline constexpr float_ operator-(float_ right) const {
			return float_( _Myfp - right._Myfp );
		}
		inline constexpr float_ operator*(float_ right) const {
			return float_( _Myfp * right._Myfp );
		}
		inline constexpr float_ operator/(float_ right) const {
			return float_( _Myfp / right._Myfp );
		}

		inline float_ operator+=(float_ right) {
			reinterpret_cast<float&>(*this) += reinterpret_cast<const float&>(right);
			return *this;
		}
		inline float_ operator-=(float_ right) {
			reinterpret_cast<float&>(*this) -= reinterpret_cast<const float&>(right);
			return *this;
		}
		inline float_ operator*=(float_ right) {
			reinterpret_cast<float&>(*this) *= reinterpret_cast<const float&>(right);
			return *this;
		}
		inline float_ operator/=(float_ right) {
			reinterpret_cast<float&>(*this) /= reinterpret_cast<const float&>(right);
			return *this;
		}

		static inline bool iszero(float_ left) {
			return _Mybase::iszero(reinterpret_cast<const _Mybase&>(left));
		}
		static inline bool isinf(float_ left) {
			return _Mybase::isinf(reinterpret_cast<const _Mybase&>(left));
		}

		static inline float_ abs(float_ left) {
			return float_( _CSTD fabsf(reinterpret_cast<const float&>(left)) );
		}
		static inline float_ floor(float_ left) {
			return float_( _CSTD floorf(reinterpret_cast<const float&>(left)) );
		}
		static inline float_ fract(float_ left) {
			return float_( reinterpret_cast<const float&>(left) - _CSTD floorf(reinterpret_cast<const float&>(left)) );
		}
		static inline float_ ceil(float_ left) {
			return float_( _CSTD ceilf(reinterpret_cast<const float&>(left)) );
		}
		static inline float_ round(float_ left) {
			return float_( _CSTD roundf(reinterpret_cast<const float&>(left)) );
		}
		static inline float_ exp(float_ left) {
			return float_( _CSTD expf(reinterpret_cast<const float&>(left)) );
		}
		static inline float_ log(float_ left) {
			return float_( _CSTD logf(reinterpret_cast<const float&>(left)) );
		}
		static inline float_ sqrt(float_ left) {
			return float_( _CSTD sqrtf(reinterpret_cast<const float&>(left)) );
		}
		static inline float_ cbrt(float_ left) {
			return float_( _CSTD cbrtf(reinterpret_cast<const float&>(left)) );
		}
		static inline float_ pow(float_ left, float_ right) {
			return float_( _CSTD powf(reinterpret_cast<const float&>(left), reinterpret_cast<const float&>(right)) );
		}
		static inline float_ mod(float_ left, float_ right) {
			return float_( _CSTD fmodf(reinterpret_cast<const float&>(left), reinterpret_cast<const float&>(right)) );
		}
		
		static inline float_ sin(float_ left) {
			return float_( _CSTD sinf(reinterpret_cast<const float&>(left)) );
		}
		static inline float_ cos(float_ left) {
			return float_( _CSTD cosf(reinterpret_cast<const float&>(left)) );
		}
		static inline float_ tan(float_ left) {
			return float_( _CSTD tanf(reinterpret_cast<const float&>(left)) );
		}
		static inline float_ asin(float_ left) {
			return float_( _CSTD asinf(reinterpret_cast<const float&>(left)) );
		}
		static inline float_ acos(float_ left) {
			return float_( _CSTD acosf(reinterpret_cast<const float&>(left)) );
		}
		static inline float_ atan(float_ left) {
			return float_( _CSTD atanf(reinterpret_cast<const float&>(left)) );
		}
	};

	// { IEEE754 double-precision }
	template<>
	class float_<52,11> {
		using _Mybase = float_<52, 11, true>;
	public:
		static constexpr size_t mantissa_bits = 52;
		static constexpr size_t exponent_bits = 11;
		static constexpr size_t sign_bits = 1;
		static constexpr size_t bits = 64;

		static constexpr size_t mantissa_offset_bits = 0;
		static constexpr size_t exponent_offset_bits = mantissa_offset_bits + mantissa_bits;
		static constexpr size_t sign_offset_bits = exponent_offset_bits + exponent_bits;

		static constexpr std::bitset<bits> sign_mask() {
			return std::bitset<bits>(0b1000000000000000000000000000000000000000000000000000000000000000);
		}
		static constexpr std::bitset<bits> exponent_mask() {
			return std::bitset<bits>(0b0111111111110000000000000000000000000000000000000000000000000000);
		}
		static constexpr std::bitset<bits> mantissa_mask() {
			return std::bitset<bits>(0b0000000000001111111111111111111111111111111111111111111111111111);
		}
		
		static constexpr std::bitset<bits> hidden_significant() {
			return std::bitset<bits>(0b0000000000010000000000000000000000000000000000000000000000000000);
		}
		static constexpr std::bitset<bits> exponent_bias() {
			return std::bitset<bits>(0b0011111111110000000000000000000000000000000000000000000000000000);
		}
		
		static constexpr std::bitset<bits> zero() {
			return std::bitset<bits>(0b0000000000000000000000000000000000000000000000000000000000000000);
		}
		static constexpr std::bitset<bits> infinite() {
			return exponent_mask();
		}
	public:
		union {
			std::bitset<bits> _Mybitset;
			double _Myfp;
		};
		std::bitset<bits>& bitset() { return _Mybitset; }
		const std::bitset<bits>& bitset() const { return _Mybitset; }
		
		float_() : _Myfp() {}
		float_(const float_&) = default;
		float_(float_&&) = default;
		template<size_t Mn,size_t En>
		float_(const float_<Mn, En>& other) : _Mybitset(reinterpret_cast<const std::bitset<bits>&>( _Mybase(other) )) {}
		float_& operator=(const float_&) = default;

		constexpr float_(double other) : _Myfp(other) {}
		constexpr float_(float other) : float_(static_cast<double>(other)) {}
		constexpr float_(int other) : float_(static_cast<double>(other)) {}
		constexpr float_(long long other) : float_(static_cast<double>(other)) {}
		constexpr float_(unsigned int other) : float_(static_cast<double>(other)) {}
		constexpr float_(unsigned long long other) : float_(static_cast<double>(other)) {}
		constexpr operator double() const { return _Myfp; }
		constexpr explicit operator float() const { return  static_cast<float>(this->operator double()); }
		constexpr explicit operator int() const { return static_cast<int>(this->operator double()); }
		constexpr explicit operator long long() const { return static_cast<long long>(this->operator double()); }
		constexpr explicit operator unsigned int() const { return static_cast<unsigned int>(this->operator double()); }
		constexpr explicit operator unsigned long long() const { return static_cast<unsigned long long>(this->operator double()); }

		inline constexpr bool operator==(float_ right) const {
			return _Myfp == right._Myfp;
		}
		inline constexpr bool operator!=(float_ right) const {
			return !(*this == right);
		}
		inline constexpr bool operator<(float_ right) const {
			return _Myfp < right._Myfp;
		}
		inline constexpr bool operator>(float_ right) const {
			return right < *this;
		}
		inline constexpr bool operator<=(float_ right) const {
			return !(*this > right);
		}
		inline constexpr bool operator>=(float_ right) const {
			return !(*this < right);
		}

		inline constexpr float_ operator-() const {
			return float_( -_Myfp );
		}
		inline constexpr float_ operator+(float_ right) const {
			return float_( _Myfp + right._Myfp );
		}
		inline constexpr float_ operator-(float_ right) const {
			return float_( _Myfp - right._Myfp );
		}
		inline constexpr float_ operator*(float_ right) const {
			return float_( _Myfp * right._Myfp );
		}
		inline constexpr float_ operator/(float_ right) const {
			return float_( _Myfp / right._Myfp );
		}

		inline float_ operator+=(float_ right) {
			reinterpret_cast<double&>(*this) += reinterpret_cast<const double&>(right);
			return *this;
		}
		inline float_ operator-=(float_ right) {
			reinterpret_cast<double&>(*this) -= reinterpret_cast<const double&>(right);
			return *this;
		}
		inline float_ operator*=(float_ right) {
			reinterpret_cast<double&>(*this) *= reinterpret_cast<const double&>(right);
			return *this;
		}
		inline float_ operator/=(float_ right) {
			reinterpret_cast<double&>(*this) /= reinterpret_cast<const double&>(right);
			return *this;
		}

		static inline bool iszero(float_ left) {
			return _Mybase::iszero(reinterpret_cast<const _Mybase&>(left));
		}
		static inline bool isinf(float_ left) {
			return _Mybase::isinf(reinterpret_cast<const _Mybase&>(left));
		}

		static inline float_ abs(float_ left) {
			return float_( _CSTD fabs(reinterpret_cast<const double&>(left)) );
		}
		static inline float_ floor(float_ left) {
			return float_( _CSTD floor(reinterpret_cast<const double&>(left)) );
		}
		static inline float_ fract(float_ left) {
			return float_( reinterpret_cast<const double&>(left) - _CSTD floor(reinterpret_cast<const double&>(left)) );
		}
		static inline float_ ceil(float_ left) {
			return float_( _CSTD ceil(reinterpret_cast<const double&>(left)) );
		}
		static inline float_ round(float_ left) {
			return float_( _CSTD round(reinterpret_cast<const double&>(left)) );
		}
		static inline float_ exp(float_ left) {
			return float_( _CSTD exp(reinterpret_cast<const double&>(left)) );
		}
		static inline float_ log(float_ left) {
			return float_( _CSTD log(reinterpret_cast<const double&>(left)) );
		}
		static inline float_ sqrt(float_ left) {
			return float_( _CSTD sqrt(reinterpret_cast<const double&>(left)) );
		}
		static inline float_ cbrt(float_ left) {
			return float_( _CSTD cbrt(reinterpret_cast<const double&>(left)) );
		}
		static inline float_ pow(float_ left, float_ right) {
			return float_( _CSTD pow(reinterpret_cast<const double&>(left), reinterpret_cast<const double&>(right)) );
		}
		static inline float_ mod(float_ left, float_ right) {
			return float_( _CSTD fmod(reinterpret_cast<const double&>(left), reinterpret_cast<const double&>(right)) );
		}
		
		static inline float_ sin(float_ left) {
			return float_( _CSTD sin(reinterpret_cast<const double&>(left)) );
		}
		static inline float_ cos(float_ left) {
			return float_( _CSTD cos(reinterpret_cast<const double&>(left)) );
		}
		static inline float_ tan(float_ left) {
			return float_( _CSTD tan(reinterpret_cast<const double&>(left)) );
		}
		static inline float_ asin(float_ left) {
			return float_( _CSTD asin(reinterpret_cast<const double&>(left)) );
		}
		static inline float_ acos(float_ left) {
			return float_( _CSTD acos(reinterpret_cast<const double&>(left)) );
		}
		static inline float_ atan(float_ left) {
			return float_( _CSTD atan(reinterpret_cast<const double&>(left)) );
		}
	};

	// next_work: { XXXX quadruple-precision }


	template<size_t Mn, size_t En> inline
	float_<Mn,En> abs(const float_<Mn,En>& x) {
		return float_<Mn,En>::abs(x);
	}

	template<size_t Mn, size_t En> inline
	float_<Mn,En> floor(const float_<Mn,En>& x) {
		return float_<Mn,En>::floor(x);
	}
	
	template<size_t Mn, size_t En> inline
	float_<Mn,En> fract(const float_<Mn,En>& x) {
		return float_<Mn,En>::fract(x);
	}

	template<size_t Mn, size_t En> inline
	float_<Mn,En> ceil(const float_<Mn,En>& x) {
		return float_<Mn,En>::ceil(x);
	}

	template<size_t Mn, size_t En> inline
	float_<Mn,En> round(const float_<Mn,En>& x) {
		return float_<Mn,En>::round(x);
	}

	template<size_t Mn, size_t En> inline
	float_<Mn,En> exp(const float_<Mn,En>& x) {
		return float_<Mn,En>::exp(x);
	}

	template<size_t Mn, size_t En> inline
	float_<Mn,En> log(const float_<Mn,En>& x) {
		return float_<Mn,En>::log(x);
	}

	template<size_t Mn, size_t En> inline
	float_<Mn,En> sqrt(const float_<Mn,En>& x) {
		return float_<Mn,En>::sqrt(x);
	}

	template<size_t Mn, size_t En> inline
	float_<Mn,En> cbrt(const float_<Mn,En>& x) {
		return float_<Mn,En>::cbrt(x);
	}

	template<size_t Mn, size_t En> inline
	float_<Mn,En> pow(const float_<Mn,En>& left, const float_<Mn,En>& right) {
		return float_<Mn,En>::ceil(left, right);
	}

	template<size_t Mn, size_t En> inline
	float_<Mn,En> mod(const float_<Mn,En>& left, const float_<Mn,En>& right) {
		return float_<Mn,En>::mod(left, right);
	}

	/* Helper:{
		MathSymbol: { [int, float, double, ...], These should not be regarded as types, but as symbols },
		MathType: [float_<...>, integer_<...>, ration<...>, ...]
	} */

	template<size_t Mn, size_t En> inline
	float_<Mn,En> sin(const float_<Mn,En>& x) {
		return float_<Mn,En>::sin(x);
	}

	template<size_t Mn, size_t En> inline
	float_<Mn,En> cos(const float_<Mn,En>& x) {
		return float_<Mn,En>::cos(x);
	}

	template<size_t Mn, size_t En> inline
	float_<Mn,En> tan(const float_<Mn,En>& x) {
		return float_<Mn,En>::tan(x);
	}

	template<size_t Mn, size_t En> inline
	float_<Mn,En> asin(const float_<Mn,En>& x) {
		return float_<Mn,En>::asin(x);
	}

	template<size_t Mn, size_t En> inline
	float_<Mn,En> acos(const float_<Mn,En>& x) {
		return float_<Mn,En>::acos(x);
	}

	template<size_t Mn, size_t En> inline
	float_<Mn,En> atan(const float_<Mn,En>& x) {
		return float_<Mn,En>::atan(x);
	}

	// { dependence on calculation::decimal }
	template<size_t Mn, size_t En>
	std::string to_string(float_<Mn,En> _Source) {
		if (float_<Mn, En>::iszero(_Source)) {
			return "0.";
		} 
		if (float_<Mn, En>::isinf(_Source)) {
			return (_Source < float_<Mn,En>(0) ? "-inf" : "inf"); 
		}
		
		decimal _Dest = decimal("0.");

		// _Dest.numbers = significant * 2^(mantissa_bits)
		decimal _Exp2_n = decimal("1.");
		auto _Mymantissa = _Source.bitset() & _Source.mantissa_mask();
		for (size_t i = 0; i != _Source.mantissa_bits; ++i, _Exp2_n.mulev2()) {
			if (_Mymantissa.test(i)) {
				_Dest += _Exp2_n;
			}
		}
		_Dest += _Exp2_n;// hidden-bit

		// correct exponent
		auto _Zero = std::bitset<Mn+En+1>(0);
		auto _One = std::bitset<Mn+En+1>(1) << _Source.exponent_offset_bits;
		auto _Exponent = _Source.bitset() & _Source.exponent_mask();
		_Exponent -= std::bitset<Mn + En + 1>(_Source.mantissa_bits) << _Source.exponent_offset_bits;// errro if isinf
		if (_Exponent > _Source.exponent_bias()) {
			do {
				_Dest.mulev2();
				_Exponent -= _One;
			} while (_Exponent != _Source.exponent_bias());
		} else if (_Exponent < _Source.exponent_bias()) {
			do {
				_Dest.divev2();
				_Exponent += _One;
			} while (_Exponent != _Source.exponent_bias());
		}

		// set sign
		_Dest.negative = (_Source.bitset() & _Source.sign_mask()) != _Zero;

		// set precision
		size_t desired_significant_count = 
			static_cast<size_t>( ::floor(::log(2)/::log(10) * (_Source.mantissa_bits+1)) );
		_Dest.reset_precision( desired_significant_count );
		return _Dest.to_string();
	}

	template<size_t Mn, size_t En> inline
	std::ostream& operator<<(std::ostream& _Ostr, float_<Mn,En> _Fp) {
		return (_Ostr << to_string(_Fp));
	}

	// symbol operators
#define __float_and_symbol_compare_operator(SYMBOL) \
	template<size_t Mn, size_t En> inline           \
	bool operator==(const float_<Mn,En>& left, SYMBOL right) { \
		return left == float_<Mn,En>(right);        \
	}                                               \
	template<size_t Mn, size_t En> inline           \
	bool operator!=(const float_<Mn,En>& left, SYMBOL right) { \
		return left != float_<Mn,En>(right);        \
	}                                               \
	template<size_t Mn, size_t En> inline           \
	bool operator<(const float_<Mn,En>& left, SYMBOL right) { \
		return left < float_<Mn,En>(right);         \
	}                                               \
	template<size_t Mn, size_t En> inline           \
	bool operator>(const float_<Mn,En>& left, SYMBOL right) { \
		return left > float_<Mn,En>(right);         \
	}                                               \
	template<size_t Mn, size_t En> inline           \
	bool operator<=(const float_<Mn,En>& left, SYMBOL right) { \
		return left <= float_<Mn,En>(right);        \
	}                                               \
	template<size_t Mn, size_t En> inline           \
	bool operator>=(const float_<Mn,En>& left, SYMBOL right) { \
		return left >= float_<Mn,En>(right);        \
	}                                               \
	template<size_t Mn, size_t En> inline           \
	bool operator==(SYMBOL left, const float_<Mn,En>& right) { \
		return float_<Mn,En>(left) == right;        \
	}                                               \
	template<size_t Mn, size_t En> inline           \
	bool operator!=(SYMBOL left, const float_<Mn,En>& right) { \
		return float_<Mn,En>(left) != right;        \
	}                                               \
	template<size_t Mn, size_t En> inline           \
	bool operator<(SYMBOL left, const float_<Mn,En>& right) { \
		return float_<Mn,En>(left) < right;         \
	}                                               \
	template<size_t Mn, size_t En> inline           \
	bool operator>(SYMBOL left, const float_<Mn,En>& right) { \
		return float_<Mn,En>(left) > right;         \
	}                                               \
	template<size_t Mn, size_t En> inline           \
	bool operator<=(SYMBOL left, const float_<Mn,En>& right) { \
		return float_<Mn,En>(left) <= right;        \
	}                                               \
	template<size_t Mn, size_t En> inline           \
	bool operator>=(SYMBOL left, const float_<Mn,En>& right) { \
		return float_<Mn,En>(left) >= right;        \
	}

	__float_and_symbol_compare_operator(int)
	__float_and_symbol_compare_operator(long long)
	__float_and_symbol_compare_operator(unsigned int)
	__float_and_symbol_compare_operator(unsigned long long)
	__float_and_symbol_compare_operator(float)
	__float_and_symbol_compare_operator(double)

#define __float_and_symbol_arithmetic_operator(SYMBOL) \
	template<size_t Mn, size_t En> inline              \
	float_<Mn,En> operator+(const float_<Mn,En>& left, SYMBOL right) { \
		return left + float_<Mn,En>(right);            \
	}                                                  \
	template<size_t Mn, size_t En> inline              \
	float_<Mn,En> operator-(const float_<Mn,En>& left, SYMBOL right) { \
		return left - float_<Mn,En>(right);            \
	}                                                  \
	template<size_t Mn, size_t En> inline              \
	float_<Mn,En> operator*(const float_<Mn,En>& left, SYMBOL right) { \
		return left * float_<Mn,En>(right);            \
	}                                                  \
	template<size_t Mn, size_t En> inline              \
	float_<Mn,En> operator/(const float_<Mn,En>& left, SYMBOL right) { \
		return left / float_<Mn,En>(right);            \
	}                                                  \
	template<size_t Mn, size_t En> inline              \
	float_<Mn,En> operator+(SYMBOL left, const float_<Mn,En>& right) { \
		return float_<Mn,En>(left) / right;            \
	}                                                  \
	template<size_t Mn, size_t En> inline              \
	float_<Mn,En> operator-(SYMBOL left, const float_<Mn,En>& right) { \
		return float_<Mn,En>(left) - right;            \
	}                                                  \
	template<size_t Mn, size_t En> inline              \
	float_<Mn,En> operator*(SYMBOL left, const float_<Mn,En>& right) { \
		return float_<Mn,En>(left) * right;            \
	}                                                  \
	template<size_t Mn, size_t En> inline              \
	float_<Mn,En> operator/(SYMBOL left, const float_<Mn,En>& right) { \
		return float_<Mn,En>(left) / right;            \
	}

	__float_and_symbol_arithmetic_operator(int)
	__float_and_symbol_arithmetic_operator(long long)
	__float_and_symbol_arithmetic_operator(unsigned int)
	__float_and_symbol_arithmetic_operator(unsigned long long)
	__float_and_symbol_arithmetic_operator(float)
	__float_and_symbol_arithmetic_operator(double)


	/* Helper:{
		MathSymbol: { [int, float, double, ...], These should not be regarded as types, but as symbols },
		MathType: [float_<...>, integer_<...>, ration<...>, ...]
	} */

	// { (1 + 0.Mantissa) * 2^Exponent * (-1)^Sign }
	using float32 = float_<23,8>;
	
	// { (1 + 0.Mantissa) * 2^Exponent * (-1)^Sign }
	using float64 = float_<52,11>;
	
	// { (1 + 0.Mantissa) * 2^Exponent * (-1)^Sign, GCC/quadmath/ALL }
	using float128 = float_<112,15>;

	// { ... }
	using float256 = float_<235,20>;

	// { ... }
	using float512 = float_<485,26>;

	/* Helper:{
		MathSymbol: { [int, float, double, ...], These should not be regarded as types, but as symbols },
		MathType: [float_<...>, integer_<...>, ration<...>, ...]
	} */
}// namespace float_<Mn,En>