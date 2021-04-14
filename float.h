#pragma once
/*{ "clmagic/calculation/fundamental/float":{
  "Author": "LongJiangnan",
  "Mail": "Jiang1998Nan@outlook.com",
  "Date": "2019-2021",
  "License": "Please identify Author",
  "Reference": [
    "https://www.rfwireless-world.com/Tutorials/floating-point-tutorial.html",
    "https://github.com/gcc-mirror/gcc/blob/master/gcc/real.c",
    "https://www.exploringbinary.com/binary-division/"
  ]
} }*/

#include <bitset>
#include <stdexcept>
namespace std {
  // { slower than <add> for 5 times }
  template<size_t _Bits>
  std::bitset<_Bits>& operator+=(std::bitset<_Bits>& _Left, std::bitset<_Bits> _Right) {
/*{
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
}*/
    auto _Carry = _Left & _Right;
    for ( ; _Carry.any(); _Carry = _Left & _Right) {
      _Left  = _Left ^ _Right;
      _Right = _Carry << 1;
    }

    return (_Left ^= _Right);
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
------------------
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

  template<size_t _Bits>
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
  /*infinite precision floating_point
  * digit: (1 + 0.Mantissa) * 2^Exponent * (-1)^Sign
  * __m: MantissaBits
  * __e: ExponentBits
  * __m + __e + 1 = Bits
  */
  template<size_t __m, size_t __e, bool Isbase = false>
  class floatX {
  public:
    static constexpr size_t mantissa_bits = __m;
    static constexpr size_t exponent_bits = __e;
    static constexpr size_t sign_bits = 1;
    static constexpr size_t bits = mantissa_bits + exponent_bits + sign_bits;
    static_assert(exponent_bits != 0, "assert(exponent_bits != 0)");
    static_assert(mantissa_bits != 0, "assert(mantissa_bits != 0)");

    static constexpr size_t mantissa_offset_bits = 0;
    static constexpr size_t exponent_offset_bits = mantissa_offset_bits + mantissa_bits;
    static constexpr size_t sign_offset_bits = exponent_offset_bits + exponent_bits;
		
    static std::bitset<bits> sign_mask() {
      static auto _Mask = std::bitset<bits>(1) << sign_offset_bits;
      return _Mask;
    }
		
    static std::bitset<bits> exponent_mask() {
      static auto _Mask = 
        ( (std::bitset<bits>(1) << exponent_bits) - 1 ) 
          << exponent_offset_bits;
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
      static auto _Exp2_zero = 
        ( (std::bitset<bits>(1) << (exponent_bits - 1)) - 1 )
          << exponent_offset_bits;
      return _Exp2_zero;
    }

    static std::bitset<bits> zero_mask() {
      return std::bitset<bits>(0);
    }

    static std::bitset<bits> infinity_mask() {
      // infinity = full 1 in exponent_bitset
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

    explicit floatX(std::bitset<bits> __bitset) : _Mybitset(__bitset) {}

    static floatX zero() {
      return floatX(zero_mask());
    }

    static floatX epsilon() noexcept {
      // epsilon = exp2(0 - mantissa_bits) 
      static floatX _Epsilon = floatX(
        exponent_bias()
        - (std::bitset<bits>(mantissa_bits) << exponent_offset_bits)
      );
      return _Epsilon;
    }

    static floatX infinity() {
      return floatX(infinity_mask());
    }

    static floatX quiet_NaN() {
      static floatX _Quiet_nan = floatX(
        exponent_mask() | (hidden_significant() >> 1)/* | 0 */
      );
      return _Quiet_nan;
    }
	
    static floatX signaling_NaN() {
      static floatX _Signaling_nan = floatX(
        exponent_mask() | (hidden_significant() >> 1) | std::bitset<bits>(1)
      );
      return _Signaling_nan;
    }

    static bool iszero(const floatX& x) {
      return (x.bitset() & (exponent_mask() | mantissa_mask())) == zero_mask();
    }

    static bool isinf(const floatX& x) {
      return (x.bitset() & exponent_mask()) == exponent_mask();
    }

  public:
    floatX() = default;
    
    /*{ "_Right_bitset to _This_bitset":{
      if(_This_bits < _Right_bits) {
        small first, next storage
      } else {
        storage first, next large
	  }
    }}*/
    template<size_t m2, size_t e2> explicit
    floatX(const floatX<m2, e2>& other) {
      if (other == 0) {
        // exponent is error, so ...
        _Mybitset = zero_mask();
        return;
      }

      using this_float = floatX<__m, __e>;
      using other_float = floatX<m2, e2>;
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
        other_exponent_sign = 1;//negative
        other_abs_exponent = other_float::exponent_bias() - other_abs_exponent;
      } else {
        other_exponent_sign = 0;//positive
        other_abs_exponent -= other_float::exponent_bias();
      }
			
      // other_exponent is infinity ?, abandon an exponent-upper [-125,128] to [-125,127]
      std::bitset<bits> _My_exponent;
      std::bitset<other_bits> _My_abs_infinite;
      if constexpr ( this_float::exponent_offset_bits < other_float::exponent_offset_bits ) {
        std::bitset<other_bits> casted_infinite = std::bitset_cast<other_bits>( this_float::infinity_mask() - this_float::exponent_bias() - 1 );
        std::bitset<other_bits> shifted_casted_infinite = casted_infinite << ( other_float::exponent_offset_bits - this_float::exponent_offset_bits );
        _My_abs_infinite = shifted_casted_infinite;
      } else if constexpr ( this_float::exponent_offset_bits > other_float::exponent_offset_bits ) {
        std::bitset<bits> shifted_infinite = (infinity_mask() - exponent_bias() - 1) >> (exponent_offset_bits - other_float::exponent_offset_bits);
        std::bitset<other_bits> casted_shifted_infinite = std::bitset_cast<other_bits>(shifted_infinite);
        _My_abs_infinite = casted_shifted_infinite;
      } else /*if ( other_float::exponent_offset_bits == this_float::exponent_offset_bits )*/ {
        _My_abs_infinite = std::bitset_cast<other_bits>( infinity_mask() - exponent_bias() - 1 );
      }

      if (_My_abs_infinite < other_abs_exponent) {
        _My_exponent = infinity_mask();
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
    }
		
    floatX(float ohter) : floatX(reinterpret_cast<const floatX<23,8>&>(ohter)) {}
		
    floatX(double ohter) : floatX(reinterpret_cast<const floatX<52,11>&>(ohter)) {}

    floatX(long double other) : floatX(static_cast<double>(other)) {}

    floatX(unsigned int other) {
      if (other == 0U) {
        _Mybitset = zero_mask();
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

    floatX(unsigned long long other) {
      if (other == 0ULL) {
        _Mybitset = zero_mask();
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

    floatX(bool other) {
      if( other ){
        _Mybitset = exponent_bias()/* | std::bitset<bits>(0)*/;
      } else {
        _Mybitset = zero_mask();
      }
    }

    floatX(int other) : floatX(other < 0 ? static_cast<unsigned int>(-other) : static_cast<unsigned int>(other)) {
      _Mybitset |= std::bitset<bits>(other < 0) << sign_offset_bits;
    }

    floatX(long long other) : floatX(other < 0 ? static_cast<unsigned long long>(-other) : static_cast<unsigned long long>(other)) {
      _Mybitset |= std::bitset<bits>(other < 0) << sign_offset_bits;
    }
		
    explicit operator float() const {
      auto the_float = floatX<23,8>(*this);
      return reinterpret_cast<const float&>(the_float);
    }

    explicit operator double() const {
      auto the_double = floatX<52,11>(*this);
      return reinterpret_cast<const double&>(the_double);
    }

    explicit operator long double() const {
      return static_cast<long double>(static_cast<double>(*this));
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
        throw std::overflow_error("floatX<...>::operator unsigned int() const");
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
      if (saved_bits > sizeof(unsigned long long) * 8) {
        throw std::overflow_error("floatX<...>::operator unsigned int() const");
      }

      // smaller-shift to correct integer-bits
      auto this_significant = _Mybitset & mantissa_mask() | hidden_significant();
      this_significant >>= ((mantissa_bits+1) - saved_bits);
       
      auto dest = std::bitset_cast<sizeof(unsigned long long)*8>(this_significant);
      return reinterpret_cast<const unsigned long long&>(dest);
    }
    
    explicit operator bool() const {
      return iszero(*this) ? false : true;
    }
    
    explicit operator int() const {
      return (_Mybitset & sign_mask()).none() 
        ? static_cast<int>(static_cast<unsigned int>(*this))
        : -static_cast<int>(static_cast<unsigned int>(*this));
    }

    explicit operator long long() const {
      return (_Mybitset & sign_mask()).none() 
        ? static_cast<long long>(static_cast<unsigned long long>(*this)) 
        : -static_cast<long long>(static_cast<unsigned long long>(*this));
    }

    template<typename Ty> 
    const Ty& as() const {
      return *reinterpret_cast<const Ty*>(&_Mybitset);
    }

    template<typename Ty> 
    Ty& as() { 
      return *reinterpret_cast<Ty*>(&_Mybitset); 
    }

    bool operator==(const floatX& right) const {
      return _Mybitset == right._Mybitset 
        || (iszero(*this) && iszero(right));
    }

    bool operator!=(const floatX& right) const {
      return !(*this == right);
    }

    bool operator<(const floatX& right) const {
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

    bool operator>(const floatX& right) const {
      return right < *this;
    }
		
    bool operator<=(const floatX& right) const {
      return !(*this > right);
    }

    bool operator>=(const floatX& right) const {
      return !(*this < right);
    }

    floatX operator-() const {
      return floatX{ ((_Mybitset & sign_mask()) ^ sign_mask()) | (_Mybitset & (~sign_mask())) };
    }

    floatX operator+(const floatX& right) const {
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
        return zero();
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

      return floatX{ this_sign | this_exponent | (this_significant & mantissa_mask()) };
    }

    floatX operator-(const floatX& right) const {
      return *this + (-right);
    }

    floatX operator*(const floatX& right) const {
      if (iszero(*this) || iszero(right)) {
        return zero();
      }

      // this_sign = this_sign * right_sign
      auto this_sign = _Mybitset & sign_mask();
      auto right_sign = right.bitset() & sign_mask();
      this_sign ^= right_sign;
			
      // this_exponent = this_exponent + right_exponent
      auto this_exponent = (_Mybitset & exponent_mask());
      auto right_exponent = (right.bitset() & exponent_mask());
      if (right_exponent > exponent_bias()) {
        this_exponent += (right_exponent - exponent_bias());
      } else if ( right_exponent < exponent_bias()) {
        this_exponent -= (exponent_bias() - right_exponent);
      }

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
			
      return floatX{ this_sign | this_exponent | (std::bitset_cast<bits>(this_significant) & mantissa_mask()) };
    }
		
    floatX operator/(const floatX& right) const {
      auto this_sign = _Mybitset & sign_mask();
      auto right_sign = right.bitset() & sign_mask();
      this_sign ^= right_sign;
      if (iszero(*this)) {
        return zero();
      }
      if ( iszero(right) ) {
        return floatX{ this_sign | infinity_mask() };
      }

      auto this_exponent = (_Mybitset & exponent_mask());
      auto right_exponent = (right.bitset() & exponent_mask());
      if (right_exponent > exponent_bias()) { 
        this_exponent -= (right_exponent - exponent_bias());
      } else if (right_exponent < exponent_bias()) { 
        this_exponent += (exponent_bias() - right_exponent);
      }

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
        const auto highbit_mask = ~mantissa_mask();
        while ( (this_significant & highbit_mask) != hidden_significant() ) {
          this_significant >>= 1;
          this_exponent += exp2_one;
        }
      }

      return floatX{ this_sign | this_exponent | (this_significant & mantissa_mask()) };
    }

    floatX operator%(const floatX& right) const {
      abort();
    }
		
    floatX& operator+=(const floatX& right) {
      *this = *this + right;
      return *this;
    }

    floatX& operator-=(const floatX& right) {
      *this = *this - right;
      return *this;
    }

    floatX& operator*=(const floatX& right) {
      *this = *this * right;
      return *this;
    }

    floatX& operator/=(const floatX& right) {
      *this = *this / right;
      return *this;
    }

    floatX& operator%=(const floatX& right) {
      abort();
    }

    static floatX abs(const floatX& left) {
      const auto abs_mask = ~(sign_mask());
      return floatX{ left.bitset() & abs_mask };
    }

    static floatX floor(const floatX& left) {
      const auto zero_exponent = exponent_bias();

      // only fraction, 1.010101010111... * exp2(0)
      auto left_exponent = left.bitset() & exponent_mask();
      if ( left_exponent < zero_exponent ) {
        return zero();
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
      return floatX{ left.bitset() & trunc_mask };
    }

    static floatX fract(const floatX& left) {
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
        return zero();
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
      return floatX{ (left.bitset() & sign_mask()) | left_exponent | (left_significant & mantissa_mask()) };
    }		
  };

	// { IEEE754 single-precision }
	template<>
	class floatX<23,8> {
		using _Mybase = floatX<23, 8, true>;
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
		static constexpr std::bitset<bits> zero_mask() {
			return std::bitset<bits>(0b00000000000000000000000000000000);
		}
		static constexpr std::bitset<bits> infinity_mask() {
			return exponent_mask();
		}
		
		static constexpr floatX zero() {
			return floatX{ 0.0f };
		}
		static constexpr floatX epsilon() {
			return floatX{ 1.192092896e-07F };
			// 0b00110100000000000000000000000000
		}
		static constexpr floatX infinity() {
			return floatX{ __builtin_huge_valf() };
			//0b01111111100000000000000000000000
		}
		static constexpr floatX quiet_NaN() {
			return floatX{ __builtin_nanf("0") };
			//0b01111111110000000000000000000000
		}
		static constexpr floatX signaling_NaN() {
			return floatX{ __builtin_nanf("1") };
			//0b01111111110000000000000000000001
		}
	public:
		float _Myfp;
		std::bitset<bits>& bitset() { return reinterpret_cast<std::bitset<bits>&>(*this); }
		const std::bitset<bits>& bitset() const { return reinterpret_cast<const std::bitset<bits>&>(*this); }
		
		floatX() = default;
		floatX(const floatX&) = default;
		floatX(floatX&&) = default;
		template<size_t m,size_t e>
		floatX(const floatX<m,e>& other) : _Myfp(reinterpret_cast<const float&>( _Mybase(other) )) {}
		floatX(const std::bitset<bits>& other) : _Myfp(reinterpret_cast<const float&>( other )) {}
		floatX& operator=(const floatX&) = default;

		constexpr floatX(float other) : _Myfp(other) {}
		constexpr floatX(double other) : floatX(static_cast<float>(other)) {}
		constexpr floatX(long double other) : floatX(static_cast<float>(other)) {}
		constexpr floatX(bool other) : floatX(static_cast<float>(other)) {}
		constexpr floatX(int other) : floatX(static_cast<float>(other)) {}
		constexpr floatX(long long other) : floatX(static_cast<float>(other)) {}
		constexpr floatX(unsigned int other) : floatX(static_cast<float>(other)) {}
		constexpr floatX(unsigned long long other) : floatX(static_cast<float>(other)) {}
		constexpr operator float() const { return _Myfp; }
		constexpr explicit operator double() const { return static_cast<double>(this->operator float()); }
		constexpr explicit operator long double() const { return static_cast<long double>(this->operator float()); }
		constexpr explicit operator bool() const { return static_cast<bool>(this->operator float()); }
		constexpr explicit operator int() const { return static_cast<int>(this->operator float()); }
		constexpr explicit operator long long() const { return static_cast<long long>(this->operator float()); }
		constexpr explicit operator unsigned int() const { return static_cast<unsigned int>(this->operator float()); }
		constexpr explicit operator unsigned long long() const { return static_cast<unsigned long long>(this->operator float()); }

		inline constexpr bool operator==(floatX right) const {
			return _Myfp == right._Myfp;
		}
		inline constexpr bool operator!=(floatX right) const {
			return _Myfp != right._Myfp;
		}
		inline constexpr bool operator<(floatX right) const {
			return _Myfp < right._Myfp;
		}
		inline constexpr bool operator>(floatX right) const {
			return _Myfp > right._Myfp;
		}
		inline constexpr bool operator<=(floatX right) const {
			return _Myfp <= right._Myfp;
		}
		inline constexpr bool operator>=(floatX right) const {
			return _Myfp >= right._Myfp;
		}

		inline constexpr floatX operator-() const {
			return floatX{ -_Myfp };
		}
		inline constexpr floatX operator+(floatX right) const {
			return floatX{ _Myfp + right._Myfp };
		}
		inline constexpr floatX operator-(floatX right) const {
			return floatX{ _Myfp - right._Myfp };
		}
		inline constexpr floatX operator*(floatX right) const {
			return floatX{ _Myfp * right._Myfp };
		}
		inline constexpr floatX operator/(floatX right) const {
			return floatX{ _Myfp / right._Myfp };
		}
		inline floatX operator%(floatX right) const {
			return floatX{ _Myfp - _CSTD floor(_Myfp / right._Myfp) * right._Myfp };
			//return _CSTD fmodf(_Myfp, right._Myfp);
		}

		inline floatX& operator+=(floatX right) {
			_Myfp += right._Myfp;
			return *this;
		}
		inline floatX& operator-=(floatX right) {
			_Myfp -= right._Myfp;
			return *this;
		}
		inline floatX& operator*=(floatX right) {
			_Myfp *= right._Myfp;
			return *this;
		}
		inline floatX& operator/=(floatX right) {
			_Myfp /= right._Myfp;
			return *this;
		}
		inline floatX& operator%=(floatX right) {
			_Myfp = _CSTD fmodf(_Myfp, right._Myfp);
			return *this;
		}
	};

	// { IEEE754 double-precision }
	template<>
	class floatX<52,11> {
		using _Mybase = floatX<52, 11, true>;
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
		static constexpr std::bitset<bits> zero_mask() {
			return std::bitset<bits>(0b0000000000000000000000000000000000000000000000000000000000000000);
		}
		static constexpr std::bitset<bits> infinity_mask() {
			return exponent_mask();
		}
	
		static constexpr floatX zero() {
			return floatX(0.0);
		}
		static constexpr floatX epsilon() noexcept {
			return floatX(2.2204460492503131e-016);
			//0b0011110010110000000000000000000000000000000000000000000000000000
		}
		static constexpr floatX infinity() {
			return __builtin_huge_val();
			//0b0111111111110000000000000000000000000000000000000000000000000000
		}
		static constexpr floatX quiet_NaN() {
			return floatX(__builtin_nan("0"));
			//0b0111111111111000000000000000000000000000000000000000000000000000
		}
		static constexpr floatX signaling_NaN() {
			return floatX(__builtin_nan("1"));
			//0b0111111111111000000000000000000000000000000000000000000000000001
		}
	public:
		double _Myfp;
		std::bitset<bits>& bitset() { return reinterpret_cast<std::bitset<bits>&>(*this); }
		const std::bitset<bits>& bitset() const { return reinterpret_cast<const std::bitset<bits>&>(*this); }
		
		floatX() = default;
		floatX(const floatX&) = default;
		floatX(floatX&&) = default;
		template<size_t m,size_t e>
		floatX(const floatX<m,e>& other) : _Myfp(reinterpret_cast<const double&>( _Mybase(other) )) {}
		floatX(const std::bitset<bits>& other) : _Myfp(reinterpret_cast<const double&>( other )) {}
		floatX& operator=(const floatX&) = default;

		constexpr floatX(double other) : _Myfp(other) {}
		constexpr floatX(float other) : floatX(static_cast<double>(other)) {}
		constexpr floatX(long double other) : floatX(static_cast<double>(other)) {}
		constexpr floatX(bool other) : floatX(static_cast<double>(other)) {}
		constexpr floatX(int other) : floatX(static_cast<double>(other)) {}
		constexpr floatX(long long other) : floatX(static_cast<double>(other)) {}
		constexpr floatX(unsigned int other) : floatX(static_cast<double>(other)) {}
		constexpr floatX(unsigned long long other) : floatX(static_cast<double>(other)) {}
		constexpr operator double() const { return _Myfp; }
		constexpr explicit operator float() const { return  static_cast<float>(this->operator double()); }
		constexpr explicit operator long double() const { return static_cast<long double>(this->operator double()); }
		constexpr explicit operator bool() const { return  static_cast<bool>(this->operator double()); }
		constexpr explicit operator int() const { return static_cast<int>(this->operator double()); }
		constexpr explicit operator long long() const { return static_cast<long long>(this->operator double()); }
		constexpr explicit operator unsigned int() const { return static_cast<unsigned int>(this->operator double()); }
		constexpr explicit operator unsigned long long() const { return static_cast<unsigned long long>(this->operator double()); }

		inline constexpr bool operator==(floatX right) const {
			return _Myfp == right._Myfp;
		}
		inline constexpr bool operator!=(floatX right) const {
			return _Myfp != right._Myfp;
		}
		inline constexpr bool operator<(floatX right) const {
			return _Myfp < right._Myfp;
		}
		inline constexpr bool operator>(floatX right) const {
			return _Myfp > right._Myfp;
		}
		inline constexpr bool operator<=(floatX right) const {
			return _Myfp <= right._Myfp;
		}
		inline constexpr bool operator>=(floatX right) const {
			return _Myfp >= right._Myfp;
		}

		inline constexpr floatX operator-() const {
			return floatX( -_Myfp );
		}
		inline constexpr floatX operator+(floatX right) const {
			return floatX( _Myfp + right._Myfp );
		}
		inline constexpr floatX operator-(floatX right) const {
			return floatX( _Myfp - right._Myfp );
		}
		inline constexpr floatX operator*(floatX right) const {
			return floatX( _Myfp * right._Myfp );
		}
		inline constexpr floatX operator/(floatX right) const {
			return floatX( _Myfp / right._Myfp );
		}
		inline floatX operator%(floatX right) const {
			return _CSTD fmod(_Myfp, right._Myfp);
		}

		inline floatX& operator+=(floatX right) {
			_Myfp += right._Myfp;
			return *this;
		}
		inline floatX& operator-=(floatX right) {
			_Myfp -= right._Myfp;
			return *this;
		}
		inline floatX& operator*=(floatX right) {
			_Myfp *= right._Myfp;
			return *this;
		}
		inline floatX& operator/=(floatX right) {
			_Myfp /= right._Myfp;
			return *this;
		}
		inline floatX& operator%=(floatX right) {
			_Myfp = _CSTD fmod(_Myfp, right._Myfp);
			return *this;
		}
	};

	// next_work: { XXXX quadruple-precision-templatespectialation }

#define __calculation_float_operator_with_literal(_OP_, _LITERAL_TYPE_)           \
	template<size_t m, size_t e> inline                                           \
	floatX<m,e> operator##_OP_##(const floatX<m,e>& left, _LITERAL_TYPE_ right) { \
		return left _OP_ static_cast< floatX<m,e> >(right);                       \
	}   

#define __calculation_float_operator_with_literal2(_OP_, _LITERAL_TYPE_)          \
	template<size_t m, size_t e> inline                                           \
	floatX<m,e> operator##_OP_##(_LITERAL_TYPE_ left, const floatX<m,e>& right) { \
		return static_cast< floatX<m,e> >(left) _OP_ right;                       \
	}   

#define __calculation_float_operator_with_literal_commutatibity(_OP_, _LITERAL_TYPE_) \
	template<size_t m, size_t e> inline                                           \
	floatX<m,e> operator##_OP_##(const floatX<m,e>& left, _LITERAL_TYPE_ right) { \
		return left _OP_ static_cast< floatX<m,e> >(right);                       \
	}                                                                             \
	template<size_t m, size_t e> inline                                           \
	floatX<m,e> operator##_OP_##(_LITERAL_TYPE_ left, const floatX<m,e>& right) { \
		return static_cast< floatX<m,e> >(left) _OP_ right;                       \
	} 

#define __calculation_float_lvalueoperator_with_literal(_OP_, _LITERAL_TYPE_)  \
	template<size_t m, size_t e> inline                                        \
	floatX<m,e>& operator##_OP_##(floatX<m,e>& left, _LITERAL_TYPE_ right) {   \
		return left _OP_ static_cast< floatX<m,e> >(right);                    \
	}   

#define __calculation_float_comparison_with_literal_commutatibity(_OP_, _LITERAL_TYPE_) \
	template<size_t m, size_t e> inline                                    \
	bool operator##_OP_##(const floatX<m,e>& left, _LITERAL_TYPE_ right) { \
		return left _OP_ static_cast< floatX<m,e> >(right);                \
	}                                                                      \
	template<size_t m, size_t e> inline                                    \
	bool operator##_OP_##(_LITERAL_TYPE_ left, const floatX<m,e>& right) { \
		return static_cast< floatX<m,e> >(left) _OP_ right;                \
	} 

	__calculation_float_lvalueoperator_with_literal(+=, float)
	__calculation_float_lvalueoperator_with_literal(+=, double)
	__calculation_float_lvalueoperator_with_literal(+=, long double)
	__calculation_float_lvalueoperator_with_literal(+=, bool)
	__calculation_float_lvalueoperator_with_literal(+=, int)
	__calculation_float_lvalueoperator_with_literal(+=, long long)
	__calculation_float_lvalueoperator_with_literal(+=, unsigned int)
	__calculation_float_lvalueoperator_with_literal(+=, unsigned long long)
	__calculation_float_operator_with_literal_commutatibity(+, float)
	__calculation_float_operator_with_literal_commutatibity(+, double)
	__calculation_float_operator_with_literal_commutatibity(+, long double)
	__calculation_float_operator_with_literal_commutatibity(+, bool)
	__calculation_float_operator_with_literal_commutatibity(+, int)
	__calculation_float_operator_with_literal_commutatibity(+, long long)
	__calculation_float_operator_with_literal_commutatibity(+, unsigned int)
	__calculation_float_operator_with_literal_commutatibity(+, unsigned long long)
	
	__calculation_float_lvalueoperator_with_literal(-=, float)
	__calculation_float_lvalueoperator_with_literal(-=, double)
	__calculation_float_lvalueoperator_with_literal(-=, long double)
	__calculation_float_lvalueoperator_with_literal(-=, bool)
	__calculation_float_lvalueoperator_with_literal(-=, int)
	__calculation_float_lvalueoperator_with_literal(-=, long long)
	__calculation_float_lvalueoperator_with_literal(-=, unsigned int)
	__calculation_float_lvalueoperator_with_literal(-=, unsigned long long)
	__calculation_float_operator_with_literal(-, float)
	__calculation_float_operator_with_literal2(-, float)
	__calculation_float_operator_with_literal(-, double)
	__calculation_float_operator_with_literal2(-, double)
	__calculation_float_operator_with_literal(-, long double)
	__calculation_float_operator_with_literal2(-, long double)
	__calculation_float_operator_with_literal(-, bool)
	__calculation_float_operator_with_literal2(-, bool)
	__calculation_float_operator_with_literal(-, int)
	__calculation_float_operator_with_literal2(-, int)
	__calculation_float_operator_with_literal(-, long long)
	__calculation_float_operator_with_literal2(-, long long)
	__calculation_float_operator_with_literal(-, unsigned int)
	__calculation_float_operator_with_literal2(-, unsigned int)
	__calculation_float_operator_with_literal(-, unsigned long long)
	__calculation_float_operator_with_literal2(-, unsigned long long)

	__calculation_float_lvalueoperator_with_literal(*=, float)
	__calculation_float_lvalueoperator_with_literal(*=, double)
	__calculation_float_lvalueoperator_with_literal(*=, long double)
	__calculation_float_lvalueoperator_with_literal(*=, bool)
	__calculation_float_lvalueoperator_with_literal(*=, int)
	__calculation_float_lvalueoperator_with_literal(*=, long long)
	__calculation_float_lvalueoperator_with_literal(*=, unsigned int)
	__calculation_float_lvalueoperator_with_literal(*=, unsigned long long)
	__calculation_float_operator_with_literal_commutatibity(*, float)
	__calculation_float_operator_with_literal_commutatibity(*, double)
	__calculation_float_operator_with_literal_commutatibity(*, long double)
	__calculation_float_operator_with_literal_commutatibity(*, bool)
	__calculation_float_operator_with_literal_commutatibity(*, int)
	__calculation_float_operator_with_literal_commutatibity(*, long long)
	__calculation_float_operator_with_literal_commutatibity(*, unsigned int)
	__calculation_float_operator_with_literal_commutatibity(*, unsigned long long)

	__calculation_float_lvalueoperator_with_literal(/=, float)
	__calculation_float_lvalueoperator_with_literal(/=, double)
	__calculation_float_lvalueoperator_with_literal(/=, long double)
	__calculation_float_lvalueoperator_with_literal(/=, bool)
	__calculation_float_lvalueoperator_with_literal(/=, int)
	__calculation_float_lvalueoperator_with_literal(/=, long long)
	__calculation_float_lvalueoperator_with_literal(/=, unsigned int)
	__calculation_float_lvalueoperator_with_literal(/=, unsigned long long)
	__calculation_float_operator_with_literal(/, float)
	__calculation_float_operator_with_literal2(/, float)
	__calculation_float_operator_with_literal(/, double)
	__calculation_float_operator_with_literal2(/, double)
	__calculation_float_operator_with_literal(/, long double)
	__calculation_float_operator_with_literal2(/, long double)
	__calculation_float_operator_with_literal(/, bool)
	__calculation_float_operator_with_literal2(/, bool)
	__calculation_float_operator_with_literal(/, int)
	__calculation_float_operator_with_literal2(/, int)
	__calculation_float_operator_with_literal(/, long long)
	__calculation_float_operator_with_literal2(/, long long)
	__calculation_float_operator_with_literal(/, unsigned int)
	__calculation_float_operator_with_literal2(/, unsigned int)
	__calculation_float_operator_with_literal(/, unsigned long long)
	__calculation_float_operator_with_literal2(/, unsigned long long)

	__calculation_float_lvalueoperator_with_literal(%=, float)
	__calculation_float_lvalueoperator_with_literal(%=, double)
	__calculation_float_lvalueoperator_with_literal(%=, long double)
	__calculation_float_lvalueoperator_with_literal(%=, bool)
	__calculation_float_lvalueoperator_with_literal(%=, int)
	__calculation_float_lvalueoperator_with_literal(%=, long long)
	__calculation_float_lvalueoperator_with_literal(%=, unsigned int)
	__calculation_float_lvalueoperator_with_literal(%=, unsigned long long)
	__calculation_float_operator_with_literal_commutatibity(%, float)
	__calculation_float_operator_with_literal_commutatibity(%, double)
	__calculation_float_operator_with_literal_commutatibity(%, long double)
	__calculation_float_operator_with_literal_commutatibity(%, bool)
	__calculation_float_operator_with_literal_commutatibity(%, int)
	__calculation_float_operator_with_literal_commutatibity(%, long long)
	__calculation_float_operator_with_literal_commutatibity(%, unsigned int)
	__calculation_float_operator_with_literal_commutatibity(%, unsigned long long)

	__calculation_float_comparison_with_literal_commutatibity(==, float)
	__calculation_float_comparison_with_literal_commutatibity(==, double)
	__calculation_float_comparison_with_literal_commutatibity(==, long double)
	__calculation_float_comparison_with_literal_commutatibity(==, bool)
	__calculation_float_comparison_with_literal_commutatibity(==, int)
	__calculation_float_comparison_with_literal_commutatibity(==, long long)
	__calculation_float_comparison_with_literal_commutatibity(==, unsigned int)
	__calculation_float_comparison_with_literal_commutatibity(==, unsigned long long)

	__calculation_float_comparison_with_literal_commutatibity(!=, float)
	__calculation_float_comparison_with_literal_commutatibity(!=, double)
	__calculation_float_comparison_with_literal_commutatibity(!=, long double)
	__calculation_float_comparison_with_literal_commutatibity(!=, bool)
	__calculation_float_comparison_with_literal_commutatibity(!=, int)
	__calculation_float_comparison_with_literal_commutatibity(!=, long long)
	__calculation_float_comparison_with_literal_commutatibity(!=, unsigned int)
	__calculation_float_comparison_with_literal_commutatibity(!=, unsigned long long)

	__calculation_float_comparison_with_literal_commutatibity(<, float)
	__calculation_float_comparison_with_literal_commutatibity(<, double)
	__calculation_float_comparison_with_literal_commutatibity(<, long double)
	__calculation_float_comparison_with_literal_commutatibity(<, bool)
	__calculation_float_comparison_with_literal_commutatibity(<, int)
	__calculation_float_comparison_with_literal_commutatibity(<, long long)
	__calculation_float_comparison_with_literal_commutatibity(<, unsigned int)
	__calculation_float_comparison_with_literal_commutatibity(<, unsigned long long)

	__calculation_float_comparison_with_literal_commutatibity(>, float)
	__calculation_float_comparison_with_literal_commutatibity(>, double)
	__calculation_float_comparison_with_literal_commutatibity(>, long double)
	__calculation_float_comparison_with_literal_commutatibity(>, bool)
	__calculation_float_comparison_with_literal_commutatibity(>, int)
	__calculation_float_comparison_with_literal_commutatibity(>, long long)
	__calculation_float_comparison_with_literal_commutatibity(>, unsigned int)
	__calculation_float_comparison_with_literal_commutatibity(>, unsigned long long)

	__calculation_float_comparison_with_literal_commutatibity(<=, float)
	__calculation_float_comparison_with_literal_commutatibity(<=, double)
	__calculation_float_comparison_with_literal_commutatibity(<=, long double)
	__calculation_float_comparison_with_literal_commutatibity(<=, bool)
	__calculation_float_comparison_with_literal_commutatibity(<=, int)
	__calculation_float_comparison_with_literal_commutatibity(<=, long long)
	__calculation_float_comparison_with_literal_commutatibity(<=, unsigned int)
	__calculation_float_comparison_with_literal_commutatibity(<=, unsigned long long)

	__calculation_float_comparison_with_literal_commutatibity(>=, float)
	__calculation_float_comparison_with_literal_commutatibity(>=, double)
	__calculation_float_comparison_with_literal_commutatibity(>=, long double)
	__calculation_float_comparison_with_literal_commutatibity(>=, bool)
	__calculation_float_comparison_with_literal_commutatibity(>=, int)
	__calculation_float_comparison_with_literal_commutatibity(>=, long long)
	__calculation_float_comparison_with_literal_commutatibity(>=, unsigned int)
	__calculation_float_comparison_with_literal_commutatibity(>=, unsigned long long)

#undef __calculation_float_operator_with_literal
#undef __calculation_float_operator_with_literal_commutatibity
#undef __calculation_float_lvalueoperator_with_literal
#undef __calculation_float_comparison_with_literal_commutatibity


	// { (1 + 0.Mantissa) * 2^Exponent * (-1)^Sign }
	using float32 = floatX<23,8>;
	
	// { (1 + 0.Mantissa) * 2^Exponent * (-1)^Sign }
	using float64 = floatX<52,11>;
	
	// { (1 + 0.Mantissa) * 2^Exponent * (-1)^Sign, GCC/quadmath/ALL }
	using float128 = floatX<112,15>;

	using float256 = floatX<235,20>;

	using float512 = floatX<485,26>;

	template<size_t m, size_t e> inline
	bool isinf(const floatX<m,e>& x) {
		return (x.bitset() & floatX<m,e>::exponent_mask()) == floatX<m,e>::infinity_mask();
	}
	inline bool isinf(float32 x) {
		return _CSTD isinf(x);
	}
	inline bool isinf(float64 x) {
		return _CSTD isinf(x);
	}

	template<size_t m, size_t e> inline
	bool isnan(const floatX<m,e>& x) {
		return (x.bitset() & floatX<m,e>::exponent_mask()) == floatX<m,e>::infinity_mask()
			&& (x.bitset() & floatX<m, e>::mantissa_mask()) != floatX<m, e>::zero_mask();
	}
	inline bool isnan(float32 x) {
		return _CSTD isnan(x);
	}
	inline bool isnan(float64 x) {
		return _CSTD isnan(x);
	}

	template<size_t m, size_t e> inline
	floatX<m,e> abs(const floatX<m,e>& x) {
		auto abs_mask = ~( floatX<m,e>::sign_mask() );
		return floatX<m,e>{ x.bitset() & abs_mask };
	}
	inline float32 abs(float32 x) {
		return float32{ _CSTD fabsf(x) };
	}
	inline float64 abs(float64 x) {
		return _CSTD fabs(x);
	}

	template<size_t m, size_t e> inline
	floatX<m,e> floor(const floatX<m,e>& x) {
		return floatX<m,e>::floor(x);
	}
	inline float32 floor(float32 x) {
		return float32{ _CSTD floorf(x) };
	}
	inline float64 floor(float64 x) {
		return _CSTD floor(x);
	}

	template<size_t m, size_t e> inline
	floatX<m,e> fract(const floatX<m,e>& x) {
		return floatX<m,e>::fract(x);
	}
	inline float32 fract(float32 x) {
		return x - _CSTD floorf(x);
	}
	inline float64 fract(float64 x) {
		return x -_CSTD floor(x);
	}

	template<size_t m, size_t e> inline
	floatX<m,e> ceil(const floatX<m,e>& x) {
		return floor(x) + static_cast<floatX<m,e>>(1);
	}
	inline float32 ceil(float32 x) {
		return float32{ _CSTD ceilf(x) };
	}
	inline float64 ceil(float64 x) {
		return _CSTD ceil(x);
	}

	template<size_t m, size_t e> inline
	floatX<m,e> round(const floatX<m,e>& x) {
		return floor(x + static_cast<floatX<m,e>>(0.5));
	}
	inline float32 round(float32 x) {
		return float32{ _CSTD roundf(x) };
	}
	inline float64 round(float64 x) {
		return _CSTD round(x);
	}

	template<size_t m, size_t e>
	floatX<m,e> pow(const floatX<m,e>& x, int power) {
		if (power > 1) {
			floatX<m,e> result = pow(x, power >> 1);
            result *= result;
                
            if (power & 1) {
                result *= x; // n odd
            }
                
            return result;
        } else if (power == 1) {
            return x;
        } else if (power == 0) {
            return 1;
        } else  /* power < 0 */ {
            return pow(1/x,-power);
        }
	}
	inline float32 pow(float32 x, int power) {
		if (power == 2) {
			return x * x;
		}
		return float32{ _CSTD powf(static_cast<float>(x), static_cast<float>(power)) };
	}
	inline float64 pow(float64 x, int power) {
		if (power == 2) {
			return x * x;
		}
		return _CSTD pow(x, static_cast<double>(power));
	}

	template<size_t m, size_t e> inline
	floatX<m,e> pow(const floatX<m,e>& left, const floatX<m,e>& right) {
		abort();
	}
	inline float32 pow(float32 x, float32 power) {
		return float32{ _CSTD powf(static_cast<float>(x), static_cast<float>(power)) };
	}
	inline float64 pow(float64 x, float64 power) {
		return _CSTD pow(x, power);
	}

	template<size_t m, size_t e> inline
	floatX<m,e> exp(const floatX<m,e>& x) {
		abort();
	}
	inline float32 exp(float32 x) {
		return float32{ _CSTD expf(static_cast<float>(x)) };
	}
	inline float64 exp(float64 x) {
		return _CSTD exp(x);
	}

	template<size_t m, size_t e> inline
	floatX<m,e> log(const floatX<m,e>& x) {
		abort();
	}
	inline float32 log(float32 x) {
		return float32{ _CSTD logf(static_cast<float>(x)) };
	}
	inline float64 log(float64 x) {
		return _CSTD log(x);
	}

	template<size_t m, size_t e>
	floatX<m,e> sqrt(const floatX<m,e>& x) {
		using Float = floatX<m,e>;
		if ( x > 0 ) {

			// (.../wiki/Methods_of_computing_square_roots)
			// Newton-iterate
			// Performance: 100000000[sqrt]/8550[ms] with float64 in Intel-i5-8500(3.00G[Hz])
			auto xi = x.bitset();
				// substract pow(2, mantissa_bits)
			xi -= std::bitset<Float::bits>(1) << Float::mantissa_bits;
				// divide by 2
			xi >>= 1;
				// add (exponent_bias + 1) / 2 * pow(2, mantissa_bits)
			xi += std::bitset<Float::bits>(1) << (Float::exponent_bits - 2 + Float::mantissa_bits);
			
			size_t counter = 100;// needs to be improved
			Float y = Float(xi);
			Float ym1;
			do {
				ym1 = y;
				y = (y + x / y) / 2;// y_next = y - (y-x/y)/2
			} while (--counter && abs(y-ym1) > Float::epsilon()*max(abs(y),Float::epsilon()));

			return std::move(y);
		
		} else if ( x == 0 ) {
			return Float::zero();
		} else   /* x <  0 */{
			return Float::quiet_NaN();
		}
	}
	inline float32 sqrt(float32 x) {
		return float32{ _CSTD sqrtf(x) };
	}
	inline float64 sqrt(float64 x) {
		return _CSTD sqrt(x);
	}

	template<size_t m, size_t e>
	floatX<m,e> cbrt(const floatX<m,e>& x) {
		using Float = floatX<m,e>;
		
		// (.../wiki/Cube_root)
		// Mathematician: Edmond-Halley
		// Performance: 100000000[cbrt]/10500[ms] with float64 in Intel-i5-8500(3.00G[Hz])
		size_t counter = 100;// needs to be improved
		Float y = x;
		Float factor;
		do {
			Float yyy = y * y * y;
			factor = (yyy + x + x) / (yyy + yyy + x);
			y *= factor;// y_next = y * (pow(y,3) + 2*x) / (2*pow(y,3) + x)
		} while (--counter && abs(factor - 1) > Float::epsilon());

		return std::move(y);
	}
	inline float32 cbrt(float32 x) {
		return float32{ _CSTD cbrtf(x) };
	}
	inline float64 cbrt(float64 x) {
		return _CSTD cbrt(x);
	}

	template<size_t m, size_t e> inline
	floatX<m,e> sin(const floatX<m,e>& x) {
		return static_cast<floatX<m,e>>(sin(static_cast<float64>(x)));
		//abort();
	}
	inline float32 sin(float32 x) {
		return float32{ _CSTD sinf(x) };
	}
	inline float64 sin(float64 x) {
		return _CSTD sin(x);
	}

	template<size_t m, size_t e> inline
	floatX<m,e> cos(const floatX<m,e>& x) {
		return static_cast<floatX<m,e>>(cos(static_cast<float64>(x)));
		//abort();
	}
	inline float32 cos(float32 x) {
		return float32{ _CSTD cosf(x) };
	}
	inline float64 cos(float64 x) {
		return _CSTD cos(x);
	}

	template<size_t m, size_t e> inline
	floatX<m,e> tan(const floatX<m,e>& x) {
		abort();
	}
	inline float32 tan(float32 x) {
		return float32{ _CSTD tanf(x) };
	}
	inline float64 tan(float64 x) {
		return _CSTD tan(x);
	}

	template<size_t m, size_t e> inline
	floatX<m,e> asin(const floatX<m,e>& x) {
		abort();
	}
	inline float32 asin(float32 x) {
		return float32{ _CSTD asinf(x) };
	}
	inline float64 asin(float64 x) {
		return _CSTD asin(x);
	}

	template<size_t m, size_t e> inline
	floatX<m,e> acos(const floatX<m,e>& x) {
		abort();
	}
	inline float32 acos(float32 x) {
		return float32{ _CSTD acosf(x) };
	}
	inline float64 acos(float64 x) {
		return _CSTD acos(x);
	}

	template<size_t m, size_t e> inline
	floatX<m,e> atan(const floatX<m,e>& x) {
		abort();
	}
	inline float32 atan(float32 x) {
		return float32{ _CSTD atanf(x) };
	}
	inline float64 atan(float64 x) {
		return _CSTD atan(x);
	}

	template<size_t m, size_t e> inline
	floatX<m,e> atan2(const floatX<m,e>& y, const floatX<m,e> x) {
		abort();
	}
	inline float32 atan2(float32 y, float32 x) {
		return float32{ _CSTD atan2f(y, x) };
	}
	inline float64 atan2(float64 y, float64 x) {
		return _CSTD atan2(y, x);
	}
}// namespace float_<m,e>

#include <limits>
namespace std {
	template<size_t m, size_t e>
	class numeric_limits<calculation::floatX<m,e>> : public _Num_float_base{ // numeric limits for arbitrary type _Ty (say little or nothing)
	public:
		_NODISCARD static constexpr calculation::floatX<m,e>(min)() noexcept {
			return calculation::floatX<m,e>();
		}

		_NODISCARD static calculation::floatX<m,e>(max)() noexcept {
			if constexpr (std::is_same_v<calculation::floatX<m,e>, calculation::float32>) {
				return std::numeric_limits<float>::max();
			}
			return calculation::floatX<m,e>();
		}

		_NODISCARD static constexpr calculation::floatX<m,e> lowest() noexcept {
			return calculation::floatX<m,e>();
		}

		_NODISCARD static calculation::floatX<m,e> epsilon() noexcept {
			return calculation::floatX<m,e>::epsilon();
		}

		_NODISCARD static constexpr calculation::floatX<m,e> round_error() noexcept {
			return calculation::floatX<m,e>();
		}

		_NODISCARD static constexpr calculation::floatX<m,e> denorm_min() noexcept {
			return calculation::floatX<m,e>();
		}

		_NODISCARD static calculation::floatX<m,e> infinity() noexcept {
			return calculation::floatX<m, e>::infinity();
		}

		_NODISCARD static calculation::floatX<m,e> quiet_NaN() noexcept {
			return calculation::floatX<m, e>::quiet_NaN();
		}

		_NODISCARD static calculation::floatX<m,e> signaling_NaN() noexcept {
			return calculation::floatX<m, e>::signaling_NaN();
		}
	};
}


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

	// { dependence on calculation::decimal }
	template<size_t m, size_t e>
	std::string to_string(floatX<m,e> _Source) {
		if (isnan(_Source)) {
			return "nan";
		}
		if (isinf(_Source)) {
			return (_Source < floatX<m,e>(0) ? "-inf" : "inf"); 
		}
		if (_Source == 0) {
			return "0.";
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
		using Float = floatX<m, e>;
		auto _Zero = std::bitset<Float::bits>(0);
		auto _One = std::bitset<Float::bits>(1) << _Source.exponent_offset_bits;
		auto _Exponent = _Source.bitset() & _Source.exponent_mask();
		_Exponent -= std::bitset<Float::bits>(_Source.mantissa_bits) << _Source.exponent_offset_bits;
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
			static_cast<size_t>( ::ceil(::log(2)/::log(10) * (_Source.mantissa_bits+1)) );
		_Dest.reset_precision( desired_significant_count );
		return _Dest.to_string();
	}

	template<size_t m, size_t e> inline
	std::ostream& operator<<(std::ostream& _Ostr, floatX<m,e> _Fp) {
		return (_Ostr << to_string(_Fp));
	}
}

/*Example:{
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
}*/