#pragma once
/*{ "clmagic/calculation/fundamental/float":{
  "Description": "infinite precision calculation"
  "License": "Please identify Author",
  "Author": "LongJiangnan",
  "Date": "2019-2021",
  "Mail": "Jiang1998Nan@outlook.com",
  "Reference": [
    "https://www.rfwireless-world.com/Tutorials/floating-point-tutorial.html",
    "https://www.slideshare.net/prochwani95/06-floating-point"
    "https://github.com/gcc-mirror/gcc/blob/master/gcc/real.c",
    "https://www.exploringbinary.com/binary-division/"
  ]
} }*/


#include <stdexcept>
#include <bitset>
namespace std 
  {
// { slower than <add> for 5 times }
template<size_t _Bits>
std::bitset<_Bits>& operator+=(std::bitset<_Bits>& _Left, std::bitset<_Bits> _Right) {
/**
 *   00110101                      
 * + 00000001                    
 * = 00110110                     
 * carry_bit = A & B
 * carry_bit <<= 1
 * result = A | carry_bit
 *
 *   10101111 A
 * + 00110110 B
 * = 11100101 Result
 *		 
 *   10011001 A ^ B        = A2
 * + 01001100 (A & B) << 1 = B2
 * = 11100101
 *
 *   11010101 A2 ^ B2        = A3
 * + 00010000 (A2 & B2) << 1 = B3
 * = 11100101
 *
 *   11000101 A3 ^ B3        = A4
 * + 00100000 (A3 & B3) << 1 = B4
 * = 11100101
 *
 *   00000000 A4 & B4
 *   11100101 A4 ^ B4 = Result
*/
  std::bitset<_Bits> _Carry = _Left & _Right;
  for ( ; _Carry.any(); _Carry = _Left & _Right) {
    _Left  = _Left ^ _Right;
    _Right = _Carry << 1;
  }

  return (_Left ^= _Right);
}

// { slower than <sub> for 3.1 times }
template<size_t _Bits>
std::bitset<_Bits>& operator-=(std::bitset<_Bits>& _Left, std::bitset<_Bits> _Right) {
/**
 *   10100111 A
 * - 10001111 B
 * = 00011000 Result
		
 *   00101000 A ^ B         = A1
 * - 00010000 (A1 & B) << 1 = B1
 * = 24(10digit)

 *   00111000 A1 ^ B1        = A2
 * - 00100000 (A2 & B1) << 1 = B2
 * = 24(10digit)

 *   00011000 A2 ^ B2 = A3 = Result
 *   00000000 A3 & B2
*/
/**
 *   00000000 A
 * - 01001000 B
 * = -72(10digit)
	   
 *   01001000 A ^ B         = A1
 * - 10010000 (A1 & B) << 1 = B1
 * = 72-144(10digit)

 *   11011000 A1 ^ B1        = A2
 * - 00100000 (A2 & B1) << 1 = B2
 * = 216-32(10digit) error
 * 
 * A^B = A1, A1 must less A if A >= B
 * A^B must eliminate leftmost-bit if A.leftmost-bit == B.leftmost-bit == 1, so B1.test(farleft) == 1 must B > A, 
*/
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
}

template<size_t _Bits>
std::bitset<_Bits>& operator*=(std::bitset<_Bits>& _Left, std::bitset<_Bits> _Right) {
/**
 *          10101001 = A
 *        * 10101111 = B
 * ------------------
 *         +10101001 = A<<0
 *        +10101001  = A<<1
 *       +10101001   = A<<2
 *      +10101001    = A<<3
 *     +00000000
 *    +10101001      = A<<5
 *   +00000000
 *  +10101001        = A<<7
 * ------------------
 * =0111001110000111
*/
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
std::bitset<_Bits>& operator+=(std::bitset<_Bits>& _Left, const unsigned int _Right) noexcept {
  return _Left += std::bitset<_Bits>(_Right);
}

template<size_t _Bits> inline
std::bitset<_Bits>& operator-=(std::bitset<_Bits>& _Left, const unsigned int _Right) noexcept {
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
std::bitset<_Bits> operator+(const std::bitset<_Bits>& _Left, const unsigned int _Right) noexcept {
  std::bitset<_Bits> _Ans = _Left;
  return _Ans += std::bitset<_Bits>(_Right);
}

template<size_t _Bits> inline
std::bitset<_Bits> operator-(const std::bitset<_Bits>& _Left, const unsigned int _Right) noexcept {
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

template<size_t _OutBits, size_t _InBits> inline
std::bitset<_OutBits> bitset_cast(const std::bitset<_InBits>& _Source) {
  if constexpr ( _OutBits == _InBits ) {
	  return _Source;
  } else if constexpr ( _OutBits < _InBits ) {
    return reinterpret_cast<const std::bitset<_OutBits>&>(_Source);
  } else /* constexpr   _InBits < _OutBits */ {
    std::bitset<_OutBits> _Destination;
    std::copy(reinterpret_cast<const char*>(&_Source), 
              reinterpret_cast<const char*>(&_Source) + sizeof(_Source), 
              reinterpret_cast<char*>(&_Destination));
    return _Destination;
  }
}
  }// namespace std::bitset extension



#include <math.h>
#include <assert.h>
namespace calculation 
{

#define floating_operation_special_case(PROC_ZERO, PROC_INF, PROC_NAN) \
if ( (floating & ~sign_mask) == zero_bitset ) {       \
	/** special case 'zero' */                         \
  PROC_ZERO;                                         \
}                                                    \
if ( (floating & inf_bitset) == inf_bitset ) {       \
	if ( (floating & significand_mask) != 0 ) {        \
		/** special case 'NaN' */                        \
		PROC_NAN;                                        \
	}                                                  \
                                                    \
	/** special case 'inf' */                       \
	PROC_INF;                              \
}                               \

template<typename Bitset> 
void floating_significand(
  const Bitset& floating, 
        Bitset& significand,
	const Bitset zero_bitset,
	const Bitset inf_bitset,
	const Bitset nan_bitset,
	const Bitset significand_mask,
	const Bitset sign_mask,
	const Bitset exponent_bias_bitset,
	const Bitset possign_bitset)
{
  floating_operation_special_case(
    significand = zero_bitset;
    return,
    significand = nan_bitset;
    return,
    significand = nan_bitset;
    return  
    );

  significand = floating & significand_mask;
  significand |= exponent_bias_bitset;
  significand |= possign_bitset;
}

template<typename Bitset>
void floating_exponent(
   const  Bitset& floating, 
  /*OUT*/ Bitset& exponent, 
  /*Reg1*/Bitset& exponent_part, 
  /*Reg2*/Bitset& sign_part,
  const Bitset zero_bitset,
  const Bitset inf_bitset,
  const Bitset nan_bitset,
  const Bitset significand_mask,
  const Bitset exponent_mask, 
  const Bitset sign_mask,
  const Bitset exponent_bias_bitset, 
  const size_t last_bit_offset,
  const size_t hidden_significant_offset,
  const size_t exponent_offset,
  const Bitset possign_bitset,
  const Bitset negsign_bitset)
{
  floating_operation_special_case(
    exponent = zero_bitset;
    return,
    exponent = inf_bitset;
    return,
    exponent = nan_bitset;
    return  
    );

  // get unnormalized-exponent and exponent-sign
  exponent = floating & exponent_mask;
  if ( exponent == exponent_bias_bitset ) {
    exponent = zero_bitset;
    return;
  } else if ( exponent > exponent_bias_bitset ) {
    exponent -= exponent_bias_bitset;
    sign_part = possign_bitset;
  } else   /* exponent < exponent_bias_bitset */{
    exponent = exponent_bias_bitset - exponent;
    sign_part = negsign_bitset;
  }
				
  // normalize select shift
  size_t smaller_shift = 0;
  while ( !exponent.test(last_bit_offset - smaller_shift) ) {
    ++smaller_shift;
  }
  smaller_shift = last_bit_offset - hidden_significant_offset - smaller_shift;
  // normalize significand
  exponent >>= smaller_shift;
  exponent &= significand_mask;
  // normalize exponent
  exponent_part = smaller_shift;
  exponent_part <<= exponent_offset;
  exponent_part += exponent_bias_bitset;

  // combine
  exponent |= exponent_part;
  exponent |= sign_part;
}
		
template<typename Bitset>
void floating_sign(
  const Bitset& floating, 
        Bitset& sign,
  const Bitset zero_bitset,
  const Bitset exponent_bias_bitset,
  const Bitset sign_mask)
{
  sign = zero_bitset;
  sign |= exponent_bias_bitset;
  sign |= (floating & sign_mask);
}


template<typename Bitset>
bool need_roundup_for_nearest_even(
  const Bitset& significand, 
  Bitset& sticky_mask, 
  size_t smaller_shift) 
{
/**
 * the rounding has two kinds, 
 *   one is round-up, 
 *   another is round-down.
 * 
 * so round-error has two kinds,
 *   one is abs(round_up_number - origin_number)
 *   another is abs(round_down_number - origin_number)
 * 
 * @nearest-even method purpose is
 *   min( abs(round_up_number - origin_number), abs(round_down_number - origin_number) )
 * 
 * @understanding structure
 *              1.1111111 11111111 11111111  :IEEE-32bit significand
 *      ...
 *      ... do some arthmetic
 *   = 1.0001011 00011111 00000111 10000100  :need smaller_shift 8 places
 * 
 *                                           sticky_bits
 *                                   back_bit   |
 *                                        |  |--+-...
 *              1.0001011 00011111 00000111 10000100  :shifted
 *                                          |
 *                                       round-bit
 * 
 * @analysis
 *   round-up-number   = 1.0001011 00011111 00001000
 *   round-down-number = 1.0001011 00011111 00000111
 * 
 *   abs(round-up-number - origin-number)
 *   = 1.0001011 00011111 00001000
 *    -1.0001011 00011111 00000111 10000100
 *   = 0.0000000 00000000 00000000 01111100
 * 
 *   abs(round-down-number - origin-number) 
 *   = 1.0001011 00011111 00000111 10000100 
 *    -1.0001011 00011111 00000111
 *   = 0.0000000 00000000 00000000 10000100
 *   
 *   so we round-up.
 *   we think back on this calculation, we can see that.
 * 
 *   if round-bit is '1', 
 *     round-up-error must<= epsilon/2
 *     round-down-error must>= epsilon/2
 *   then if sticky-bits not 'all(0)'
 *     round-up-error must< epsilon/2
 *     round-down-error must> epsilon/2
 *   so we round-up.
*/
  if ( smaller_shift != 0/* && smaller_shift - 1 < significand.size()*/ ) {
    bool round_bit = significand.test( smaller_shift - 1 );
    bool sticky_bit = false;
    if ( smaller_shift != 1 ) {
      sticky_mask = 1;
      sticky_mask <<= (smaller_shift - 1);
      sticky_mask -= 1;
      sticky_bit = (significand & sticky_mask) != 0;
    }
  
    if ( round_bit && (sticky_bit | significand.test(smaller_shift)) ) {
      return true;
    }
  }

  return false;
}

template<typename Bitset>
void add_floating(
  /*INOUT*/Bitset& this_floating,
  /*Reg1*/ Bitset& this_significant, 
  /*Reg2*/ Bitset& this_exponent,
  /*Reg3*/ Bitset& this_sign,
    const  Bitset& right_floating, 
  /*Reg4*/ Bitset& right_significant,
  /*Reg5*/ Bitset& right_exponent,
  /*Reg6*/ Bitset& right_sign,
  /*Reg7*/ Bitset& exponent_difference,
  /*Reg8*/ Bitset& sticky_mask,
  const Bitset zero_bitset,
  const Bitset inf_bitset,
  const Bitset nan_bitset, 
  const Bitset significand_mask,
  const Bitset hidden_significand_mask,
  const Bitset exponent_mask, 
  const Bitset sign_mask,
  const size_t exponent_offset,
  const size_t hidden_significant_offset,
  const size_t lastbit_offset) 
{
  /** 0. check special case 'zero' */
  if ( (this_floating & (~sign_mask)) == zero_bitset ) {
    this_floating = right_floating;
    return;
  }
  if ( (right_floating & (~sign_mask)) == zero_bitset ) {
    return;
  }

  /** 0. check special case 'inf' and 'nan' */
  if ( (this_floating & inf_bitset) == inf_bitset ) {
    this_floating = nan_bitset;
    return;
  }
  if ( (right_floating & inf_bitset) == inf_bitset ) {
    this_floating = nan_bitset;
    return;
  }

  /** 1. decode into significand-bitset, exponent-bitset, sign-bitset */
  this_significant = (this_floating&significand_mask)|hidden_significand_mask;
  this_exponent = this_floating&exponent_mask;
  this_sign = this_floating&sign_mask;
  right_significant = (right_floating&significand_mask)|hidden_significand_mask;
  right_exponent = right_floating&exponent_mask;
  right_sign = right_floating&sign_mask;
  
  /** 2. extend more accuracy */
  assert( hidden_significant_offset <= lastbit_offset );
  size_t shift = 0;
  if ( hidden_significant_offset < lastbit_offset ) {
    /**
     *    11111111
     * +  11111111
     * = 111111110
     * 
     * so we only reserved '1' bit
    */
    const size_t reserved_bit = 1;
    shift = lastbit_offset - hidden_significant_offset;
    shift -= reserved_bit;
    this_significant <<= shift;
    right_significant <<= shift;
  }

  /** 3. sync exponent-bitset, to greater */
  if ( this_exponent < right_exponent ) {

    exponent_difference = right_exponent - this_exponent;
    size_t smaller_shift = static_cast<size_t>( (exponent_difference >> exponent_offset).to_ullong() );
    
    if ( smaller_shift - 1 <= lastbit_offset && need_roundup_for_nearest_even(this_significant, sticky_mask, smaller_shift) ){
      this_significant >>= smaller_shift;
      this_significant += 1;
    } else {
      this_significant >>= smaller_shift;
    }
    
    // assert( inf_bitset - exponent_difference < this_exponent );
    this_exponent += exponent_difference;
  
  } 
  else if ( this_exponent != right_exponent ) {

    exponent_difference = this_exponent - right_exponent;
    size_t smaller_shift = static_cast<size_t>( (exponent_difference >> exponent_offset).to_ullong() );

    if ( smaller_shift - 1 <= lastbit_offset && need_roundup_for_nearest_even(right_significant, sticky_mask, smaller_shift) ){
      right_significant >>= smaller_shift;
      right_significant += 1;
    } else {
      right_significant >>= smaller_shift;
    }

    // assert( inf_bitset - exponent_difference < right_exponent );
    right_exponent += exponent_difference;

  }

  /** 4. add significand-bitset, update sign-bits */
  if ( this_sign == right_sign ) {
    this_significant += right_significant;
  } 
  else {
    if (this_significant == right_significant) {
      this_floating = zero_bitset;
      return;
    } else if (this_significant < right_significant) {
      this_significant = right_significant - this_significant;
      this_sign = right_sign;
    } else  /* this_significand > right_significand */{
      this_significant -= right_significant;
      //this_sign = this_sign;
    }
  }

  /** 5. normalize significand-bits and exponent-bits */
  if ( (this_significant&(~significand_mask)) != 0 ) {

    // normalize for select 'shift'
    size_t smaller_shift = 0;
    while ( !this_significant.test(lastbit_offset - smaller_shift) ) {
      ++smaller_shift;
    }
    smaller_shift = lastbit_offset - hidden_significant_offset - smaller_shift;

    // normalize for 'significand' with 'round-up'
    if ( need_roundup_for_nearest_even(this_significant, sticky_mask, smaller_shift) ){
      this_significant >>= smaller_shift;
      this_significant += 1;
      if ( this_significant.test(hidden_significant_offset + 1) ) {
        this_significant >>= 1;
        smaller_shift += 1;
        /**
         *     1111 1111
         * +   0000 0001
         * = 1 0000 0000
         *            |
         *         this is test-bit
         *            |
         *         it must be zero, if carry
        */
      }
    } else {
      this_significant >>= smaller_shift;
    }

    // normalize for 'exponent'
    if ( smaller_shift > shift ) {
      exponent_difference = smaller_shift - shift;
      exponent_difference <<= exponent_offset;
      if ( inf_bitset - exponent_difference <= this_exponent ) {
        /** special case 'infinite large' */
        this_floating = inf_bitset|this_sign;
        return;
      }
      this_exponent += exponent_difference;
    } else {
      exponent_difference = shift - smaller_shift;
      exponent_difference <<= exponent_offset;
      if ( this_exponent < exponent_difference ) {
        /** special case 'infinite small' */
        this_floating = zero_bitset;
        return;
      }
      this_exponent -= exponent_difference;
    }
    
  } 
  else {

    // normalize for select 'shift'
    size_t larger_shift = 0;
    while ( ! this_significant.test(hidden_significant_offset - larger_shift) ) {
      ++larger_shift;
    }

    // normalize for 'significand'
    this_significant <<= larger_shift;

    // normalize for 'exponent'
    exponent_difference = shift + larger_shift;
    exponent_difference <<= exponent_offset;
    if ( this_exponent < exponent_difference ) {
      /** special case 'infinite small' */
      this_floating = zero_bitset;
      return;
    }
    this_exponent -= exponent_difference;

  }

  /** 6. encode ... */
  this_floating = this_sign|this_exponent|(this_significant&significand_mask);
}

template<typename Bitset, typename BigBitset>
void mul_floating(
  /*INOUT*/Bitset&    this_floating,
  /*Reg1*/ BigBitset& this_significand,
  /*Reg2*/ Bitset&    this_exponent,
    const  Bitset&    right_floating,
  /*Reg3*/ BigBitset& right_significand,
  /*Reg4*/ Bitset&    right_exponent,
  /*Reg5*/ Bitset&    exponent_difference,
  /*Reg6*/ BigBitset& sticky_mask,
  const Bitset zero_bitset,
  const Bitset inf_bitset,
  const Bitset nan_bitset,
  const Bitset significand_mask,
  const Bitset hidden_significand_mask,
  const Bitset exponent_mask,
  const Bitset sign_mask,
  const Bitset exponent_bias_bitset,
  const size_t significand_bits,
  const size_t bits)
{
  /** 0. check special case 'zero' */
  if ( (this_floating & (~sign_mask)) == zero_bitset || (right_floating & (~sign_mask)) == zero_bitset ) {
    this_floating = zero_bitset;
    return;
  }

  /** 0. check special case 'inf' and 'nan' */
  if ( (this_floating&inf_bitset) == inf_bitset || (right_floating&inf_bitset) == inf_bitset ) {
    if( (this_floating&significand_mask) != 0 || (right_floating&significand_mask) != 0 ) {
      this_floating = nan_bitset;
      return;
    }

    this_floating = inf_bitset | ((this_floating&sign_mask)^(right_floating&sign_mask));
    return;
  }

  /** 1. decode ... */
  this_significand = 0;
  reinterpret_cast<Bitset&>(this_significand) = (this_floating&significand_mask)|hidden_significand_mask;
  this_exponent = this_floating&exponent_mask;
  right_significand = 0;
  reinterpret_cast<Bitset&>(right_significand) = (right_floating&significand_mask)|hidden_significand_mask;
  right_exponent = right_floating&exponent_mask;

  /** 2. significand multiplication 
   * this_significand<<(significand_bits*2) = (this_significand<<significand_bits * right_significand<<significand_bits)
   * this_significand<<significand_bits = (this_significand<<significand_bits * right_significand<<significand_bits) >> significand_bits
   * this_significand<<'shift'          = (this_significand<<significand_bits * right_significand<<significand_bits) >> significand_bits
  */
  this_significand *= right_significand;
  size_t shift = significand_bits;
  
  /** 3. this_exponent = this_exponent + right_exponent */
  if ( right_exponent > exponent_bias_bitset ) {  
    this_exponent += (right_exponent - exponent_bias_bitset);
  } else if ( right_exponent < exponent_bias_bitset ) {
    this_exponent -= (exponent_bias_bitset - right_exponent);
  }
  
  /** 4. normalize, assert( this_significand&(~significand_mask) != 0 ); */
  {
    // normalize for select 'shift'
    const size_t lastbit_offset = bits * 2 - 1;
    const size_t hidden_significand_offset = significand_bits;
    size_t smaller_shift = 0;
    while ( !this_significand.test(lastbit_offset - smaller_shift) ) {
      ++smaller_shift;
    }
    smaller_shift = lastbit_offset - hidden_significand_offset - smaller_shift;

    // normalize for 'significand' with 'round-up'
    if ( need_roundup_for_nearest_even(this_significand, sticky_mask, smaller_shift) ) {
      this_significand >>= smaller_shift;
      this_significand += 1;
      if ( this_significand.test(hidden_significand_offset + 1) ) {
        this_significand >>= 1;
        smaller_shift += 1;
        // assert( round-bit == 0 );
      }
    } else {
      this_significand >>= smaller_shift;
    }

    // normalize for 'exponent'
    size_t exponent_offset = significand_bits;
    if ( smaller_shift > shift ) {
      exponent_difference = (smaller_shift - shift);
      exponent_difference <<= exponent_offset;
      if ( inf_bitset - exponent_difference <= this_exponent ) {
        /** special case 'infinite large' */
        this_floating = inf_bitset|((this_floating&sign_mask)^(right_floating&sign_mask));
        return;
      }
      this_exponent += exponent_difference;
    } else {
      exponent_difference = (shift - smaller_shift);
      exponent_difference <<= exponent_offset;
      if ( this_exponent < exponent_difference ) {
        /** special case 'infinite small' */
        this_floating = zero_bitset;
        return;
      }
      this_exponent -= exponent_difference;
    }
  }

  /** 4. encode ... */
  this_floating = 
    (reinterpret_cast<const Bitset&>(this_significand)&significand_mask)
    | this_exponent
    | ( (this_floating&sign_mask)^(right_floating&sign_mask) );
}

template<typename Bitset>
void div_floating(
  /*INOUT*/Bitset& this_floating,
  /*Reg1*/ Bitset& dividend,
  /*Reg2*/ Bitset& this_significand,
  /*Reg3*/ Bitset& this_exponent,
    const  Bitset& right_floating, 
  /*Reg4*/ Bitset& divisor,
  /*Reg5*/ Bitset& right_exponent,
  /*Reg6*/ Bitset& exponent_difference,
  /*Reg7*/ Bitset& divbit,
  /*Reg8*/ Bitset& sticky_mask,
  const Bitset zero_bitset,
  const Bitset inf_bitset,
  const Bitset nan_bitset,
  const Bitset significand_mask,
  const Bitset hidden_significand_mask,
  const Bitset exponent_mask,
  const Bitset sign_mask,
  const Bitset exponent_bias_bitset,
  const size_t exponent_offset,
  const size_t hidden_significant_offset,
  const size_t lastbit_offset)
{
  /** 0. check special case 'zero' */
  if ( (this_floating & (~sign_mask)) == zero_bitset && (right_floating & (~sign_mask)) == zero_bitset ) {
    this_floating = nan_bitset;
    return;
  } else if ( (this_floating & (~sign_mask)) == zero_bitset ) {
    // not change
    return;
  } else if ( (right_floating & (~sign_mask)) == zero_bitset ) {
    this_floating = inf_bitset | ((this_floating&sign_mask)^(right_floating&sign_mask));
    return;
  }

  /** 0. check special case 'inf' and 'nan' */
  if ( (this_floating & inf_bitset) == inf_bitset || (right_floating & inf_bitset) == inf_bitset ) {
    this_floating = nan_bitset;
    return;
  }

  /** 1. decode into divide-operator, exponent-bitset */
  dividend = (this_floating&significand_mask)|hidden_significand_mask;
  this_exponent = this_floating&exponent_mask;
  divisor = (right_floating&significand_mask)|hidden_significand_mask;
  right_exponent = right_floating&exponent_mask;
  
  /** 2. extend more accuracy, significand-largest-bit
   * [1<<(significand_bits-1), 1<<significand_bits] 
   * to 
   * [1<<(significand_bits+1), 1<<(significand_bits+2)] */
  assert( hidden_significant_offset <= lastbit_offset );
  size_t shift = 0;
  if ( hidden_significant_offset < lastbit_offset ) {
    /**
     *    11111111
     * /  10000000
     * ------------
     * =  11111111
     * 
     * so we only reserved '0' + '1'(division) bit
    */
    const size_t reserved_bit = 1;
    shift = lastbit_offset - hidden_significant_offset;
    shift -= reserved_bit;
    dividend <<= shift;
    divisor <<= shift;
  }

  /** 3. significand division */
  this_significand = 0;
  size_t offset = hidden_significant_offset + shift;// contains hidden-significand
  do {
    if ( dividend >= divisor ) {
      divbit = 1;  divbit <<= offset;
      this_significand |= divbit;
      dividend -= divisor;
    }

    dividend <<= 1;

  } while ( --offset != size_t(-1) );

  /** 4. this_exponent = this_exponent - right_exponent */
  if (right_exponent > exponent_bias_bitset) { 
    this_exponent -= (right_exponent - exponent_bias_bitset);
  } else if (right_exponent < exponent_bias_bitset) { 
    this_exponent += (exponent_bias_bitset - right_exponent);
  }

  /** 5. normalize, assert( this_significand&(~significand_mask) != 0 ); see(step 2.) */
  {
    // normalize for select 'shift'
    size_t smaller_shift = 0;
    while ( !this_significand.test(lastbit_offset - smaller_shift) ) {
      ++smaller_shift;
    }
    smaller_shift = lastbit_offset - hidden_significant_offset - smaller_shift;

    // normalize for 'significand' with 'round-up'
    if ( need_roundup_for_nearest_even(this_significand, sticky_mask, smaller_shift) ) {
      this_significand >>= smaller_shift;
      this_significand += 1;
      if ( this_significand.test(hidden_significant_offset + 1) ) {
        this_significand >>= 1;
        smaller_shift += 1;
        // assert( round-bit == 0 );
      }
    } else {
      this_significand >>= smaller_shift;
    }

    // normalize for 'exponent'
    assert( smaller_shift <= shift );
    exponent_difference = shift - smaller_shift;
    exponent_difference <<= exponent_offset;
    if ( inf_bitset - exponent_difference <= this_exponent ) {
      /** special case 'infinite small' */
      this_floating = zero_bitset;
      return;
    }
    this_exponent -= exponent_difference;
  }

  /** 6. encode ... */
  this_floating = 
    (this_significand&significand_mask)
    | this_exponent
    | ( (this_floating&sign_mask)^(right_floating&sign_mask) );
}



/** 
 * TEMPLATE infinite precision floating_point
 *   digit: (1 + 0.Significand) * 2^Exponent * (-1)^Sign, 
 *   __s + __e + 1 = Bits.
 * @tparam __s SignificandBits
 * @tparam __e ExponentBits
*/
template<size_t __s, size_t __e, bool _IsOpt = false>
class floatX {
public:
  static_assert( __s != 0, "assert(significand_bits != 0)" );
  static_assert( __e != 0, "assert(exponent_bits != 0)" );
    
  static constexpr size_t significand_bits = __s;
  static constexpr size_t exponent_bits = __e;
  static constexpr size_t sign_bits = 1;
  static constexpr size_t bits = significand_bits + exponent_bits + sign_bits;
    
  static constexpr size_t significand_offset = 0;
  static constexpr size_t exponent_offset = significand_offset + significand_bits;
  static constexpr size_t sign_offset = exponent_offset + exponent_bits;
		
  static std::bitset<bits> significand_mask() {
    static auto _Mask = 
      (std::bitset<bits>(1) << significand_bits) - 1 
        /* << 0 */;
    return _Mask;
  }

  static std::bitset<bits> exponent_mask() {
    static auto _Mask = 
      ( (std::bitset<bits>(1) << exponent_bits) - 1 ) 
        << exponent_offset;
    return _Mask;
  }

  static std::bitset<bits> sign_mask() {
    static auto _Mask = 
      std::bitset<bits>(1) << sign_offset;
    return _Mask;
  }

  static std::bitset<bits> hidden_significand_mask() {
    static auto _Mask = 
      std::bitset<bits>(1) << significand_bits;
    return _Mask;
  }

  static std::bitset<bits> exponent_bias_bitset() {
    static auto _Exp2_zero = 
      ( (std::bitset<bits>(1) << (exponent_bits - 1)) - 1 )
        << exponent_offset;
    return _Exp2_zero;
  }

  static std::bitset<bits> zero_bitset() {
    return std::bitset<bits>(0);
  }

  static std::bitset<bits> inf_bitset() {
    // inf = full 1 in exponent_bitset
    return exponent_mask();
  }

  static std::bitset<bits> quiet_NaN_bitset() {
    static auto _Quiet_nan =
      exponent_mask() 
        | (std::bitset<bits>(1)<<(significand_bits - 1)) ;
    return _Quiet_nan;
  }

  static std::bitset<bits> signaling_NaN_bitset() {
    static auto _Signaling_nan = 
      exponent_mask() 
        | (std::bitset<bits>(1)<<(significand_bits - 1)) 
          | std::bitset<bits>(1);
    return _Signaling_nan;
  }

  static std::bitset<bits> epsilon_bitset() {
    // epsilon = exp2(0 - mantissa_bits) 
    static auto _Epsilon = 
      exponent_bias_bitset() - (std::bitset<bits>(significand_bits) << exponent_offset);
    return _Epsilon;
  }

  static std::bitset<bits> possign_bitset() {
    return zero_bitset();
  }

  static std::bitset<bits> negsign_bitset() {
    return sign_mask();
  }

public:
  std::bitset<bits> _Mybitset;

  const std::bitset<bits>& bitset() const {
    return _Mybitset;
  }

  std::bitset<bits>& bitset() {
    return _Mybitset;
  }

  explicit floatX(std::bitset<bits> _bitset) : _Mybitset(_bitset) {}

  static bool iszero(const floatX& x) {
    return (x.bitset() & (exponent_mask() | significand_mask())) == zero_bitset();
  }

  static bool isinf(const floatX& x) {
    return (x.bitset() & exponent_mask()) == exponent_mask();
  }

  static bool isnan(const floatX& x) {
    return ((x.bitset() & exponent_mask()) == exponent_mask()) && ((x.bitset() & significand_mask()) != 0);
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
      _Mybitset = zero_bitset();
      return;
    }

    using this_float = floatX<__s, __e>;
    using other_float = floatX<m2, e2>;
    constexpr 
    size_t other_bits = other_float::bits;
    auto other_bitset = other.bitset();

    // _My_mantissa = shifted(other_mantissa)
    std::bitset<other_bits> other_mantissa = other_bitset & other_float::significand_mask();
    std::bitset<bits> _My_mantissa;
    if constexpr ( significand_bits < other_float::significand_bits ) {
      std::bitset<other_bits> shifted_mantissa = other_mantissa >> (other_float::significand_bits - significand_bits);
      std::bitset<bits> casted_shifted_mantissa = std::bitset_cast<bits>(shifted_mantissa);
      _My_mantissa = casted_shifted_mantissa;
    } else if constexpr ( other_float::significand_bits < significand_bits ) {
      std::bitset<bits> casted_mantissa = std::bitset_cast<bits>( other_mantissa );
      std::bitset<bits> shifted_casted_mantissa = casted_mantissa << (significand_bits - other_float::significand_bits);
      _My_mantissa = shifted_casted_mantissa;
    } else {
      _My_mantissa = std::bitset_cast<bits>( other_mantissa );
    }
			
    // decompose other_exponent
    std::bitset<other_bits> other_abs_exponent = (other_bitset & other_float::exponent_mask());
    bool other_exponent_sign;
    if ( other_abs_exponent < other_float::exponent_bias_bitset() ) {
      other_exponent_sign = 1;//negative
      other_abs_exponent = other_float::exponent_bias_bitset() - other_abs_exponent;
    } else {
      other_exponent_sign = 0;//positive
      other_abs_exponent -= other_float::exponent_bias_bitset();
    }
			
    // other_exponent is infinity ?, abandon an exponent-upper [-125,128] to [-125,127]
    std::bitset<bits> _My_exponent;
    std::bitset<other_bits> _My_abs_infinite;
    if constexpr ( this_float::exponent_offset < other_float::exponent_offset ) {
      std::bitset<other_bits> casted_infinite = std::bitset_cast<other_bits>( this_float::inf_bitset() - this_float::exponent_bias_bitset() - 1 );
      std::bitset<other_bits> shifted_casted_infinite = casted_infinite << ( other_float::exponent_offset - this_float::exponent_offset );
      _My_abs_infinite = shifted_casted_infinite;
    } else if constexpr ( this_float::exponent_offset > other_float::exponent_offset ) {
      std::bitset<bits> shifted_infinite = (inf_bitset() - exponent_bias_bitset() - 1) >> (exponent_offset - other_float::exponent_offset);
      std::bitset<other_bits> casted_shifted_infinite = std::bitset_cast<other_bits>(shifted_infinite);
      _My_abs_infinite = casted_shifted_infinite;
    } else /*if ( other_float::exponent_offset == this_float::exponent_offset )*/ {
      _My_abs_infinite = std::bitset_cast<other_bits>( inf_bitset() - exponent_bias_bitset() - 1 );
    }

    if (_My_abs_infinite < other_abs_exponent) {
      _My_exponent = inf_bitset();
    } else {
      // _My_exponent = shifted(other_exponent)
      if constexpr ( exponent_offset < other_float::exponent_offset ) {
        std::bitset<other_bits> shifted_exponent = other_abs_exponent >> (other_float::exponent_offset - exponent_offset);
        std::bitset<bits> casted_shifted_exponent = std::bitset_cast<bits>(shifted_exponent);
        _My_exponent = casted_shifted_exponent;
      } else if constexpr ( other_float::exponent_offset < exponent_offset ) {
        std::bitset<bits> casted_exponent = std::bitset_cast<bits>( other_abs_exponent );
        std::bitset<bits> shifted_casted_exponent = casted_exponent << (exponent_offset - other_float::exponent_offset);
        _My_exponent = shifted_casted_exponent;
      } else {
        _My_exponent = std::bitset_cast<bits>( other_abs_exponent );
      }

      if (other_exponent_sign) {
        _My_exponent = exponent_bias_bitset() - _My_exponent;
      } else {
        _My_exponent += exponent_bias_bitset();
      }
    }
			
    bool is_negative = (other_bitset & other_float::sign_mask()) 
      == (decltype(other_bitset)(1) << other_float::sign_offset);
    auto _My_sign = std::bitset<bits>(is_negative) << sign_offset;

    _Mybitset = _My_sign | _My_exponent | _My_mantissa;
  }
		
  floatX(float ohter) : floatX(reinterpret_cast<const floatX<23,8>&>(ohter)) {}
		
  floatX(double ohter) : floatX(reinterpret_cast<const floatX<52,11>&>(ohter)) {}

  floatX(long double other) : floatX(static_cast<double>(other)) {}

  floatX(unsigned int other) {
    if ( other == 0U ) {
      _Mybitset = zero_bitset();
      return;
    }

    // 000000000.00000000000010100100010 * pow(2,significand_bits)
    auto this_exponent = exponent_bias_bitset() + (std::bitset<bits>(significand_bits) << exponent_offset);
    auto this_significant = std::bitset<bits>(other);

    // normalize significand-bits and exponent-bits
    if ( (this_significant&(~significand_mask())) != 0 ) {

      // normalize for select 'shift'
      const size_t hidden_significand_offset = significand_bits;
      const size_t lastbit_offset = bits - 1;
      size_t smaller_shift = 0;
      while ( !this_significant.test(lastbit_offset - smaller_shift) ) {
        ++smaller_shift;
      }
      smaller_shift = lastbit_offset - hidden_significand_offset - smaller_shift;

      // normalize for 'significand'
      this_significant >>= smaller_shift;
        
      // normalize for 'exponent'
      this_exponent += std::bitset<bits>(smaller_shift) << exponent_offset;

    } 
    else {
        
      // normalize for select 'shift'
      const size_t hidden_significand_offset = significand_bits;
      size_t larger_shift = 0;
      while ( ! this_significant.test(hidden_significand_offset - larger_shift) ) {
        ++larger_shift;
      }

      // normalize for 'significand'
      this_significant <<= larger_shift;

      // normalize for 'exponent'
      this_exponent -= std::bitset<bits>(larger_shift) << exponent_offset;

    }

    _Mybitset = (this_significant&significand_mask())|this_exponent;
  }

  floatX(unsigned long long other) {
    if (other == 0ULL) {
      _Mybitset = zero_bitset();
      return;
    }

    // 1expexpexp.mantissamantissa * pow(2,significand_bits)
    auto this_exponent = exponent_bias_bitset() + (std::bitset<bits>(significand_bits) << exponent_offset);
    auto this_significant = std::bitset<bits>(other);

    //  normalize significand-bits and exponent-bits
    //auto hidden_significand_mask = ~significand_mask();
    auto exponent_one = std::bitset<bits>(1) << exponent_offset;
    if ( this_significant < hidden_significand_mask() ) {
      while ( (this_significant & hidden_significand_mask()) != hidden_significand_mask() ) {
        this_significant <<= 1;
        this_exponent -= exponent_one;
      }
    } else {
      while ( (this_significant & hidden_significand_mask()) != hidden_significand_mask() ) {
      this_significant >>= 1;
      this_exponent += exponent_one;
      }
    }
      
    _Mybitset = this_exponent | (this_significant & significand_mask());
  }

  floatX(bool other) {
    if( other ){
      _Mybitset = exponent_bias_bitset()/* | std::bitset<bits>(0)*/;
    } else {
      _Mybitset = zero_bitset();
    }
  }

  floatX(int other) : floatX(other < 0 ? static_cast<unsigned int>(-other) : static_cast<unsigned int>(other)) {
    _Mybitset |= std::bitset<bits>(other < 0) << sign_offset;
  }

  floatX(long long other) : floatX(other < 0 ? static_cast<unsigned long long>(-other) : static_cast<unsigned long long>(other)) {
    _Mybitset |= std::bitset<bits>(other < 0) << sign_offset;
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
    const auto zero_exponent = exponent_bias_bitset();

    // only fraction
    auto this_exponent = (_Mybitset & exponent_mask());
    if ( this_exponent < zero_exponent ) {
      return static_cast<unsigned int>(0);
    }

    // compute saved_bits, @floor(float_)
    auto exp2_mantissa_bits = std::bitset<bits>(significand_bits) << exponent_offset;
    auto truncated_exponent = this_exponent - exp2_mantissa_bits;
    auto trunc_exponent = zero_exponent - truncated_exponent;
    size_t trunc_bits = (trunc_exponent >> exponent_offset).to_ulong();
    size_t saved_bits = (significand_bits+1) - trunc_bits;
    if (saved_bits > sizeof(unsigned int) * 8) {
      throw std::overflow_error("floatX<...>::operator unsigned int() const");
    }

    // smaller-shift to correct integer-bits
    auto this_significant = _Mybitset & significand_mask() | hidden_significand_mask();
    this_significant >>= ((significand_bits+1) - saved_bits);

    auto dest = std::bitset_cast<sizeof(unsigned int)*8>(this_significant);
    return reinterpret_cast<const unsigned int&>(dest);
  }

  explicit operator unsigned long long() const {
    const auto zero_exponent = exponent_bias_bitset();

    // only fraction
    auto this_exponent = (_Mybitset & exponent_mask());
    if ( this_exponent < zero_exponent ) {
      return static_cast<unsigned long long>(0);
    }

    // compute saved_bits, @floor(float_)
    auto exp2_mantissa_bits = std::bitset<bits>(significand_bits) << exponent_offset;
    auto truncated_exponent = this_exponent - exp2_mantissa_bits;
    auto trunc_exponent = zero_exponent - truncated_exponent;
    size_t trunc_bits = (trunc_exponent >> exponent_offset).to_ulong();
    size_t saved_bits = (significand_bits+1) - trunc_bits;
    if (saved_bits > sizeof(unsigned long long) * 8) {
      throw std::overflow_error("floatX<...>::operator unsigned int() const");
    }

    // smaller-shift to correct integer-bits
    auto this_significant = _Mybitset & significand_mask() | hidden_significand_mask();
    this_significant >>= ((significand_bits+1) - saved_bits);
       
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

    // check sign
    const auto this_sign = _Mybitset & sign_mask();
    const auto right_sign = right.bitset() & sign_mask();
    if (this_sign != right_sign) {
      if (this_sign > right_sign) {
        return true;
      } else {
        return false;
      }
    }

    // check exponent|significand
    if (this_sign == zero_bitset()) {
      // positive
      if( (_Mybitset & (~sign_mask())) < (right.bitset() & (~sign_mask())) ) {
        return true;
      } else {
        return false;
      }

    } else {
      // negative
      if( (_Mybitset & (~sign_mask())) < (right.bitset() & (~sign_mask())) ) {
        return false;
      } else {
        return true;
      }

    }
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

  floatX& operator+=(const floatX& right) {
    std::bitset<bits> this_significant;
    std::bitset<bits> this_exponent;
    std::bitset<bits> this_sign;
    std::bitset<bits> right_significant;
    std::bitset<bits> right_exponent;
    std::bitset<bits> right_sign;
    std::bitset<bits> exponent_difference;
    std::bitset<bits> sticky_mask;
    add_floating(_Mybitset, this_significant, this_exponent, this_sign,
      right.bitset(), right_significant, right_exponent, right_sign,
      exponent_difference, sticky_mask,
      zero_bitset(),
      inf_bitset(),
      signaling_NaN_bitset(),
      significand_mask(),
      hidden_significand_mask(),
      exponent_mask(),
      sign_mask(),
      exponent_offset,
      significand_bits,
      bits - 1
      );
    return *this;
  }

  floatX& operator-=(const floatX& right) {
    return *this += (-right);
  }

  floatX& operator*=(const floatX& right) {
    std::bitset<bits*2> this_significand;
    std::bitset<bits> this_exponent;
    std::bitset<bits*2> right_significand;
    std::bitset<bits> right_exponent;
    std::bitset<bits> exponent_difference;
    std::bitset<bits*2> roundup;
    mul_floating(_Mybitset, this_significand, this_exponent,
      right.bitset(), right_significand, right_exponent,
      exponent_difference, roundup,
      zero_bitset(),
      inf_bitset(),
      signaling_NaN_bitset(),
      significand_mask(),
      hidden_significand_mask(),
      exponent_mask(),
      sign_mask(),
      exponent_bias_bitset(),
      significand_bits,
      bits
      );
    return *this;
  }

  floatX& operator/=(const floatX& right) {
    std::bitset<bits> dividend;
    std::bitset<bits> this_significand;
    std::bitset<bits> this_exponent;
    std::bitset<bits> divisor;
    std::bitset<bits> right_exponent;
    std::bitset<bits> exponent_difference;
    std::bitset<bits> divbit;
    std::bitset<bits> roundup;
    div_floating(_Mybitset, dividend, this_significand, this_exponent,
      right.bitset(), divisor, right_exponent,
      exponent_difference, divbit, roundup,
      zero_bitset(),
      inf_bitset(),
      signaling_NaN_bitset(),
      significand_mask(),
      hidden_significand_mask(),
      exponent_mask(),
      sign_mask(),
      exponent_bias_bitset(),
      exponent_offset,
      significand_bits,
      bits - 1
      );
    return *this;
  }

  floatX& operator%=(const floatX& right) {
    abort();
  }
    
  floatX operator+(const floatX& right) const {
    return floatX(*this) += right;
  }

  floatX operator-(const floatX& right) const {
    return floatX(*this) -= right;
  }

  floatX operator*(const floatX& right) const {
    return floatX(*this) *= right;
  }

  floatX operator/(const floatX& right) const {
    return floatX(*this) /= right;
  }

  floatX operator%(const floatX& right) const {
    abort();
  }

  floatX significand() const {
    floatX the_significand;
    floating_significand(_Mybitset, the_significand.bitset(),
      zero_bitset(),
      inf_bitset(),
      signaling_NaN_bitset(),
      significand_mask(),
      sign_mask(),
      exponent_bias_bitset(),
      possign_bitset()
      );
    return the_significand;
  }

  floatX exponent() const {
    floatX the_exponent;
    floatX exponent_part;
    floatX sign_part;
    floating_exponent(_Mybitset, the_exponent.bitset(), exponent_part.bitset(), sign_part.bitset(),
      zero_bitset(),
      inf_bitset(),
      signaling_NaN_bitset(),
      significand_mask(),
      exponent_mask(),
      sign_mask(),
      exponent_bias_bitset(),
      bits - 1,
      significand_bits,
      exponent_offset,
      possign_bitset(),
      negsign_bitset()
    );
    return the_exponent;
  }

  floatX sign() const {
    floatX the_sign;
    floating_sign(_Mybitset, the_sign.bitset(),
      zero_bitset(),
      exponent_bias_bitset(),
      sign_mask()
    );
    return the_sign;
  }

  static floatX abs(const floatX& left) {
    const auto abs_mask = ~(sign_mask());
    return floatX{ left.bitset() & abs_mask };
  }

  static floatX trunc(const floatX& left) {
    const auto zero_exponent = exponent_bias_bitset();

    // only fraction, 1.010101010111... * exp2(0)
    auto left_exponent = left.bitset() & exponent_mask();
    if ( left_exponent < zero_exponent ) {
      return floatX{ zero_bitset() };
    }

    // check only integer, 1010101010111... * exp2(exponent - significand_bits), 
    auto exp2_mantissa_bits = std::bitset<bits>(significand_bits) << exponent_offset;
    auto truncated_exponent = left_exponent - exp2_mantissa_bits;
    if ( truncated_exponent >= zero_exponent ) {
      return left;
    }

    auto trunc_exponent = zero_exponent - truncated_exponent;
    size_t trunc_bits = (trunc_exponent >> exponent_offset).to_ulong();
    auto trunc_mask = ~( ( std::bitset<bits>(1) << trunc_bits ) - 1 );
    return floatX{ left.bitset() & trunc_mask };
  }

  static floatX frac(const floatX& left) {
    const auto zero_exponent = exponent_bias_bitset();

    // 1.010101010111... * exp2(0), check only fraction
    auto left_exponent = left.bitset() & exponent_mask();
    if ( left_exponent < zero_exponent ) {
      return left;
    }

    // 1010101010111... * exp2(exponent - significand_bits), check only integer
    auto exp2_mantissa_bits = std::bitset<bits>(significand_bits) << exponent_offset;
    auto truncated_exponent = left_exponent - exp2_mantissa_bits;
    if ( truncated_exponent >= zero_exponent ) {
      return floatX{ zero_bitset() };
    }

    auto trunc_exponent = zero_exponent - truncated_exponent;
    size_t trunc_bits = (trunc_exponent >> exponent_offset).to_ulong();
    auto fract_mask = (( std::bitset<bits>(1) << trunc_bits ) - 1);
			
    //  normalize significand-bits and exponent-bits, assert(larger_shift)
    auto left_significant = left.bitset() & fract_mask;
    auto exp2_one = std::bitset<bits>(1) << exponent_offset;
    assert ( left_significant < hidden_significand_mask() );
    while ( (left_significant & hidden_significand_mask()) != hidden_significand_mask() ) {
      left_significant <<= 1;
      left_exponent -= exp2_one;
    }
    return floatX{ (left.bitset() & sign_mask()) | left_exponent | (left_significant & significand_mask()) };
  }		
};

/** TEMPLATE SPECIAL, IEEE754 single-precision */
template<>
class floatX<23,8,true> {
  using _Mybase = floatX<23, 8, false>;
public:
  static constexpr size_t significand_bits = 23;
  static constexpr size_t exponent_bits = 8;
  static constexpr size_t sign_bits = 1;
  static constexpr size_t bits = 32;

  static constexpr size_t significand_offset = 0;
  static constexpr size_t exponent_offset = significand_offset + significand_bits;
  static constexpr size_t sign_offset = exponent_offset + exponent_bits;

  static constexpr std::bitset<bits> significand_mask() {
    return std::bitset<bits>(0b00000000011111111111111111111111);
  }
  static constexpr std::bitset<bits> exponent_mask() {
    return std::bitset<bits>(0b01111111100000000000000000000000);
  }
  static constexpr std::bitset<bits> sign_mask() {
    return std::bitset<bits>(0b10000000000000000000000000000000);
  }
  static constexpr std::bitset<bits> hidden_significand_mask() {
    return std::bitset<bits>(0b00000000100000000000000000000000);
  }
		
  static constexpr std::bitset<bits> exponent_bias_bitset() {
    return std::bitset<bits>(0b00111111100000000000000000000000);
  }
  static constexpr std::bitset<bits> zero_bitset() {
    return std::bitset<bits>(0b00000000000000000000000000000000);
  }
  static constexpr std::bitset<bits> inf_bitset() {
    return exponent_mask();
  }
  static constexpr std::bitset<bits> quiet_NaN_bitset() {
    return std::bitset<bits>(0b01111111110000000000000000000000);
  }
  static constexpr std::bitset<bits> signaling_NaN_bitset() {
    return std::bitset<bits>(0b01111111110000000000000000000001);
  }
  static constexpr std::bitset<bits> epsilon_bitset() {
    return std::bitset<bits>(0b00110100000000000000000000000000);
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

/**  TEMPLATE SPECIAL, IEEE754 double-precision */
template<>
class floatX<52,11,true> {
	using _Mybase = floatX<52, 11, false>;
public:
	static constexpr size_t significand_bits = 52;
	static constexpr size_t exponent_bits = 11;
	static constexpr size_t sign_bits = 1;
	static constexpr size_t bits = 64;

	static constexpr size_t significand_offset = 0;
	static constexpr size_t exponent_offset = significand_offset + significand_bits;
	static constexpr size_t sign_offset = exponent_offset + exponent_bits;

	static constexpr std::bitset<bits> significand_mask() {
		return std::bitset<bits>(0b0000000000001111111111111111111111111111111111111111111111111111);
	}
	static constexpr std::bitset<bits> exponent_mask() {
		return std::bitset<bits>(0b0111111111110000000000000000000000000000000000000000000000000000);
	}
	static constexpr std::bitset<bits> sign_mask() {
		return std::bitset<bits>(0b1000000000000000000000000000000000000000000000000000000000000000);
	}
	static constexpr std::bitset<bits> hidden_significand_mask() {
		return std::bitset<bits>(0b0000000000010000000000000000000000000000000000000000000000000000);
	}
		
	static constexpr std::bitset<bits> exponent_bias_bitset() {
		return std::bitset<bits>(0b0011111111110000000000000000000000000000000000000000000000000000);
	}
	static constexpr std::bitset<bits> zero_bitset() {
		return std::bitset<bits>(0b0000000000000000000000000000000000000000000000000000000000000000);
	}
	static constexpr std::bitset<bits> inf_bitset() {
		return exponent_mask();
	}
	static constexpr std::bitset<bits> quiet_NaN_bitset() {
    return std::bitset<bits>(0b0111111111111000000000000000000000000000000000000000000000000000);
	}
	static constexpr std::bitset<bits> signaling_NaN_bitset() {
    return std::bitset<bits>(0b0111111111111000000000000000000000000000000000000000000000000001);
	}
	static constexpr std::bitset<bits> epsilon_bitset() {
    return std::bitset<bits>(0b0011110010110000000000000000000000000000000000000000000000000000);
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

// next_work: { XXXX quadruple-precision-template spectialation }

#define __calculation_float_operator_with_literal(_OP_, _LITERAL_TYPE_)               \
template<size_t s, size_t e, bool opt> inline                                         \
floatX<s,e,opt> operator##_OP_##(const floatX<s,e,opt>& left, _LITERAL_TYPE_ right) { \
  return left _OP_ static_cast< floatX<s,e,opt> >(right);                             \
}

#define __calculation_float_operator_with_literal2(_OP_, _LITERAL_TYPE_)              \
template<size_t s, size_t e, bool opt> inline                                         \
floatX<s,e,opt> operator##_OP_##(_LITERAL_TYPE_ left, const floatX<s,e,opt>& right) { \
  return static_cast< floatX<s,e,opt> >(left) _OP_ right;                             \
}

#define __calculation_float_operator_with_literal_commutatibity(_OP_, _LITERAL_TYPE_) \
template<size_t s, size_t e, bool opt> inline                                         \
floatX<s,e,opt> operator##_OP_##(const floatX<s,e,opt>& left, _LITERAL_TYPE_ right) { \
  return left _OP_ static_cast< floatX<s,e,opt> >(right);                             \
}                                                                                     \
template<size_t s, size_t e, bool opt> inline                                         \
floatX<s,e,opt> operator##_OP_##(_LITERAL_TYPE_ left, const floatX<s,e,opt>& right) { \
  return static_cast< floatX<s,e,opt> >(left) _OP_ right;                             \
}

#define __calculation_float_lvalueoperator_with_literal(_OP_, _LITERAL_TYPE_)    \
template<size_t s, size_t e, bool opt> inline                                    \
floatX<s,e,opt>& operator##_OP_##(floatX<s,e,opt>& left, _LITERAL_TYPE_ right) { \
  return left _OP_ static_cast< floatX<s,e,opt> >(right);                        \
}

#define __calculation_float_comparison_with_literal_commutatibity(_OP_, _LITERAL_TYPE_) \
template<size_t s, size_t e, bool opt> inline                              \
bool operator##_OP_##(const floatX<s,e,opt>& left, _LITERAL_TYPE_ right) { \
  return left _OP_ static_cast< floatX<s,e,opt> >(right);                  \
}                                                                          \
template<size_t s, size_t e, bool opt> inline                              \
bool operator##_OP_##(_LITERAL_TYPE_ left, const floatX<s,e,opt>& right) { \
  return static_cast< floatX<s,e,opt> >(left) _OP_ right;                  \
} 

__calculation_float_lvalueoperator_with_literal(+=, float)
__calculation_float_lvalueoperator_with_literal(+=, double)
__calculation_float_lvalueoperator_with_literal(+=, bool)
__calculation_float_lvalueoperator_with_literal(+=, int)
__calculation_float_lvalueoperator_with_literal(+=, long long)
__calculation_float_lvalueoperator_with_literal(+=, unsigned int)
__calculation_float_lvalueoperator_with_literal(+=, unsigned long long)
__calculation_float_operator_with_literal_commutatibity(+, float)
__calculation_float_operator_with_literal_commutatibity(+, double)
__calculation_float_operator_with_literal_commutatibity(+, bool)
__calculation_float_operator_with_literal_commutatibity(+, int)
__calculation_float_operator_with_literal_commutatibity(+, long long)
__calculation_float_operator_with_literal_commutatibity(+, unsigned int)
__calculation_float_operator_with_literal_commutatibity(+, unsigned long long)

__calculation_float_lvalueoperator_with_literal(-=, float)
__calculation_float_lvalueoperator_with_literal(-=, double)
__calculation_float_lvalueoperator_with_literal(-=, bool)
__calculation_float_lvalueoperator_with_literal(-=, int)
__calculation_float_lvalueoperator_with_literal(-=, long long)
__calculation_float_lvalueoperator_with_literal(-=, unsigned int)
__calculation_float_lvalueoperator_with_literal(-=, unsigned long long)
__calculation_float_operator_with_literal(-, float)
__calculation_float_operator_with_literal2(-, float)
__calculation_float_operator_with_literal(-, double)
__calculation_float_operator_with_literal2(-, double)
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
__calculation_float_lvalueoperator_with_literal(*=, bool)
__calculation_float_lvalueoperator_with_literal(*=, int)
__calculation_float_lvalueoperator_with_literal(*=, long long)
__calculation_float_lvalueoperator_with_literal(*=, unsigned int)
__calculation_float_lvalueoperator_with_literal(*=, unsigned long long)
__calculation_float_operator_with_literal_commutatibity(*, float)
__calculation_float_operator_with_literal_commutatibity(*, double)
__calculation_float_operator_with_literal_commutatibity(*, bool)
__calculation_float_operator_with_literal_commutatibity(*, int)
__calculation_float_operator_with_literal_commutatibity(*, long long)
__calculation_float_operator_with_literal_commutatibity(*, unsigned int)
__calculation_float_operator_with_literal_commutatibity(*, unsigned long long)

__calculation_float_lvalueoperator_with_literal(/=, float)
__calculation_float_lvalueoperator_with_literal(/=, double)
__calculation_float_lvalueoperator_with_literal(/=, bool)
__calculation_float_lvalueoperator_with_literal(/=, int)
__calculation_float_lvalueoperator_with_literal(/=, long long)
__calculation_float_lvalueoperator_with_literal(/=, unsigned int)
__calculation_float_lvalueoperator_with_literal(/=, unsigned long long)
__calculation_float_operator_with_literal(/, float)
__calculation_float_operator_with_literal2(/, float)
__calculation_float_operator_with_literal(/, double)
__calculation_float_operator_with_literal2(/, double)
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
__calculation_float_lvalueoperator_with_literal(%=, bool)
__calculation_float_lvalueoperator_with_literal(%=, int)
__calculation_float_lvalueoperator_with_literal(%=, long long)
__calculation_float_lvalueoperator_with_literal(%=, unsigned int)
__calculation_float_lvalueoperator_with_literal(%=, unsigned long long)
__calculation_float_operator_with_literal_commutatibity(%, float)
__calculation_float_operator_with_literal_commutatibity(%, double)
__calculation_float_operator_with_literal_commutatibity(%, bool)
__calculation_float_operator_with_literal_commutatibity(%, int)
__calculation_float_operator_with_literal_commutatibity(%, long long)
__calculation_float_operator_with_literal_commutatibity(%, unsigned int)
__calculation_float_operator_with_literal_commutatibity(%, unsigned long long)

__calculation_float_comparison_with_literal_commutatibity(==, float)
__calculation_float_comparison_with_literal_commutatibity(==, double)
__calculation_float_comparison_with_literal_commutatibity(==, bool)
__calculation_float_comparison_with_literal_commutatibity(==, int)
__calculation_float_comparison_with_literal_commutatibity(==, long long)
__calculation_float_comparison_with_literal_commutatibity(==, unsigned int)
__calculation_float_comparison_with_literal_commutatibity(==, unsigned long long)

__calculation_float_comparison_with_literal_commutatibity(!=, float)
__calculation_float_comparison_with_literal_commutatibity(!=, double)
__calculation_float_comparison_with_literal_commutatibity(!=, bool)
__calculation_float_comparison_with_literal_commutatibity(!=, int)
__calculation_float_comparison_with_literal_commutatibity(!=, long long)
__calculation_float_comparison_with_literal_commutatibity(!=, unsigned int)
__calculation_float_comparison_with_literal_commutatibity(!=, unsigned long long)

__calculation_float_comparison_with_literal_commutatibity(<, float)
__calculation_float_comparison_with_literal_commutatibity(<, double)
__calculation_float_comparison_with_literal_commutatibity(<, bool)
__calculation_float_comparison_with_literal_commutatibity(<, int)
__calculation_float_comparison_with_literal_commutatibity(<, long long)
__calculation_float_comparison_with_literal_commutatibity(<, unsigned int)
__calculation_float_comparison_with_literal_commutatibity(<, unsigned long long)

__calculation_float_comparison_with_literal_commutatibity(>, float)
__calculation_float_comparison_with_literal_commutatibity(>, double)
__calculation_float_comparison_with_literal_commutatibity(>, bool)
__calculation_float_comparison_with_literal_commutatibity(>, int)
__calculation_float_comparison_with_literal_commutatibity(>, long long)
__calculation_float_comparison_with_literal_commutatibity(>, unsigned int)
__calculation_float_comparison_with_literal_commutatibity(>, unsigned long long)

__calculation_float_comparison_with_literal_commutatibity(<=, float)
__calculation_float_comparison_with_literal_commutatibity(<=, double)
__calculation_float_comparison_with_literal_commutatibity(<=, bool)
__calculation_float_comparison_with_literal_commutatibity(<=, int)
__calculation_float_comparison_with_literal_commutatibity(<=, long long)
__calculation_float_comparison_with_literal_commutatibity(<=, unsigned int)
__calculation_float_comparison_with_literal_commutatibity(<=, unsigned long long)

__calculation_float_comparison_with_literal_commutatibity(>=, float)
__calculation_float_comparison_with_literal_commutatibity(>=, double)
__calculation_float_comparison_with_literal_commutatibity(>=, bool)
__calculation_float_comparison_with_literal_commutatibity(>=, int)
__calculation_float_comparison_with_literal_commutatibity(>=, long long)
__calculation_float_comparison_with_literal_commutatibity(>=, unsigned int)
__calculation_float_comparison_with_literal_commutatibity(>=, unsigned long long)

#undef __calculation_float_operator_with_literal
#undef __calculation_float_operator_with_literal2
#undef __calculation_float_operator_with_literal_commutatibity
#undef __calculation_float_lvalueoperator_with_literal
#undef __calculation_float_comparison_with_literal_commutatibity

template<size_t m, size_t e, bool opt> inline
bool isinf(const floatX<m,e,opt>& x) {
  return (x.bitset() & floatX<m,e,opt>::exponent_mask()) == floatX<m,e,opt>::inf_bitset();
}  

template<size_t m, size_t e, bool opt> inline
bool isnan(const floatX<m,e,opt>& x) {
  return (x.bitset() & floatX<m,e,opt>::exponent_mask()) == floatX<m,e,opt>::inf_bitset()
      && (x.bitset() & floatX<m,e,opt>::significand_mask()) != floatX<m,e,opt>::zero_bitset();
}

template<size_t m, size_t e, bool opt> inline
floatX<m,e,opt> abs(const floatX<m,e,opt>& x) {
  auto abs_mask = ~( floatX<m,e,opt>::sign_mask() );
  return floatX<m,e,opt>{ x.bitset() & abs_mask };
}

template<size_t m, size_t e, bool opt> inline
floatX<m,e,opt> floor(const floatX<m,e,opt>& x) {
  if ( x >= 0 ) {
    return floatX<m,e,opt>::trunc(x);
  } else {
    return floatX<m,e,opt>::trunc(x) - 1;
  }
}

template<size_t m, size_t e, bool opt> inline
floatX<m,e,opt> ceil(const floatX<m,e,opt>& x) {
  return floor(x) + 1;
}

template<size_t m, size_t e, bool opt> inline
floatX<m,e,opt> round(const floatX<m,e,opt>& x) {
  return floor(x + 0.5);
}

template<size_t s, size_t e, bool opt> inline
floatX<s,e,opt> trunc(const floatX<s,e,opt>& x) {
  return floatX<s,e,opt>::trunc(x);
}

template<size_t s, size_t e, bool opt> inline
floatX<s,e,opt> frac(const floatX<s,e,opt>& x) {
  return floatX<s,e,opt>::frac(x);
}

template<size_t s, size_t e, bool opt> inline
floatX<s,e,opt> fmod(const floatX<s,e,opt>& x, const floatX<s,e,opt>& y) {
  return x - y * trunc(x/y);
}

/** IEEE754 single-precision */
using float32_t = floatX<23,8,true>;

/** IEEE754 double-precision */
using float64_t = floatX<52,11,true>;

/** (1 + 0.Significand) * 2^Exponent * (-1)^Sign,
 * @see GCC/quadmath/ALL  
*/
using float128_t = floatX<112,15>;

using float256_t = floatX<235,20>;

using float512_t = floatX<485,26>;

using float1024_t = floatX<992,31>;

inline bool isinf(float x) { return _CSTD isinf(x); }
inline bool isinf(double x) { return _CSTD isinf(x); }
inline bool isinf(float32_t x) { return _CSTD isinf(x); }
inline bool isinf(float64_t x) { return _CSTD isinf(x); }
inline bool isnan(float x) { return _CSTD isnan(x); }
inline bool isnan(double x) { return _CSTD isnan(x); }
inline bool isnan(float32_t x) { return _CSTD isnan(x); }
inline bool isnan(float64_t x) { return _CSTD isnan(x); }
using _CSTD abs;
inline float32_t abs(float32_t x) { return _CSTD fabsf(x); }
inline float64_t abs(float64_t x) { return _CSTD fabs(x); }
using _CSTD floor;
inline float32_t floor(float32_t x) { return _CSTD floorf(x); }
inline float64_t floor(float64_t x) { return _CSTD floor(x); }
using _CSTD ceil;
inline float32_t ceil(float32_t x) { return _CSTD ceilf(x); }
inline float64_t ceil(float64_t x) { return _CSTD ceil(x); }
using _CSTD round;
inline float32_t round(float32_t x) { return _CSTD roundf(x); }
inline float64_t round(float64_t x) { return _CSTD round(x); }
using _CSTD trunc;
inline float32_t trunc(float32_t x) { return _CSTD truncf(x); }
inline float64_t trunc(float64_t x) { return _CSTD trunc(x); }
inline float frac(float x) { return x - _CSTD truncf(x); }
inline double frac(double x) { return x - _CSTD trunc(x); }
inline float32_t frac(float32_t x) { return x - _CSTD truncf(x); }
inline float64_t frac(float64_t x) { return x - _CSTD trunc(x); }
using _CSTD fmod;
inline float32_t fmod(float32_t x, float32_t y) { return _CSTD fmodf(x, y); }
inline float64_t fmod(float64_t x, float64_t y) { return _CSTD fmod(x, y); }

}// namespace calculation


#include <cassert>
#include <deque>
#include <algorithm>
#include <xstring>
namespace calculation 
{

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
template<size_t m, size_t e, bool _Isbase>
std::string to_string(floatX<m,e,_Isbase> _Source) {
  using floatX_t = floatX<m, e, _Isbase>;

	if (isnan(_Source)) {
		return "nan";
	}
	if (isinf(_Source)) {
		return (_Source < floatX_t(0) ? "-inf" : "inf");
	}
	if (_Source == 0) {
		return "0.";
	} 
		
	decimal _Dest = decimal("0.");

	// _Dest.numbers = significand * 2^(significand_bits)
	decimal _Exp2_n = decimal("1.");
	auto _Mymantissa = _Source.bitset() & _Source.significand_mask();
	for (size_t i = 0; i != _Source.significand_bits; ++i, _Exp2_n.mulev2()) {
		if (_Mymantissa.test(i)) {
			_Dest += _Exp2_n;
		}
	}
	_Dest += _Exp2_n;// hidden-bit

	// correct exponent
	auto _Zero = std::bitset<floatX_t::bits>(0);
	auto _One = std::bitset<floatX_t::bits>(1) << _Source.exponent_offset;
	auto _Exponent = _Source.bitset() & _Source.exponent_mask();
	_Exponent -= std::bitset<floatX_t::bits>(_Source.significand_bits) << _Source.exponent_offset;
	if (_Exponent > _Source.exponent_bias_bitset()) {
		do {
			_Dest.mulev2();
			_Exponent -= _One;
		} while (_Exponent != _Source.exponent_bias_bitset());
	} else if (_Exponent < _Source.exponent_bias_bitset()) {
		do {
			_Dest.divev2();
			_Exponent += _One;
		} while (_Exponent != _Source.exponent_bias_bitset());
	}

	// set sign
	_Dest.negative = (_Source.bitset() & _Source.sign_mask()) != _Zero;

	// set precision
	size_t desired_significant_count = 
		static_cast<size_t>( ::ceil(::log(2)/::log(10) * (_Source.significand_bits+1)) );
	_Dest.reset_precision( desired_significant_count );
	return _Dest.to_string();
}

template<size_t m, size_t e, bool _Isbase> inline
std::ostream& operator<<(std::ostream& _Ostr, floatX<m,e,_Isbase> _Fp) {
	return (_Ostr << to_string(_Fp));
}

}// namespace calculation


#include <limits>
namespace std 
  {
template<size_t m, size_t e, bool _IsOpt>
class numeric_limits< calculation::floatX<m,e,_IsOpt> > : public _Num_float_base{ // numeric limits for arbitrary type _Ty (say little or nothing)
public:
  _NODISCARD static constexpr calculation::floatX<m,e,_IsOpt>(min)() noexcept {
    abort();
	}

  _NODISCARD static constexpr calculation::floatX<m,e,_IsOpt>(max)() noexcept {
    abort();
	}

  _NODISCARD static constexpr calculation::floatX<m,e,_IsOpt> lowest() noexcept {
    abort();
	}

  _NODISCARD static calculation::floatX<m,e,_IsOpt> epsilon() noexcept {
		return calculation::floatX<m,e,_IsOpt>( calculation::floatX<m,e,_IsOpt>::epsilon_bitset() );
	}

  _NODISCARD static calculation::floatX<m,e,_IsOpt> round_error() noexcept {
    return epsilon() / 2;
	}

  _NODISCARD static constexpr calculation::floatX<m,e,_IsOpt> denorm_min() noexcept {
    abort();
	}

  _NODISCARD static calculation::floatX<m,e,_IsOpt> infinity() noexcept {
    return calculation::floatX<m,e,_IsOpt>( calculation::floatX<m,e,_IsOpt>::inf_bitset() );
	}

  _NODISCARD static calculation::floatX<m,e,_IsOpt> quiet_NaN() noexcept {
		return calculation::floatX<m,e,_IsOpt>( calculation::floatX<m,e,_IsOpt>::quiet_NaN_bitset() );
	}

  _NODISCARD static calculation::floatX<m,e,_IsOpt> signaling_NaN() noexcept {
		return calculation::floatX<m,e,_IsOpt>( calculation::floatX<m,e,_IsOpt>::signaling_NaN_bitset() );
	}

  static constexpr size_t digits = calculation::floatX<m,e,_IsOpt>::significand_bits + 1;
};
  }

/**
#include "float.h"
#include <iostream>

int main() {
  using namespace::calculation;

  std::cout.precision(16);
  for (size_t i = 0; i < 10000; i++) {
    double a = rand() + rand() / double(RAND_MAX);
    double b = (rand() + rand() / double(RAND_MAX));
    double c = a / b;
    float128_t c1 = float128_t(a) / float128_t(b);
    std::cout << c << "\t" << c1 << std::endl;

    // not( relative error <= epsilon()/2 )
    if ( !(abs(c1 - c) <= (std::numeric_limits<double>::epsilon()/2)*c) ) {
      std::cout << "Error" << std::endl;
      std::cout << "error = " << abs(c1 - c) << std::endl;
      std::cout << "eps   = " << c * float128_t(std::numeric_limits<double>::epsilon()/2) << std::endl;
      std::cin.get();
    }
  }
  
  return 0;
}
*/

/*Example:{
#include "clmagic/calculation/fundamental/float.h"
#include <iostream>

int main(int, char**) {
    using namespace::calculation;

    std::cout << to_string(float256_t(1.5)) << std::endl;
    std::cout << to_string(float256_t(1.5) * float256_t(2.0)) << std::endl;
    std::cout << to_string(float256_t(5.5) * float256_t(129.0)) << std::endl << std::endl;

    std::cout << to_string(float128_t(1.0) / float128_t(3.0)) << std::endl;
    std::cout << to_string(float128_t(1.0) / float128_t(9.0)) << std::endl;
    std::cout.precision(16);
    std::cout << double(float128_t(1.0) / (float128_t(9) + (float128_t(9) / float128_t(10)))) << std::endl;
    std::cout << to_string(float128_t(1.0) / (float128_t(9) + (float128_t(9) / float128_t(10)))) << std::endl << std::endl;

    std::cout << to_string(float256_t(1.0) / float256_t(3.0)) << std::endl;
    std::cout << to_string(float256_t(1.0) / float256_t(9.0)) << std::endl;
    std::cout << to_string(float256_t(1.0) / (float256_t(9) + (float256_t(9) / float256_t(10)))) << std::endl << std::endl;

    std::cout << to_string(float512_t(1.0) / float512_t(3.0)) << std::endl;
    std::cout << to_string(float512_t(1.0) / float512_t(9.0)) << std::endl;
    std::cout << to_string(float512_t(1.0) / (float512_t(9) + (float512_t(9) / float512_t(10)))) << std::endl;

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
    Spectrum<float128_t, 10> rod = {
    float128_t(1.0), float128_t(9)/10, float128_t(6)/10, float128_t(4)/10, float128_t(9)/10,
    float128_t(9)/10, float128_t(4)/10, float128_t(4)/10, float128_t(2)/10, float128_t(1)/10
    };
      
    for (float i = 385.f; i < 740.0f; i += 10.0f) {
    std::cout << i << " = " << to_string(rod(float128_t(i))) << std::endl;
    }
	
    return 0;
}
}*/