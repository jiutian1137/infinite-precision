/**
 * @license
 *   Please identify Author, 2020 - 2021
 * @author 
 *   LongJiangnan, Jiang1998Nan@outlook.com
 * @breif
 *   infinite precision calculation for floating number
 * @readme
 *   FloatX<..>::frexp and ldexp are not perfect, because we don't have memory-operation functions.
 *   In next optimization, we will add dynamic something and will have memory-operation functions.
 * 
 *   FloatX<..>::cvtui and operator Integer() are not perfect.
*/
#pragma once

#include <stdexcept>
#include <cassert>

#include <bitset>
namespace std {
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
} // namespace std::bitset extension

#include <cmath>

namespace calculation {
/**
 * the rounding has only two kinds, 
 *   one is round-up, 
 *   another is round-down.
 * 
 * so round-error has two kinds,
 *   one is abs(round_up_number - origin_number)
 *   another is abs(round_down_number - origin_number)
 * 
 * @nearest-even
 *   min( abs(round_up_number - origin_number), abs(round_down_number - origin_number) )
 * 
 * @structure
 *     1.1111111 11111111 11111111           :IEEE-32bit significand
 *      ...
 *      ... do some arthmetic
 *   = 1.0001011 00011111 00000111 10000100  :need smaller_shift 8 places
 * 
 *                                  sticky_bits
 *                            back_bit |
 *                               |  |--+-...
 *     1.0001011 00011111 00000111 10000100  :shifted
 *                                 |
 *                             round-bit
 * 
 * @example
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
 * 
 * @analysis
 *   we think back on this calculation, we can see that
 *     if round-bit is '1'
 * 
 *       round-up-error <= epsilon/2
 *       abs(1000'0000 - 0111'1000) <= 0001'0000/2
 *                        0000'1000 <= 0000'1000
 * 
 *       round-down-error >= epsilon/2
 *       abs(0111'0000 - 0111'1000) >= 0001'0000/2
 *                        0000'1000 >= 0000'1000
 * 
 *       we can't choose
 * 
 *       if sticky-bits not 'all(0)'
 *         round-up-error < epsilon/2
 *         abs(1000'0000 - 0111'1000'0100) < 0001'0000/2
 *                          0000'0111'1100 < 0000'1000
 * 
 *         round-down-error > epsilon/2
 *         abs(0111'0000 - 0111'1000'0100) > 0001'0000/2
 *                          0000'1000'0100 > 0000'1000
 *         so we round-up.
 *   other
 *     we round-down
 * 
 * @reference "https://www.slideshare.net/prochwani95/06-floating-point"
*/
template<typename Bitset>
bool need_roundup_for_nearest_even(const Bitset& significand, size_t smaller_shift,
  Bitset& sticky_mask) 
{
  if ( smaller_shift != 0/* && smaller_shift - 1 < significand.size()*/ ) {
    bool round_bit = significand.test( smaller_shift - 1 );
    bool sticky_bit = false;
    if ( smaller_shift != 1 ) {
      sticky_mask = 1;
      sticky_mask <<= (smaller_shift - 1);
      sticky_mask -= 1;
      sticky_bit = (significand & sticky_mask) != 0;
    }
  
    if ( round_bit && (sticky_bit/* | significand.test(smaller_shift)*/) ) {
      return true;
    }
  }

  return false;
}

/** 
 * @brief TEMPLATE infinite precision floating_point
 *   digit: (1 + 0.Significand) * 2^Exponent * (-1)^Sign, 
 *   Bits: Sigbits + Expbits + 1 = Bits.
 * @reference 
 *   "https://www.rfwireless-world.com/Tutorials/floating-point-tutorial.html",
 *   "https://github.com/gcc-mirror/gcc/blob/master/gcc/real.c",
 *   "https://www.exploringbinary.com/binary-division/"
*/
template<size_t s, size_t e, typename Bitset = std::bitset<s + e + 1>, 
  bool IsOpt = false>
class FloatX {
public:
  static_assert(2 <= s);
  static_assert(2 <= e && e <= sizeof(size_t)*8);

  using bitset_type = Bitset;

  static constexpr size_t significand_bits = s;
  static constexpr size_t exponent_bits = e;
  static constexpr size_t sign_bits = 1;
  static constexpr size_t bits = significand_bits + exponent_bits + sign_bits;
    
  static constexpr size_t significand_offset = 0;
  static constexpr size_t exponent_offset = significand_offset + significand_bits;
  static constexpr size_t sign_offset = exponent_offset + exponent_bits;

  static Bitset significand_mask() {
    static Bitset _Mask =
      (Bitset(1) << significand_bits) - 1
        /* << 0 */;
    return _Mask;
  }

  static Bitset exponent_mask() {
    static Bitset _Mask =
      ( (Bitset(1) << exponent_bits) - 1 )
        << exponent_offset;
    return _Mask;
  }

  static Bitset sign_mask() {
    static Bitset _Mask =
      Bitset(1) << sign_offset;
    return _Mask;
  }

  static Bitset hidden_significand_mask() {
    static Bitset _Mask =
      Bitset(1) << significand_bits;
    return _Mask;
  }

  static Bitset exponent_bias_bitset() {
    static Bitset _Exp2_zero =
      ( (Bitset(1) << (exponent_bits - 1)) - 1 )
        << exponent_offset;
    return _Exp2_zero;
  }

  static Bitset possign_bitset() {
    return zero_bitset();
  }

  static Bitset negsign_bitset() {
    return sign_mask();
  }

  static Bitset zero_bitset() {
    return Bitset(0);
  }

  static Bitset inf_bitset() {
    // inf = full 1 in exponent_bitset
    return exponent_mask();
  }

  static Bitset quiet_NaN_bitset() {
    static Bitset _Quiet_nan =
      exponent_mask() 
        | (Bitset(1) << (significand_bits - 1)) ;
    return _Quiet_nan;
  }

  static Bitset signaling_NaN_bitset() {
    static Bitset _Signaling_nan =
      exponent_mask() 
        | (Bitset(1) << (significand_bits - 1))
          | Bitset(1);
    return _Signaling_nan;
  }

  static Bitset epsilon_bitset() {
    // epsilon = exp2(0 - mantissa_bits) 
    static Bitset _Epsilon =
      exponent_bias_bitset() - (Bitset(significand_bits) << exponent_offset);
    return _Epsilon;
  }

  static bool do_iszero(const Bitset& x) {
    return (x&(exponent_mask()|significand_mask())) == zero_bitset();
  }

  static bool do_isinf(const Bitset& x) {
    return (x&exponent_mask()) == exponent_mask();
  }

  static bool do_isnan(const Bitset& x) {
    return ((x&exponent_mask()) == exponent_mask()) 
      && ((x&significand_mask()) != zero_bitset());
  }

  static bool do_equal(const Bitset& left, const Bitset& right) {
    return (left == right) || (do_iszero(left) && do_iszero(right));
  }

  static bool do_less(const Bitset& left, const Bitset& right) {
    const Bitset SignMask = sign_mask();
    const Bitset SigExpMask = ~SignMask;
    const Bitset Zero = zero_bitset();
    // Check equal
    Bitset left_sigexp = left & SigExpMask;
    Bitset right_sigexp = right & SigExpMask;
    if ((left == right) || (left_sigexp == Zero && right_sigexp == Zero)) {
      return false;
    } else {

      const Bitset Positive = possign_bitset();
      const Bitset Negative = negsign_bitset();
      // Check sign
      Bitset left_sign = left & SignMask;
      Bitset right_sign = right & SignMask;
      if (left_sign != right_sign) {
        return left_sign == Negative && right_sign == Positive;
      } else {

        // Check exponent|significand
        return (left_sign == Positive) ? left_sigexp < right_sigexp
          : right_sigexp < left_sigexp;

      }

    }
  }

  static void do_plus(const Bitset& left, const Bitset& right, Bitset& result) {
    const Bitset SignMask = sign_mask();
    const Bitset SigExpMask = ~SignMask;
    const Bitset Zero = zero_bitset();
    // 0. Check special case '0'
    if ((left & SigExpMask) == Zero) {
      result = right;
      return;
    } else if ((right & SigExpMask) == Zero) {
      result = left;
      return;
    }

    const Bitset ExpMask = exponent_mask();
    const Bitset Inf = inf_bitset();
    const Bitset NaN = quiet_NaN_bitset();
    // 0. Check special case 'Inf' and 'NaN'
    if ((left & ExpMask) == Inf || (right & ExpMask) == Inf) {
      result = NaN;
      return;
    }

    const Bitset SigMask = significand_mask();
    const Bitset HdMask = hidden_significand_mask();
    // 1. Decode into significand-bitset, exponent-bitset, sign-bitset
    Bitset left_significand = (left & SigMask) | HdMask;
    Bitset left_exponent = left & ExpMask;
    Bitset left_sign = left & SignMask;
    Bitset right_significand = (right & SigMask) | HdMask;
    Bitset right_exponent = right & ExpMask;
    Bitset right_sign = right & SignMask;
  
    /**
     *  11111111 = MAX(8bits)
     *  11111111 = MAX(8bits)
     * 111111110 = MAX(8bits) + MAX(8bits)
     * So we only reserved '1' bits
    */
    const size_t ReservedBits = 1;
    const size_t HdsigBits = significand_bits + 1;
    const size_t Bits = bits;
    assert(HdsigBits + ReservedBits <= Bits);
    // 2. Extend to more accuracy for trunc_error
    size_t shift = 0;
    if ( HdsigBits + ReservedBits < Bits ) {
      shift = Bits - (HdsigBits + ReservedBits);
      left_significand <<= shift;
      right_significand <<= shift;
    }

    const size_t ExpOffset = exponent_offset;
    // 3. Sync exponent-bitset to greater
    if ( left_exponent < right_exponent ) 
    {
      Bitset exponent_difference = right_exponent - left_exponent;
      Bitset temp = exponent_difference >> ExpOffset;
      size_t smaller_shift = reinterpret_cast<const size_t&>(temp);
      assert(smaller_shift != 0);
      if (smaller_shift >= HdsigBits+1 || left_exponent < exponent_difference) {
        /**
         * 00000000 1.1111111 11111111 11111111 = (left&SigMask)|HdsigMask
         * 00000000 0 0000000 00000000 00000000 1.11111 = ((right&SigMask)|HdsigMask) >> HdsigBits
         * We also need to compute and rounding
         * 
         * 00000000 1.1111111 11111111 11111111 = (left&SigMask)|HdsigMask
         * 00000000 0 0000000 00000000 00000000 01.11111 = ((right&SigMask)|HdsigMask) >> (HdsigBits+1)
         * Not need to compute
         * 
         * 00000000 1.1111111 11111111 11111111 = (left&SigMask)|HdsigMask
         *                             INFSmall = ((right&SigMask)|HdsigMask) >> smaller_shift
         * Not need to compute
        */
        result = right;
        return;
      }

      if (need_roundup_for_nearest_even(left_significand, smaller_shift, temp)) {
        left_significand >>= smaller_shift;
        left_significand += 1;
      } else {
        left_significand >>= smaller_shift;
      }

      left_exponent += exponent_difference;
    } 
    else if ( left_exponent != right_exponent ) 
    {
      Bitset exponent_difference = left_exponent - right_exponent;
      Bitset temp = exponent_difference >> ExpOffset;
      size_t smaller_shift = reinterpret_cast<const size_t&>(temp);
      assert(smaller_shift != 0);
      if (smaller_shift >= HdsigBits+1 || right_exponent < exponent_difference) {
        result = left;
        return;
      }

      if (need_roundup_for_nearest_even(right_significand, smaller_shift, temp)){
        right_significand >>= smaller_shift;
        right_significand += 1;
      } else {
        right_significand >>= smaller_shift;
      }

      right_exponent += exponent_difference;
    }

    // 4. Plus significand-bitset, update sign-bits
    if ( left_sign == right_sign ) 
    {
      left_significand += right_significand;
    } 
    else 
    {
      if (left_significand == right_significand) {
        result = Zero;
        return;
      } else if (left_significand < right_significand) {
        left_significand = right_significand - left_significand;
        left_sign = right_sign;
      } else  /* left_significand > right_significand */{
        left_significand -= right_significand;
        /*left_sign = left_sign;*/
      }
    }

    // 5. Normalize significand-bits and exponent-bits
    if ( (left_significand & (~SigMask)) != 0 )
    {
      // normalize for select 'shift'
      size_t smaller_shift = 0;
      while ( !left_significand.test((Bits - 1) - smaller_shift) )
        ++smaller_shift;
      assert(Bits - HdsigBits >= smaller_shift);
      smaller_shift = (Bits - HdsigBits) - smaller_shift;

      // normalize for 'significand' with 'round-up'
      Bitset temp;
      if (need_roundup_for_nearest_even(left_significand, smaller_shift, temp)) {
        left_significand >>= smaller_shift;
        left_significand += 1;
        if (left_significand.test(HdsigBits)) {
          /**
           * 00000000 1.1111111 11111111 11111111 = MAX(float32_significand), only case adding one bit after carried
           * 00000001 0.0000000 00000000 00000000 = MAX(float32_significand) carried
           *                                   ||
           *                                   |next round bit must be '0', not carry again
           *                                next first bit
          */
          left_significand >>= 1;
          smaller_shift += 1;
        }
      } else {
        left_significand >>= smaller_shift;
      }

      // normalize for 'exponent'
      if ( smaller_shift > shift ) {
        Bitset exponent_difference = Bitset(smaller_shift - shift) << ExpOffset;
        if (left_exponent >= Inf - exponent_difference) {/** special case 'infinite large' */
          result = Inf|left_sign;
          return;
        }
        left_exponent += exponent_difference;
      } else {
        Bitset exponent_difference = Bitset(shift - smaller_shift) << ExpOffset;
        if (left_exponent < exponent_difference) {/** special case 'infinite small' */
          result = Zero;
          return;
        }
        left_exponent -= exponent_difference;
      }
    } 
    else 
    {
      // normalize for select 'shift'
      size_t larger_shift = 0;
      while ( !left_significand.test((HdsigBits - 1) - larger_shift) )
        ++larger_shift;

      // normalize for 'significand'
      left_significand <<= larger_shift;

      // normalize for 'exponent'
      Bitset exponent_difference = Bitset(shift + larger_shift) << ExpOffset;
      if (left_exponent < exponent_difference) {/** special case 'infinite small' */
        result = Zero;
        return;
      }
      left_exponent -= exponent_difference;
    }

    // 6. Encode ...
    result = (left_significand&SigMask)|left_exponent|left_sign;
  }

  static void do_multiply(const Bitset& left, const Bitset& right, Bitset& result) {
    const Bitset SignMask = sign_mask();
    const Bitset SigExpMask = ~SignMask;
    const Bitset Zero = zero_bitset();
    // 0. Check special case '0'
    if ( (left & SigExpMask) == Zero || (right & SigExpMask) == Zero ) {
      result = Zero;
      return;
    }

    const Bitset SigMask = significand_mask();
    const Bitset ExpMask = exponent_mask();
    const Bitset Inf = inf_bitset();
    const Bitset NaN = quiet_NaN_bitset();
    // 0. Check special case 'Inf' and 'NaN'
    if ( (left & ExpMask) == Inf || (right & ExpMask) == Inf ) {
      if ( (left&SigMask) != 0 || (right&SigMask) != 0 ) {
        result = NaN;
        return;
      }

      result = Inf|((left&SignMask)^(right&SignMask));
      return;
    }

    const Bitset HdMask = hidden_significand_mask();
    const size_t SigBits = significand_bits;
    const size_t HdsigBits = significand_bits + 1;
    /**
     *          11111111 = MAX(8bits)
     *          11111111 = MAX(8bits)
     * 11111110 00000001 = MAX(8bits)*MAX(8bits)
     * we need double bits
    */
    const size_t DoubleHdsigBits = (significand_bits + 1)*2;
    using BigBitset = std::bitset<(significand_bits + 1)*2> ;
    // 1. Decode ...
    BigBitset left_significand = 0;
    reinterpret_cast<
    Bitset&>(left_significand) = (left&SigMask)|HdMask;
    Bitset left_exponent = left&ExpMask;
    BigBitset right_significand = 0;
    reinterpret_cast<
    Bitset&>(right_significand) = (right&SigMask)|HdMask;
    Bitset right_exponent = right&ExpMask;

    // 2. result_significand = (left_significand>>SigBits * right_significand>>SigBits) << SigBits = left_significand*right_significand >> SigBits = left_significand*right_significand >> SigBits
    left_significand *= right_significand;
    const size_t shift = SigBits;
  
    const Bitset ExpBias = exponent_bias_bitset();
    // 3. result_exponent = (left_exponent-ExpBias)+(right_exponent-ExpBias) + ExpBias = left_exponent+right_exponent - ExpBias
    if ( right_exponent > ExpBias) {
      right_exponent -= ExpBias;
      if (left_exponent >= Inf - right_exponent) {
        result = Inf|((left&SignMask)^(right&SignMask));
        return;
      }
      left_exponent += right_exponent;
    } else if ( right_exponent < ExpBias) {
      right_exponent = ExpBias - right_exponent;
      if (left_exponent < right_exponent) {
        result = Zero;
        return;
      }
      left_exponent -= right_exponent;
    }
  
    // 4. Normalize, assert( this_significand&(~significand_mask) != 0 );
    {
      // Select 'shift'
      size_t smaller_shift = 0;
      while ( !left_significand.test((DoubleHdsigBits - 1)- smaller_shift) )
        ++smaller_shift;
      assert(DoubleHdsigBits - HdsigBits >= smaller_shift);
      smaller_shift = (DoubleHdsigBits - HdsigBits) - smaller_shift;

      // Apply 'shift' to 'significand' with 'round-up'
      BigBitset tmp;
      if ( need_roundup_for_nearest_even(left_significand, smaller_shift, tmp) ) {
        left_significand >>= smaller_shift;
        left_significand += 1;
        if ( left_significand.test(HdsigBits) ) {
          left_significand >>= 1;
          smaller_shift += 1;
          /* assert(not need carry); @see do_plus(..) */
        }
      } else {
        left_significand >>= smaller_shift;
      }

      const size_t ExpOffset = exponent_offset;
      // Apply 'shift' to 'exponent'
      if ( smaller_shift > shift ) {
        Bitset exponent_difference = Bitset(smaller_shift - shift) << ExpOffset;
        if (left_exponent >= Inf - exponent_difference) {
          result = Inf|((left&SignMask)^(right&SignMask));
          return;
        }
        left_exponent += exponent_difference;
      } else {
        Bitset exponent_difference = Bitset(shift - smaller_shift) << ExpOffset;
        if (left_exponent < exponent_difference) {
          result = Zero;
          return;
        }
        left_exponent -= exponent_difference;
      }
    }

    // 5. Encode ...
    result = (reinterpret_cast<const Bitset&>(left_significand)&SigMask)|left_exponent|((left&SignMask)^(right&SignMask));
  }

  static void do_divide(const Bitset& left, const Bitset& right, Bitset& result) {
    const Bitset SignMask = sign_mask();
    const Bitset SigExpMask = ~SignMask;
    const Bitset Zero = zero_bitset();
    const Bitset Inf = inf_bitset();
    const Bitset NaN = quiet_NaN_bitset();
    // 0. Check special case '0'
    if ((left & SigExpMask) == Zero && (right & SigExpMask) == Zero) {/* 0/0 = NaN */
      result = NaN;
      return;
    } else if ( (left & SigExpMask) == Zero ) {/* 0/X = 0 */
      result = Zero;
      return;
    } else if ( (right & SigExpMask) == Zero ) {/* X/0 = Inf */
      result = Inf|((left&SignMask)^(right&SignMask));
      return;
    }

    const Bitset ExpMask = exponent_mask();
    // 0. Check special case 'Inf' and 'Nan'
    if ((left & ExpMask) == Inf || (right & ExpMask) == Inf) {
      result = NaN;
      return;
    }

    const Bitset SigMask = significand_mask();
    const Bitset HdMask = hidden_significand_mask();
    // 1. Decode into divide-operator, exponent-bitset
    Bitset dividend = (left & SigMask) | HdMask;
    Bitset left_exponent = left & ExpMask;
    Bitset divisor = (right & SigMask) | HdMask;
    Bitset right_exponent = right & ExpMask;
  
    const size_t ReservedBits = 1;/* @see (step 3.) */
    const size_t Bits = bits;
    const size_t HdsigBits = significand_bits + 1;
    // 2. Extend more accuracy
    assert(HdsigBits + ReservedBits <= Bits);
    size_t shift = 0;
    if (HdsigBits + ReservedBits < Bits) {
      shift = Bits - (HdsigBits + ReservedBits);
      dividend <<= shift;
      divisor <<= shift;
    }

    // 3. Divide significand
    Bitset left_significand = 0;
    /*Bitset divstep = Bitset(1) << ((HdsigBits - 1) + shift);*/
    Bitset divstep = Bitset(1) << (Bits-ReservedBits - 1);
    do {
      if ( dividend >= divisor ) {
        left_significand |= divstep;
        dividend -= divisor;
      }
      dividend <<= 1;/* if(dividend >= divisor == false), here need reserved 1 bits. */
      divstep >>= 1;
    } while (divstep != 0);

    const Bitset ExpBias = exponent_bias_bitset();
    // 4. result_exponent = (left_exponent-ExpBias)-(right_exponent-ExpBias) + ExpBias = left_exponent-right_exponent + ExpBias
    if (right_exponent > ExpBias) {
      right_exponent -= ExpBias;
      if (left_exponent < right_exponent) {
        result = Zero;
        return;
      }
      left_exponent -= right_exponent;
    } else if (right_exponent < ExpBias) {
      right_exponent = ExpBias - right_exponent;
      if (left_exponent >= Inf - right_exponent) {
        result = Inf|((left&SignMask)^(right&SignMask));
        return;
      }
      left_exponent += right_exponent;
    }

    // 5. Normalize, see(step 2.)
    if (shift >= 1 || (left_significand&(~SigMask)) != 0) 
    {
      // Select 'shift'
      size_t smaller_shift = 0;
      while ( !left_significand.test((Bits - 1) - smaller_shift) )
        ++smaller_shift;
      assert(Bits - HdsigBits >= smaller_shift);
      smaller_shift = (Bits - HdsigBits) - smaller_shift;


      // Apply 'shift' to 'significand' with 'round-up'
      Bitset tmp;
      if ( need_roundup_for_nearest_even(left_significand, smaller_shift, tmp) ) {
        left_significand >>= smaller_shift;
        left_significand += 1;
        if ( left_significand.test(HdsigBits) ) {
          left_significand >>= 1;
          smaller_shift += 1;
          /* assert(not need carry); @see do_plus(..) */
        }
      } else {
        left_significand >>= smaller_shift;
      }

      // Apply 'shift' to 'exponent'
      assert( smaller_shift <= shift );/* smaller_shift = shift-1|shift */
      Bitset exponent_difference = Bitset(shift - smaller_shift) << exponent_offset;
      if ( left_exponent >= Inf - exponent_difference ) {
        result = Zero;
        return;
      }
      left_exponent -= exponent_difference;
    } 
    else 
    {
      abort();
    }

    // 6. Encode ...
    result = (left_significand&SigMask) | left_exponent | ((left&SignMask)^(right&SignMask));
  }

  static void do_modf(const Bitset& source, Bitset& fpart, Bitset& ipart) {
    const Bitset ExpBias = exponent_bias_bitset();
    const Bitset ExpMask = exponent_mask();
    const Bitset Zero = zero_bitset();
    // Check only fraction, 1.010101010111... * exp2(0), 
    Bitset exponent = source & ExpMask;
    if ( exponent < ExpBias ) {
      fpart = source;
      ipart = Zero;
      return;
    }

    /**
     *                   00000000 1.0000000 00000000 00000000 = source
     * 10000000 00000000 00000000 0                           = source << 23
     * significand are all integer parts
    */
    const size_t ExpOffset = exponent_offset;
    const Bitset exp2_SigBits = Bitset(significand_bits) << exponent_offset;
    // Check only integer, 1010101010111... * exp2(exponent - significand_bits)
    Bitset truncated_exponent = exponent - exp2_SigBits;
    if ( truncated_exponent >= ExpBias ) {
      ipart = source;
      fpart = Zero;
      return;
    }

    // Compute Mask
    Bitset trunc_exponent = ExpBias - truncated_exponent;
    trunc_exponent >>= ExpOffset;
    size_t trunc_bits = reinterpret_cast<const size_t&>(trunc_exponent);
    Bitset frac_mask = (( Bitset(1) << trunc_bits ) - 1);
    Bitset trunc_mask = ~frac_mask;
    
    // Apply Mask
    ipart = source & trunc_mask;
    fpart = source & frac_mask;
    if (fpart == 0) {
      return;
    } else {
      size_t larger_shift = 0;
      while ( !fpart.test((trunc_bits - 1) - larger_shift) )
        ++larger_shift;
      larger_shift += (significand_bits+1 - trunc_bits);

      const Bitset SigMask = significand_mask();
      const Bitset SignMask = sign_mask();
      fpart <<= larger_shift;
      fpart &= SigMask;
      fpart |= (exponent - (Bitset(larger_shift)<<ExpOffset));
      fpart |= (source & SignMask);
    }
  }

  template<typename Integer>
  static void do_frexp(const Bitset& source, Bitset& base, Integer& exp2) {
    static_assert(sizeof(Integer)*8 >= e);

    const Bitset SignMask = sign_mask();
    const Bitset SigExpMask = ~SignMask;
    const Bitset Zero = zero_bitset();
    // Check special case '0'
    if ((source & SigExpMask) == Zero) {/* 0*exp2(0) = 0 */
      base = Zero;
      exp2 = 0;
      return;
    }

    const size_t ExpOffset = exponent_offset;
    const Bitset SigMask = significand_mask();
    const Bitset ExpMask = exponent_mask();
    const Bitset Inf = inf_bitset();
    const Bitset NaN = quiet_NaN_bitset();
    // Check special case 'Inf' and 'NaN'
    if ( (source & ExpMask) == Inf) {
      if ( (source & SigMask) != 0 ) {/* NaN*exp2(undefined) = NaN */
        base = NaN;
        return;
      }

      base = NaN;
      exp2 = 0;
      Bitset tmp = Inf >> ExpOffset;
      std::copy_n((const char*)(&tmp), e/8 + size_t(e%8!=0?1:0), (char*)(&exp2));
      return;
    }

    // Decode to 'source' = 'base' * exp2('exp2')
    const Bitset ExpBias = exponent_bias_bitset();
    base = source & SigMask;
    base |= ExpBias;
    base |= (source & SignMask);
    Bitset tmp = source & ExpMask;
    if (tmp >= ExpBias) {
      tmp -= ExpBias;
      tmp >>= ExpOffset;
      exp2 = 0;
      std::copy_n((const char*)(&tmp), e/8 + size_t(e%8!=0?1:0), (char*)(&exp2));
    } else {
      tmp = ExpBias - tmp;
      tmp >>= ExpOffset;
      exp2 = 0;
      std::copy_n((const char*)(&tmp), e/8 + size_t(e%8!=0?1:0), (char*)(&exp2));
      exp2 = -exp2;
    }
  }

  template<typename Integer>
  static void do_ldexp(const Bitset& base, const Integer& exp2, Bitset& destination) {
    const Bitset SignMask = sign_mask();
    const Bitset SigExpMask = ~SignMask;
    const Bitset Zero = zero_bitset();
    // Check special case '0'
    if ((base & SigExpMask) == Zero) {
      destination = Zero;
      return;
    }

    const Bitset ExpMask = exponent_mask();
    const Bitset Inf = inf_bitset();
    const Bitset NaN = quiet_NaN_bitset();
    // Check special case 'Inf' and 'NaN'
    if ((base & ExpMask) == Inf) {
      destination = NaN;
      return;
    }

    // Check special case "exp2(0)"
    if (exp2 == 0) {
      destination = base;
      return;
    }

    Bitset exponent;

    const size_t ExpOffset = exponent_offset;
    // Extract unsigned exponent
    Integer uexp2 = exp2 > 0 ? exp2 : -exp2;
    if constexpr (sizeof(Integer) < sizeof(Bitset)) {
      exponent = 0;
      std::copy_n(((const char*)&uexp2), sizeof(Integer), (char*)(&exponent));
      exponent <<= ExpOffset;
    } else {
      Integer tmp = uexp2 << ExpOffset;
      // assert(tmp < Inf/2);
      exponent = reinterpret_cast<const Bitset&>(tmp);
    }

    // Extract sign of exponent, and combination
    if (exp2 >= 0) {
      destination = base + exponent;
    } else {
      destination = base - exponent;
    }
  }

  template<typename InBitset>
  static void do_cvtui(const InBitset& source, Bitset& destination) {
    const Bitset Zero = zero_bitset();
    // 0. Check '0'
    if ( source == 0 ) {
      destination = Zero;
      return;
    }

    const size_t InBits = sizeof(InBitset) * 8;
    const size_t ExpBits = exponent_bits;
    const Bitset Inf = inf_bitset();
    // 0. Check 'super large' (beyond exponent_bits)
    if (InBits > ExpBits && (source >> (size_t(1)<<ExpBits)) != 0) {
      destination = Inf;
      return;
    }

    const Bitset ExpBias = exponent_bias_bitset();
    const size_t ExpOffset = exponent_offset;
    const size_t SigBits = significand_bits;
    // 1. source = 10100100010.00000000000000000 = 000000000.00000000000010100100010 << 23 = source.as<Significand>()|(23<<ExpOffset + ExpBias)
    size_t sshift = SigBits;
    Bitset dest_exponent = ExpBias;
    Bitset dest_significand = 0;

    const Bitset SigMask = significand_mask();
    const size_t HdsigBits = significand_bits + 1;
    const size_t HdsigBytes = HdsigBits/8 + (HdsigBits%8!=0?1:0);
    // 2. Normalize
    Bitset tmp = ~SigMask;
    if (InBits >= HdsigBits && (source & reinterpret_cast<const InBitset&>(tmp)) != 0) 
    {
      // normalize for select 'shift'
      size_t smaller_shift = 0;
      while ( !source.test(InBits-1 - smaller_shift) )
        ++smaller_shift;
      assert(InBits - HdsigBits >= smaller_shift);
      smaller_shift = (InBits - HdsigBits) - smaller_shift;

      // normalize for 'significand' with 'round-up'
      InBitset tmp2;
      if ( need_roundup_for_nearest_even(source, smaller_shift, tmp2) ) {
        assert(InBits > HdsigBits);
        tmp2 = source >> smaller_shift;
        std::memcpy(&dest_significand, &tmp2, std::min(sizeof(InBitset), HdsigBytes));
        dest_significand += 1;
        if (dest_significand.test(HdsigBits)) {
          dest_significand >>= 1;
          smaller_shift += 1;
          /* assert(not need carry), @see do_plus(..) */
        }
      } else {
        tmp2 = source >> smaller_shift;
        std::memcpy(&dest_significand, &tmp2, std::min(sizeof(InBitset), HdsigBytes));
      }

      // normalize for 'exponent'
      Bitset exponent_difference = Bitset(smaller_shift+sshift) << ExpOffset;
      if (exponent_difference > Inf || dest_exponent >= Inf - exponent_difference) {/** special case 'infinite large' */
        destination = Inf;
        return;
      }
      dest_exponent += exponent_difference;
    } 
    else 
    {
      // normalize for select 'shift'
      size_t larger_shift = 0;
      if (InBits < HdsigBits)
        larger_shift += HdsigBits - InBits;
      while ( !source.test((HdsigBits - 1) - larger_shift) )
        ++larger_shift;

      // normalize for 'significand'
      dest_significand = 0;
      std::memcpy(&dest_significand, &source, std::min(sizeof(InBitset), HdsigBytes));
      dest_significand <<= larger_shift;

      // normalize for 'exponent', integer not have "infinite small"
      if (larger_shift >= sshift) {
        Bitset exponent_difference = Bitset(larger_shift - sshift) << ExpOffset;
        dest_exponent -= exponent_difference;
      } else {
        Bitset exponent_difference = Bitset(sshift - larger_shift) << ExpOffset;
        dest_exponent += exponent_difference;
      }
    }

    destination = (dest_significand&SigMask)|dest_exponent;
  }

  template<size_t s2, size_t e2, typename InBitset>
  static void do_cvtf(const InBitset& source, Bitset& destination) {
    typedef FloatX<s2, e2, InBitset> InFloat;
    typedef FloatX<s, e, Bitset> Float;

    const InBitset InSignMask = InFloat::sign_mask();
    const InBitset InZero = InFloat::zero_bitset();
    const Bitset Zero = Float::zero_bitset();
    // Check '0'
    if ((source&(~InSignMask)) == InZero) {
      destination = Zero;
      return;
    }
    
    const Bitset Inf = Float::inf_bitset();
    const Bitset NaN = Float::quiet_NaN_bitset();
    const Bitset Positive = Float::possign_bitset();
    const Bitset Negative = Float::negsign_bitset();
    const InBitset InSigMask = InFloat::significand_mask();
    const InBitset InExpMask = InFloat::exponent_mask();
    const InBitset InInf = InFloat::inf_bitset();
    const InBitset InPositive = InFloat::possign_bitset();
    // Check 'Inf' of 'NaN'
    if ((source&InExpMask) == InInf) {
      if ((source&InSigMask) != 0) {
        destination = NaN;
        return;
      }

      destination = Inf|(((source&InSignMask) == InPositive) ? Positive : Negative);
      return;
    }

    const size_t SigBits = Float::significand_bits;
    const size_t SigBytes = s/8 + (s%8!=0?1:0);
    const size_t InSigBits = InFloat::significand_bits;
    const size_t InSigBytes = InSigBits/8 + (InSigBits%8!=0?1:0);
    // Copy significand_bits except hidden_bit (named mantissa_bits)
    InBitset source_mantissa = source & InSigMask;
    Bitset destination_mantissa = 0;
    if (SigBits < InSigBits) 
    {/* float64 to float32 */
      InBitset source_shifted_mantissa = source_mantissa >> (InSigBits - SigBits);
      std::memcpy(&destination_mantissa, &source_shifted_mantissa, SigBytes);
    } 
    else if (SigBits != InSigBits) 
    {/* float32 to float64 */
      std::memcpy(&destination_mantissa, &source_mantissa, InSigBytes);
      destination_mantissa <<= (SigBits - InSigBits);
    } 
    else 
    {/* float32 to float32 */
      std::memcpy(&destination_mantissa, &source_mantissa, SigBytes);
    }
    
    const size_t ExpBits = Float::exponent_bits;
    const size_t ExpBytes = ExpBits/8 + (ExpBits%8!=0?1:0);
    const size_t InExpBits = InFloat::exponent_bits;
    const size_t InExpBytes = InExpBits/8 + (InExpBits%8!=0?1:0);
    const size_t ExpOffset = Float::exponent_offset;
    const size_t InExpOffset = InFloat::exponent_offset;
    const Bitset ExpBias = Float::exponent_bias_bitset();
    const InBitset InExpBias = InFloat::exponent_bias_bitset();
    // Copy exponent_bits
    InBitset source_exponent = (source & InExpMask);
    Bitset destination_exponent = 0;
    if (ExpBits < InExpBits) 
    {/* float64 to float32 */
      if (source_exponent >= InExpBias) 
      {
        source_exponent -= InExpBias;
        source_exponent >>= InExpOffset;
        /** 
         * (Inf - Bias)>>ExpOffset
         * = ((1<<ExpBits) - 1) - ((1<<(ExpBits-1)) - 1)
         * = (1<<ExpBits) - (1<<(ExpBits-1))
         * = 1 << (ExpBits-1)
        */
        const InBitset inf = InBitset(1)<<(ExpBits - 1);
        if (source_exponent >= inf) {
          destination = Inf|(((source & InSignMask) == InPositive) ? Positive : Negative);
          return;
        }
        std::memcpy(&destination_exponent, &source_exponent, ExpBytes);
        destination_exponent <<= ExpOffset;
        destination_exponent += ExpBias;
      } 
      else 
      {
        source_exponent = InExpBias - source_exponent;
        source_exponent >>= InExpOffset;
        const InBitset expbias = InBitset(1)<<(ExpBits - 1);
        if (source_exponent > expbias) {
          destination = Zero;
          return;
        }
        std::memcpy(&destination_exponent, &source_exponent, ExpBytes);
        destination_exponent <<= ExpOffset;
        destination_exponent = ExpBias - destination_exponent;
      }
    } 
    else if (ExpBits != InExpBits) 
    {/* float32 to float64 */
      if (source_exponent >= InExpBias) 
      {
        source_exponent -= InExpBias;
        source_exponent >>= InExpOffset;
        std::memcpy(&destination_exponent, &source_exponent, InExpBytes);
        destination_exponent <<= ExpOffset;
        destination_exponent += ExpBias;
      } 
      else 
      {
        source_exponent = InExpBias - source_exponent;
        source_exponent >>= InExpOffset;
        std::memcpy(&destination_exponent, &source_exponent, InExpBytes);
        destination_exponent <<= ExpOffset;
        destination_exponent = ExpBias - destination_exponent;
      }
    } 
    else 
    {/* float32 to float32 */
      source_exponent >>= InExpOffset;
      std::memcpy(&destination_exponent, &source_exponent, ExpBytes);
      destination_exponent <<= ExpOffset;
    }
    
    const Bitset SigMask = Float::significand_mask();
    const Bitset ExpMask = Float::exponent_mask();
    destination = (destination_mantissa&SigMask)|(destination_exponent&ExpMask)|(((source & InSignMask) == InPositive) ? Positive : Negative);
  }

public:
  Bitset _Mybitset;

  Bitset& bitset() {
    return _Mybitset;
  }

  const Bitset& bitset() const {
    return _Mybitset;
  }

  explicit FloatX(const Bitset& _bitset) : _Mybitset(_bitset) {}

public:
  FloatX() : _Mybitset() {}
    
  template<size_t s2, size_t e2, typename InBitset> 
  explicit FloatX(const FloatX<s2,e2,InBitset>& other) {
    do_cvtf<s2, e2>(other.bitset(), this->bitset());
  }
    
  FloatX(float other) {
    do_cvtf<23, 8>(reinterpret_cast<const std::bitset<32>&>(other), this->bitset());
  }
    
  FloatX(double other) {
    do_cvtf<52, 11>(reinterpret_cast<const std::bitset<64>&>(other), this->bitset());
  }

  FloatX(long double other) : FloatX(static_cast<double>(other)) {}

  FloatX(unsigned int other) {
    do_cvtui(std::bitset<sizeof(other)*8>(other), this->bitset());
  }

  FloatX(unsigned long long other) {
    do_cvtui(std::bitset<sizeof(other)*8>(other), this->bitset());
  }

  FloatX(bool other) {
    if( other ){
      _Mybitset = exponent_bias_bitset()/* | std::bitset<bits>(0)*/;
    } else {
      _Mybitset = zero_bitset();
    }
  }

  FloatX(int other) {
    if (other < 0) {
      do_cvtui(std::bitset<32>(-other), this->bitset());
      this->bitset() |= negsign_bitset();
    } else {
      do_cvtui(std::bitset<32>(other), this->bitset());
    }
  }

  FloatX(long long other) {
    if (other < 0) {
      do_cvtui(std::bitset<64>(-other), this->bitset());
      this->bitset() |= negsign_bitset();
    } else {
      do_cvtui(std::bitset<64>(other), this->bitset());
    }
  }
    
  explicit operator float() const {
    auto the_float = FloatX<23,8>(*this);
    return reinterpret_cast<const float&>(the_float);
  }

  explicit operator double() const {
    auto the_double = FloatX<52,11>(*this);
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
    Bitset exp2_mantissa_bits = Bitset(significand_bits) << exponent_offset;
    /*Bitset truncated_exponent = this_exponent - exp2_mantissa_bits;
    Bitset trunc_exponent = zero_exponent - truncated_exponent;*/
    /*Bitset trunc_exponent = (this_exponent>>exponent_offset) >= (exp2_mantissa_bits>>exponent_offset)
      ? zero_exponent - (this_exponent - exp2_mantissa_bits)
      : zero_exponent + (exp2_mantissa_bits - this_exponent);*/
    bool offset_exponent_sign = (this_exponent>>exponent_offset) >= (exp2_mantissa_bits>>exponent_offset);
    Bitset offset_exponent = offset_exponent_sign ? (this_exponent - exp2_mantissa_bits) : (exp2_mantissa_bits - this_exponent);
    
    Bitset this_significant = (_Mybitset & significand_mask()) | hidden_significand_mask();
    if ((!offset_exponent_sign) || (offset_exponent_sign && zero_exponent >= offset_exponent)) {
      Bitset trunc_exponent = offset_exponent_sign 
        ? zero_exponent - offset_exponent 
        : zero_exponent + offset_exponent;

      size_t trunc_bits = (trunc_exponent >> exponent_offset).to_ulong();
      size_t saved_bits = (significand_bits+1) - trunc_bits;
      if (saved_bits > sizeof(unsigned int) * 8) {
        throw std::overflow_error("floatX<...>::operator unsigned int() const");
      }

      // smaller-shift to correct integer-bits
      this_significant >>= ((significand_bits+1) - saved_bits);
    } else {
      size_t larger_shift = ((offset_exponent - zero_exponent) >> exponent_offset).to_ulong();
      this_significant <<= larger_shift;
    }

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
    Bitset exp2_mantissa_bits = Bitset(significand_bits) << exponent_offset;
    /*auto truncated_exponent = this_exponent - exp2_mantissa_bits;
    auto trunc_exponent = zero_exponent - truncated_exponent;*/
    Bitset trunc_exponent = (this_exponent>>exponent_offset) >= (exp2_mantissa_bits>>exponent_offset)
      ? zero_exponent - (this_exponent - exp2_mantissa_bits)
      : zero_exponent + (exp2_mantissa_bits - this_exponent);
    size_t trunc_bits = (trunc_exponent >> exponent_offset).to_ulong();
    size_t saved_bits = (significand_bits+1) - trunc_bits;
    if (saved_bits > sizeof(unsigned long long) * 8) {
      throw std::overflow_error("FloatX<...>::operator unsigned int() const");
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

  bool operator==(const FloatX& right) const {
    return do_equal(this->bitset(), right.bitset());
  }

  bool operator!=(const FloatX& right) const {
    return !do_equal(this->bitset(), right.bitset());
  }

  bool operator<(const FloatX& right) const {
    return do_less(this->bitset(), right.bitset());
  }

  bool operator>(const FloatX& right) const {
    return right < *this;
  }
    
  bool operator<=(const FloatX& right) const {
    return !(*this > right);
  }

  bool operator>=(const FloatX& right) const {
    return !(*this < right);
  }

  FloatX operator+() const {
    return *this;
  }

  FloatX operator-() const {
    return FloatX{ ((_Mybitset & sign_mask()) ^ sign_mask()) | (_Mybitset & (~sign_mask())) };
  }

  FloatX& operator+=(const FloatX& right) {
    /*std::bitset<bits> this_significant;
    std::bitset<bits> this_exponent;
    std::bitset<bits> this_sign;
    std::bitset<bits> right_significant;
    std::bitset<bits> right_exponent;
    std::bitset<bits> right_sign;
    std::bitset<bits> exponent_difference;
    std::bitset<bits> sticky_mask;
    floating_plus(_Mybitset, this_significant, this_exponent, this_sign,
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
      );*/
    do_plus(this->bitset(), right.bitset(), this->bitset());
    return *this;
  }

  FloatX& operator-=(const FloatX& right) {
    return *this += (-right);
  }

  FloatX& operator*=(const FloatX& right) {
    /*std::bitset<bits*2> this_significand;
    std::bitset<bits> this_exponent;
    std::bitset<bits*2> right_significand;
    std::bitset<bits> right_exponent;
    std::bitset<bits> exponent_difference;
    std::bitset<bits*2> roundup;
    floating_multiply(_Mybitset, this_significand, this_exponent,
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
      );*/
    do_multiply(this->bitset(), right.bitset(), this->bitset());
    return *this;
  }

  FloatX& operator/=(const FloatX& right) {
    /*std::bitset<bits> dividend;
    std::bitset<bits> this_significand;
    std::bitset<bits> this_exponent;
    std::bitset<bits> divisor;
    std::bitset<bits> right_exponent;
    std::bitset<bits> exponent_difference;
    std::bitset<bits> divbit;
    std::bitset<bits> roundup;
    floating_divide(_Mybitset, dividend, this_significand, this_exponent,
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
      );*/
    do_divide(this->bitset(), right.bitset(), this->bitset());
    return *this;
  }

  FloatX& operator%=(const FloatX& right) {
    abort();
  }
    
  FloatX operator+(const FloatX& right) const {
    return FloatX(*this) += right;
  }

  FloatX operator-(const FloatX& right) const {
    return FloatX(*this) -= right;
  }

  FloatX operator*(const FloatX& right) const {
    return FloatX(*this) *= right;
  }

  FloatX operator/(const FloatX& right) const {
    return FloatX(*this) /= right;
  }

  FloatX operator%(const FloatX& right) const {
    abort();
  }
};

/** TEMPLATE SPECIAL, IEEE754 single-precision */
template<typename Bitset>
class FloatX<23,8,Bitset,true> {
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
    
  FloatX() = default;
  FloatX(FloatX&&) = default;
  FloatX(const FloatX&) = default;
  template<size_t s2, size_t e2, typename InBitset, bool IsOpt>
  FloatX(const FloatX<s2,e2,InBitset,IsOpt>& other) : _Myfp(FloatX<23,8,Bitset,false>(other).as<float>()) {}
  FloatX& operator=(FloatX&&) = default;
  FloatX& operator=(const FloatX&) = default;
  template<size_t s2, size_t e2, typename InBitset, bool IsOpt>
  FloatX& operator=(const FloatX<s2,e2,InBitset,IsOpt>& other) { 
    _Myfp = FloatX<23,8,Bitset,false>(other).as<float>(); 
    return *this; 
  }

  constexpr FloatX(float other) : _Myfp(other) {}
  constexpr FloatX(double other) : FloatX(static_cast<float>(other)) {}
  constexpr FloatX(long double other) : FloatX(static_cast<float>(other)) {}
  constexpr FloatX(bool other) : FloatX(static_cast<float>(other)) {}
  constexpr FloatX(int other) : FloatX(static_cast<float>(other)) {}
  constexpr FloatX(long long other) : FloatX(static_cast<float>(other)) {}
  constexpr FloatX(unsigned int other) : FloatX(static_cast<float>(other)) {}
  constexpr FloatX(unsigned long long other) : FloatX(static_cast<float>(other)) {}
  constexpr operator float() const { return _Myfp; }
  constexpr explicit operator double() const { return static_cast<double>(this->operator float()); }
  constexpr explicit operator long double() const { return static_cast<long double>(this->operator float()); }
  constexpr explicit operator bool() const { return static_cast<bool>(this->operator float()); }
  constexpr explicit operator int() const { return static_cast<int>(this->operator float()); }
  constexpr explicit operator long long() const { return static_cast<long long>(this->operator float()); }
  constexpr explicit operator unsigned int() const { return static_cast<unsigned int>(this->operator float()); }
  constexpr explicit operator unsigned long long() const { return static_cast<unsigned long long>(this->operator float()); }

  inline constexpr bool operator==(FloatX right) const {
    return _Myfp == right._Myfp;
  }
  inline constexpr bool operator!=(FloatX right) const {
    return _Myfp != right._Myfp;
  }
  inline constexpr bool operator<(FloatX right) const {
    return _Myfp < right._Myfp;
  }
  inline constexpr bool operator>(FloatX right) const {
    return _Myfp > right._Myfp;
  }
  inline constexpr bool operator<=(FloatX right) const {
    return _Myfp <= right._Myfp;
  }
  inline constexpr bool operator>=(FloatX right) const {
    return _Myfp >= right._Myfp;
  }

  inline constexpr FloatX operator+() const {
    return FloatX{ _Myfp };
  }
  inline constexpr FloatX operator-() const {
    return FloatX{ -_Myfp };
  }
  inline constexpr FloatX operator+(FloatX right) const {
    return FloatX{ _Myfp + right._Myfp };
  }
  inline constexpr FloatX operator-(FloatX right) const {
    return FloatX{ _Myfp - right._Myfp };
  }
  inline constexpr FloatX operator*(FloatX right) const {
    return FloatX{ _Myfp * right._Myfp };
  }
  inline constexpr FloatX operator/(FloatX right) const {
    return FloatX{ _Myfp / right._Myfp };
  }
  inline FloatX operator%(FloatX right) const {
    return FloatX{ _Myfp - _CSTD floor(_Myfp / right._Myfp) * right._Myfp };
    //return _CSTD fmodf(_Myfp, right._Myfp);
  }

  inline FloatX& operator+=(FloatX right) {
    _Myfp += right._Myfp;
    return *this;
  }
  inline FloatX& operator-=(FloatX right) {
    _Myfp -= right._Myfp;
    return *this;
  }
  inline FloatX& operator*=(FloatX right) {
    _Myfp *= right._Myfp;
    return *this;
  }
  inline FloatX& operator/=(FloatX right) {
    _Myfp /= right._Myfp;
    return *this;
  }
  inline FloatX& operator%=(FloatX right) {
    _Myfp = _CSTD fmodf(_Myfp, right._Myfp);
    return *this;
  }
};

/** TEMPLATE SPECIAL, IEEE754 double-precision */
template<typename Bitset>
class FloatX<52,11,Bitset,true> {
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
    
  FloatX() = default;
  FloatX(FloatX&&) = default;
  FloatX(const FloatX&) = default;
  template<size_t s2, size_t e2, typename InBitset, bool IsOpt>
  FloatX(const FloatX<s2,e2,InBitset,IsOpt>& other) : _Myfp( FloatX<52,11,Bitset,false>(other).as<double>() ) {}
  FloatX& operator=(FloatX&&) = default;
  FloatX& operator=(const FloatX&) = default;
  template<size_t s2, size_t e2, typename InBitset, bool IsOpt>
  FloatX& operator=(const FloatX<s2,e2,InBitset,IsOpt>& other) {
    _Myfp = FloatX<52,11,Bitset,false>(other).as<double>();
    return *this;
  }

  constexpr FloatX(double other) : _Myfp(other) {}
  constexpr FloatX(float other) : FloatX(static_cast<double>(other)) {}
  constexpr FloatX(long double other) : FloatX(static_cast<double>(other)) {}
  constexpr FloatX(bool other) : FloatX(static_cast<double>(other)) {}
  constexpr FloatX(int other) : FloatX(static_cast<double>(other)) {}
  constexpr FloatX(long long other) : FloatX(static_cast<double>(other)) {}
  constexpr FloatX(unsigned int other) : FloatX(static_cast<double>(other)) {}
  constexpr FloatX(unsigned long long other) : FloatX(static_cast<double>(other)) {}
  constexpr operator double() const { return _Myfp; }
  constexpr explicit operator float() const { return  static_cast<float>(this->operator double()); }
  constexpr explicit operator long double() const { return static_cast<long double>(this->operator double()); }
  constexpr explicit operator bool() const { return  static_cast<bool>(this->operator double()); }
  constexpr explicit operator int() const { return static_cast<int>(this->operator double()); }
  constexpr explicit operator long long() const { return static_cast<long long>(this->operator double()); }
  constexpr explicit operator unsigned int() const { return static_cast<unsigned int>(this->operator double()); }
  constexpr explicit operator unsigned long long() const { return static_cast<unsigned long long>(this->operator double()); }

  inline constexpr bool operator==(FloatX right) const {
    return _Myfp == right._Myfp;
  }
  inline constexpr bool operator!=(FloatX right) const {
    return _Myfp != right._Myfp;
  }
  inline constexpr bool operator<(FloatX right) const {
    return _Myfp < right._Myfp;
  }
  inline constexpr bool operator>(FloatX right) const {
    return _Myfp > right._Myfp;
  }
  inline constexpr bool operator<=(FloatX right) const {
    return _Myfp <= right._Myfp;
  }
  inline constexpr bool operator>=(FloatX right) const {
    return _Myfp >= right._Myfp;
  }

  inline constexpr FloatX operator+() const {
    return FloatX{ _Myfp };
  }
  inline constexpr FloatX operator-() const {
    return FloatX( -_Myfp );
  }
  inline constexpr FloatX operator+(FloatX right) const {
    return FloatX( _Myfp + right._Myfp );
  }
  inline constexpr FloatX operator-(FloatX right) const {
    return FloatX( _Myfp - right._Myfp );
  }
  inline constexpr FloatX operator*(FloatX right) const {
    return FloatX( _Myfp * right._Myfp );
  }
  inline constexpr FloatX operator/(FloatX right) const {
    return FloatX( _Myfp / right._Myfp );
  }
  inline FloatX operator%(FloatX right) const {
    return _CSTD fmod(_Myfp, right._Myfp);
  }

  inline FloatX& operator+=(FloatX right) {
    _Myfp += right._Myfp;
    return *this;
  }
  inline FloatX& operator-=(FloatX right) {
    _Myfp -= right._Myfp;
    return *this;
  }
  inline FloatX& operator*=(FloatX right) {
    _Myfp *= right._Myfp;
    return *this;
  }
  inline FloatX& operator/=(FloatX right) {
    _Myfp /= right._Myfp;
    return *this;
  }
  inline FloatX& operator%=(FloatX right) {
    _Myfp = _CSTD fmod(_Myfp, right._Myfp);
    return *this;
  }
};

// next_work: { XXXX quadruple-precision-template spectialation }

#define __calculation_float_operator_with_literal(_OP_, _LITERAL_TYPE_) \
template<size_t s, size_t e, typename b, bool opt> inline \
FloatX<s,e,b,opt> operator##_OP_##(const FloatX<s,e,b,opt>& left, _LITERAL_TYPE_ right) { \
  return left _OP_ static_cast< FloatX<s,e,b,opt> >(right); \
}

#define __calculation_float_operator_with_literal2(_OP_, _LITERAL_TYPE_) \
template<size_t s, size_t e, typename b, bool opt> inline \
FloatX<s,e,b,opt> operator##_OP_##(_LITERAL_TYPE_ left, const FloatX<s,e,b,opt>& right) { \
  return static_cast< FloatX<s,e,b,opt> >(left) _OP_ right; \
}

#define __calculation_float_operator_with_literal_commutatibity(_OP_, _LITERAL_TYPE_) \
template<size_t s, size_t e, typename b, bool opt> inline \
FloatX<s,e,b,opt> operator##_OP_##(const FloatX<s,e,b,opt>& left, _LITERAL_TYPE_ right) { \
  return left _OP_ static_cast< FloatX<s,e,b,opt> >(right); \
} \
template<size_t s, size_t e, typename b, bool opt> inline \
FloatX<s,e,b,opt> operator##_OP_##(_LITERAL_TYPE_ left, const FloatX<s,e,b,opt>& right) { \
  return static_cast< FloatX<s,e,b,opt> >(left) _OP_ right; \
}

#define __calculation_float_lvalueoperator_with_literal(_OP_, _LITERAL_TYPE_) \
template<size_t s, size_t e, typename b, bool opt> inline \
FloatX<s,e,b,opt>& operator##_OP_##(FloatX<s,e,b,opt>& left, _LITERAL_TYPE_ right) { \
  return left _OP_ static_cast< FloatX<s,e,b,opt> >(right); \
}

#define __calculation_float_comparison_with_literal_commutatibity(_OP_, _LITERAL_TYPE_) \
template<size_t s, size_t e, typename b, bool opt> inline \
bool operator##_OP_##(const FloatX<s,e,b,opt>& left, _LITERAL_TYPE_ right) { \
  return left _OP_ static_cast< FloatX<s,e,b,opt> >(right); \
} \
template<size_t s, size_t e, typename b, bool opt> inline \
bool operator##_OP_##(_LITERAL_TYPE_ left, const FloatX<s,e,b,opt>& right) { \
  return static_cast< FloatX<s,e,b,opt> >(left) _OP_ right; \
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


using std::isinf;
using std::isnan;
using std::abs;
using std::trunc;
using std::modf;
using std::frexp;
using std::ldexp;
using std::floor;
using std::ceil;
using std::round;
using std::fmod;

template<size_t m, size_t e, typename b> inline
bool isinf(const FloatX<m,e,b>& x) {
  return FloatX<m,e,b>::do_isinf(x.bitset());
}  

template<size_t m, size_t e, typename b> inline
bool isnan(const FloatX<m,e,b>& x) {
  return FloatX<m,e,b>::do_isnan(x.bitset());
}

template<size_t m, size_t e, typename b> inline
FloatX<m,e,b> abs(const FloatX<m,e,b>& x) {
  auto abs_mask = ~(FloatX<m,e,b>::sign_mask());
  return FloatX<m,e,b>{ x.bitset() & abs_mask };
}

template<size_t s, size_t e, typename b> inline
FloatX<s,e,b> trunc(const FloatX<s,e,b>& x) {
  FloatX<s,e,b> ipart;
  FloatX<s,e,b> fpart;
  FloatX<s,e,b>::do_modf(x.bitset(), fpart.bitset(), ipart.bitset());
  return ipart;
}

template<size_t s, size_t e, typename b> inline
FloatX<s,e,b> modf(const FloatX<s,e,b>& x, FloatX<s,e,b>* ipart) {
  FloatX<s,e,b> fpart;
  FloatX<s,e,b>::do_modf(x.bitset(), fpart.bitset(), ipart->bitset());
  return fpart;
}

template<size_t s, size_t e, typename b, typename Integer> inline
FloatX<s,e,b> frexp(const FloatX<s,e,b>& x, Integer* exp2) {
  FloatX<s,e,b> base;
  FloatX<s,e,b>::do_frexp(x.bitset(), base.bitset(), *exp2);
  return base;
}

template<size_t s, size_t e, typename b, typename Integer> inline
FloatX<s,e,b> ldexp(const FloatX<s,e,b>& base, const Integer& exp2) {
  FloatX<s,e,b> y;
  FloatX<s,e,b>::do_ldexp(base.bitset(), exp2, y.bitset());
  return y;
}

template<size_t m, size_t e, typename b> inline
FloatX<m,e,b> floor(const FloatX<m,e,b>& x) {
  if ( x >= 0 ) {
    return trunc(x);
  } else {
    return trunc(x) - 1;
  }
}

template<size_t m, size_t e, typename b> inline
FloatX<m,e,b> ceil(const FloatX<m,e,b>& x) {
  return floor(x) + 1;
}

template<size_t m, size_t e, typename b> inline
FloatX<m,e,b> round(const FloatX<m,e,b>& x) {
  return floor(x + 0.5);
}

template<size_t s, size_t e, typename b> inline
FloatX<s,e,b> fmod(const FloatX<s,e,b>& x, const FloatX<s,e,b>& y) {
  return x - y * trunc(x/y);
}

using std::sqrt;
using std::cbrt;
using std::exp;
using std::log;
using std::pow;
using std::sin;
using std::cos;
using std::tan;
using std::asin;
using std::acos;
using std::atan;
using std::atan2;

/** @see 'function.h' */


/** IEEE754 single-precision */
using Float32 = FloatX<23, 8, std::bitset<32>, true>;
/** IEEE754 double-precision */
using Float64 = FloatX<52, 11, std::bitset<64>, true>;
/** @see GCC/quadmath/ALL */
using Float128 = FloatX<112, 15>;
using Float256 = FloatX<235, 20>;
using Float512 = FloatX<485, 26>;
using Float1024 = FloatX<992, 31>;

inline bool isinf(Float32 x) { return std::isinf(x); }
inline bool isinf(Float64 x) { return std::isinf(x); }
inline bool isnan(Float32 x) { return std::isnan(x); }
inline bool isnan(Float64 x) { return std::isnan(x); }
inline Float32 abs(Float32 x) { return std::fabsf(x); }
inline Float64 abs(Float64 x) { return std::fabs(x); }
inline Float32 trunc(Float32 x) { return std::truncf(x); }
inline Float64 trunc(Float64 x) { return std::trunc(x); }
inline Float32 modf(Float32 x, Float32* i) { return std::modf(x,(float*)i); }
inline Float64 modf(Float64 x, Float64* i) { return std::modf(x,(double*)i); }
inline Float32 frexp(Float32 x, int* ep) { return std::frexp(x, ep); }
inline Float64 frexp(Float64 x, int* ep) { return std::frexp(x, ep); }
inline Float32 ldexp(Float32 x, int ep) { return std::ldexp(x, ep); }
inline Float64 ldexp(Float64 x, int ep) { return std::ldexp(x, ep); }
inline Float32 floor(Float32 x) { return std::floorf(x); }
inline Float64 floor(Float64 x) { return std::floor(x); }
inline Float32 ceil(Float32 x) { return std::ceilf(x); }
inline Float64 ceil(Float64 x) { return std::ceil(x); }
inline Float32 round(Float32 x) { return std::roundf(x); }
inline Float64 round(Float64 x) { return std::round(x); }
inline Float32 fmod(Float32 x, Float32 y) { return std::fmodf(x, y); }
inline Float64 fmod(Float64 x, Float64 y) { return std::fmod(x, y); }

inline Float32 sqrt(Float32 x) { return std::sqrtf(x); }
inline Float64 sqrt(Float64 x) { return std::sqrt(x); }
inline Float32 cbrt(Float32 x) { return std::cbrtf(x); }
inline Float64 cbrt(Float64 x) { return std::cbrt(x); }
inline Float32 exp(Float32 x) { return std::expf(x); }
inline Float64 exp(Float64 x) { return std::exp(x); }
inline Float32 log(Float32 x) { return std::logf(x); }
inline Float64 log(Float64 x) { return std::log(x); }
inline Float32 pow(Float32 x, Float32 power) { return std::powf(x, power); }
inline Float64 pow(Float64 x, Float64 power) { return std::pow(x, power); }
inline Float32 sin(Float32 x) { return std::sinf(x); }
inline Float64 sin(Float64 x) { return std::sin(x); }
inline Float32 cos(Float32 x) { return std::cosf(x); }
inline Float64 cos(Float64 x) { return std::cos(x); }
inline Float32 tan(Float32 x) { return std::tanf(x); }
inline Float64 tan(Float64 x) { return std::tan(x); }
inline Float32 asin(Float32 x) { return std::asinf(x); }
inline Float64 asin(Float64 x) { return std::asin(x); }
inline Float32 acos(Float32 x) { return std::acosf(x); }
inline Float64 acos(Float64 x) { return std::acos(x); }
inline Float32 atan(Float32 x) { return std::atanf(x); }
inline Float64 atan(Float64 x) { return std::atan(x); }
inline Float32 atan2(Float32 y, Float32 x) { return std::atan2f(y,x); }
inline Float64 atan2(Float64 y, Float64 x) { return std::atan2(y,x); }
}// end of namespace calculation


#include <cassert>
#include <deque>
#include <algorithm>
#include <xstring>
namespace calculation {
// { used for FloatX<X,X> to decimal to decimal_string }
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
template<size_t m, size_t e, typename b, bool isopt>
std::string to_string(FloatX<m,e,b,isopt> _Source) {
  using floatX_t = FloatX<m, e, b, isopt>;

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
  for (size_t i = 0; i != _Source.significand_bits; ++i) {
    _Dest.divev2();
  }

  // correct exponent
  auto _Zero = std::bitset<floatX_t::bits>(0);
  auto _One = std::bitset<floatX_t::bits>(1) << _Source.exponent_offset;
  auto _Exponent = _Source.bitset() & _Source.exponent_mask();
  //_Exponent -= std::bitset<floatX_t::bits>(_Source.significand_bits) << _Source.exponent_offset;
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

template<size_t m, size_t e, typename Bitset, bool IsOpt> inline
std::ostream& operator<<(std::ostream& _Ostr, FloatX<m,e,Bitset,IsOpt> _Fp) {
  return (_Ostr << to_string(_Fp));
}
}// end of namespace calculation


#if defined __has_include && __has_include(<limits>) && defined __has_include && __has_include(<type_traits>)
#include <limits>
_STD_BEGIN
template<size_t m, size_t e, typename b, bool opt>
class numeric_limits< calculation::FloatX<m,e,b,opt> > : public _Num_float_base{
public:
  _NODISCARD static constexpr calculation::FloatX<m,e,b,opt>(min)() noexcept {
    abort();
  }

  _NODISCARD static constexpr calculation::FloatX<m,e,b,opt>(max)() noexcept {
    abort();
  }

  _NODISCARD static constexpr calculation::FloatX<m,e,b,opt> lowest() noexcept {
    abort();
  }

  _NODISCARD static calculation::FloatX<m,e,b,opt> epsilon() noexcept {
    return calculation::FloatX<m,e,b,opt>( calculation::FloatX<m,e,b,opt>::epsilon_bitset() );
  }

  _NODISCARD static calculation::FloatX<m,e,b,opt> round_error() noexcept {
    return epsilon() / 2;
  }

  _NODISCARD static constexpr calculation::FloatX<m,e,b,opt> denorm_min() noexcept {
    abort();
  }

  _NODISCARD static calculation::FloatX<m,e,b,opt> infinity() noexcept {
    return calculation::FloatX<m,e,b,opt>( calculation::FloatX<m,e,b,opt>::inf_bitset() );
  }

  _NODISCARD static calculation::FloatX<m,e,b,opt> quiet_NaN() noexcept {
    return calculation::FloatX<m,e,b,opt>( calculation::FloatX<m,e,b,opt>::quiet_NaN_bitset() );
  }

  _NODISCARD static calculation::FloatX<m,e,b,opt> signaling_NaN() noexcept {
    return calculation::FloatX<m,e,b,opt>( calculation::FloatX<m,e,b,opt>::signaling_NaN_bitset() );
  }

  static constexpr size_t digits = calculation::FloatX<m,e,b,opt>::significand_bits + 1;
};
template<size_t s, size_t e, typename b, bool opt> constexpr 
  bool is_floating_point_v<calculation::FloatX<s,e,b,opt>> = true;
_STD_END
#endif


//template<typename Bitset>
//void floating_plus(
//  /*INOUT*/Bitset& this_floating,
//  /*Reg1*/ Bitset& this_significant, 
//  /*Reg2*/ Bitset& this_exponent,
//  /*Reg3*/ Bitset& this_sign,
//    const  Bitset& right_floating, 
//  /*Reg4*/ Bitset& right_significant,
//  /*Reg5*/ Bitset& right_exponent,
//  /*Reg6*/ Bitset& right_sign,
//  /*Reg7*/ Bitset& exponent_difference,
//  /*Reg8*/ Bitset& sticky_mask,
//  const Bitset Zero,
//  const Bitset Inf,
//  const Bitset NaN, 
//  const Bitset SigMask,
//  const Bitset hidden_significand_mask,
//  const Bitset exponent_mask, 
//  const Bitset sign_mask,
//  const size_t exponent_offset,
//  const size_t hidden_significant_offset,
//  const size_t lastbit_offset) 
//{
//  /** 0. check special case 'zero' */
//  if ( (this_floating & (~sign_mask)) == Zero ) {
//    this_floating = right_floating;
//    return;
//  }
//  if ( (right_floating & (~sign_mask)) == Zero ) {
//    return;
//  }
//
//  /** 0. check special case 'inf' and 'nan' */
//  if ( (this_floating & Inf) == Inf ) {
//    this_floating = NaN;
//    return;
//  }
//  if ( (right_floating & Inf) == Inf ) {
//    this_floating = NaN;
//    return;
//  }
//
//  /** 1. decode into significand-bitset, exponent-bitset, sign-bitset */
//  this_significant = (this_floating&SigMask)|hidden_significand_mask;
//  this_exponent = this_floating&exponent_mask;
//  this_sign = this_floating&sign_mask;
//  right_significant = (right_floating&SigMask)|hidden_significand_mask;
//  right_exponent = right_floating&exponent_mask;
//  right_sign = right_floating&sign_mask;
//  
//  /** 2. extend more accuracy */
//  assert( hidden_significant_offset <= lastbit_offset );
//  size_t shift = 0;
//  if ( hidden_significant_offset < lastbit_offset ) {
//    /**
//     *    11111111
//     * +  11111111
//     * = 111111110
//     * 
//     * so we only reserved '1' bit
//    */
//    const size_t reserved_bit = 1;
//    shift = lastbit_offset - hidden_significant_offset;
//    shift -= reserved_bit;
//    this_significant <<= shift;
//    right_significant <<= shift;
//  }
//
//  /** 3. sync exponent-bitset, to greater */
//  if ( this_exponent < right_exponent ) {
//
//    exponent_difference = right_exponent - this_exponent;
//    size_t smaller_shift = static_cast<size_t>( (exponent_difference >> exponent_offset).to_ullong() );
//    
//    if ( smaller_shift - 1 <= lastbit_offset && need_roundup_for_nearest_even(this_significant, sticky_mask, smaller_shift) ){
//      this_significant >>= smaller_shift;
//      this_significant += 1;
//    } else {
//      this_significant >>= smaller_shift;
//    }
//    
//    // assert( inf_bitset - exponent_difference < this_exponent );
//    this_exponent += exponent_difference;
//  
//  } 
//  else if ( this_exponent != right_exponent ) {
//
//    exponent_difference = this_exponent - right_exponent;
//    size_t smaller_shift = static_cast<size_t>( (exponent_difference >> exponent_offset).to_ullong() );
//
//    if ( smaller_shift - 1 <= lastbit_offset && need_roundup_for_nearest_even(right_significant, sticky_mask, smaller_shift) ){
//      right_significant >>= smaller_shift;
//      right_significant += 1;
//    } else {
//      right_significant >>= smaller_shift;
//    }
//
//    // assert( inf_bitset - exponent_difference < right_exponent );
//    right_exponent += exponent_difference;
//
//  }
//
//  /** 4. add significand-bitset, update sign-bits */
//  if ( this_sign == right_sign ) {
//    this_significant += right_significant;
//  } 
//  else {
//    if (this_significant == right_significant) {
//      this_floating = Zero;
//      return;
//    } else if (this_significant < right_significant) {
//      this_significant = right_significant - this_significant;
//      this_sign = right_sign;
//    } else  /* this_significand > right_significand */{
//      this_significant -= right_significant;
//      //this_sign = this_sign;
//    }
//  }
//
//  /** 5. normalize significand-bits and exponent-bits */
//  if ( (this_significant&(~SigMask)) != 0 ) {
//
//    // normalize for select 'shift'
//    size_t smaller_shift = 0;
//    while ( !this_significant.test(lastbit_offset - smaller_shift) ) {
//      ++smaller_shift;
//    }
//    smaller_shift = lastbit_offset - hidden_significant_offset - smaller_shift;
//
//    // normalize for 'significand' with 'round-up'
//    if ( need_roundup_for_nearest_even(this_significant, sticky_mask, smaller_shift) ){
//      this_significant >>= smaller_shift;
//      this_significant += 1;
//      if ( this_significant.test(hidden_significant_offset + 1) ) {
//        this_significant >>= 1;
//        smaller_shift += 1;
//        /**
//         *     1111 1111
//         * +   0000 0001
//         * = 1 0000 0000
//         *            |
//         *         this is test-bit
//         *            |
//         *         it must be zero, if carry
//        */
//      }
//    } else {
//      this_significant >>= smaller_shift;
//    }
//
//    // normalize for 'exponent'
//    if ( smaller_shift > shift ) {
//      exponent_difference = smaller_shift - shift;
//      exponent_difference <<= exponent_offset;
//      if ( Inf - exponent_difference <= this_exponent ) {
//        /** special case 'infinite large' */
//        this_floating = Inf|this_sign;
//        return;
//      }
//      this_exponent += exponent_difference;
//    } else {
//      exponent_difference = shift - smaller_shift;
//      exponent_difference <<= exponent_offset;
//      if ( this_exponent < exponent_difference ) {
//        /** special case 'infinite small' */
//        this_floating = Zero;
//        return;
//      }
//      this_exponent -= exponent_difference;
//    }
//    
//  } 
//  else {
//
//    // normalize for select 'shift'
//    size_t larger_shift = 0;
//    while ( ! this_significant.test(hidden_significant_offset - larger_shift) ) {
//      ++larger_shift;
//    }
//
//    // normalize for 'significand'
//    this_significant <<= larger_shift;
//
//    // normalize for 'exponent'
//    exponent_difference = shift + larger_shift;
//    exponent_difference <<= exponent_offset;
//    if ( this_exponent < exponent_difference ) {
//      /** special case 'infinite small' */
//      this_floating = Zero;
//      return;
//    }
//    this_exponent -= exponent_difference;
//
//  }
//
//  /** 6. encode ... */
//  this_floating = this_sign|this_exponent|(this_significant&SigMask);
//}

//template<typename Bitset, typename BigBitset>
//void floating_multiply(
//  /*INOUT*/Bitset&    this_floating,
//  /*Reg1*/ BigBitset& this_significand,
//  /*Reg2*/ Bitset&    this_exponent,
//    const  Bitset&    right_floating,
//  /*Reg3*/ BigBitset& right_significand,
//  /*Reg4*/ Bitset&    right_exponent,
//  /*Reg5*/ Bitset&    exponent_difference,
//  /*Reg6*/ BigBitset& sticky_mask,
//  const Bitset zero_bitset,
//  const Bitset inf_bitset,
//  const Bitset nan_bitset,
//  const Bitset significand_mask,
//  const Bitset hidden_significand_mask,
//  const Bitset exponent_mask,
//  const Bitset sign_mask,
//  const Bitset exponent_bias_bitset,
//  const size_t significand_bits,
//  const size_t bits)
//{
//  /** 0. check special case 'zero' */
//  if ( (this_floating & (~sign_mask)) == zero_bitset || (right_floating & (~sign_mask)) == zero_bitset ) {
//    this_floating = zero_bitset;
//    return;
//  }
//
//  /** 0. check special case 'inf' and 'nan' */
//  if ( (this_floating&inf_bitset) == inf_bitset || (right_floating&inf_bitset) == inf_bitset ) {
//    if( (this_floating&significand_mask) != 0 || (right_floating&significand_mask) != 0 ) {
//      this_floating = nan_bitset;
//      return;
//    }
//
//    this_floating = inf_bitset | ((this_floating&sign_mask)^(right_floating&sign_mask));
//    return;
//  }
//
//  /** 1. decode ... */
//  this_significand = 0;
//  reinterpret_cast<Bitset&>(this_significand) = (this_floating&significand_mask)|hidden_significand_mask;
//  this_exponent = this_floating&exponent_mask;
//  right_significand = 0;
//  reinterpret_cast<Bitset&>(right_significand) = (right_floating&significand_mask)|hidden_significand_mask;
//  right_exponent = right_floating&exponent_mask;
//
//  /** 2. significand multiplication 
//   * this_significand<<(significand_bits*2) = (this_significand<<significand_bits * right_significand<<significand_bits)
//   * this_significand<<significand_bits = (this_significand<<significand_bits * right_significand<<significand_bits) >> significand_bits
//   * this_significand<<'shift'          = (this_significand<<significand_bits * right_significand<<significand_bits) >> significand_bits
//  */
//  this_significand *= right_significand;
//  size_t shift = significand_bits;
//  
//  /** 3. this_exponent = this_exponent + right_exponent */
//  if ( right_exponent > exponent_bias_bitset ) {  
//    this_exponent += (right_exponent - exponent_bias_bitset);
//  } else if ( right_exponent < exponent_bias_bitset ) {
//    this_exponent -= (exponent_bias_bitset - right_exponent);
//  }
//  
//  /** 4. normalize, assert( this_significand&(~significand_mask) != 0 ); */
//  {
//    // normalize for select 'shift'
//    const size_t lastbit_offset = bits * 2 - 1;
//    const size_t hidden_significand_offset = significand_bits;
//    size_t smaller_shift = 0;
//    while ( !this_significand.test(lastbit_offset - smaller_shift) ) {
//      ++smaller_shift;
//    }
//    smaller_shift = lastbit_offset - hidden_significand_offset - smaller_shift;
//
//    // normalize for 'significand' with 'round-up'
//    if ( need_roundup_for_nearest_even(this_significand, sticky_mask, smaller_shift) ) {
//      this_significand >>= smaller_shift;
//      this_significand += 1;
//      if ( this_significand.test(hidden_significand_offset + 1) ) {
//        this_significand >>= 1;
//        smaller_shift += 1;
//        // assert( round-bit == 0 );
//      }
//    } else {
//      this_significand >>= smaller_shift;
//    }
//
//    // normalize for 'exponent'
//    size_t exponent_offset = significand_bits;
//    if ( smaller_shift > shift ) {
//      exponent_difference = (smaller_shift - shift);
//      exponent_difference <<= exponent_offset;
//      if ( inf_bitset - exponent_difference <= this_exponent ) {
//        /** special case 'infinite large' */
//        this_floating = inf_bitset|((this_floating&sign_mask)^(right_floating&sign_mask));
//        return;
//      }
//      this_exponent += exponent_difference;
//    } else {
//      exponent_difference = (shift - smaller_shift);
//      exponent_difference <<= exponent_offset;
//      if ( this_exponent < exponent_difference ) {
//        /** special case 'infinite small' */
//        this_floating = zero_bitset;
//        return;
//      }
//      this_exponent -= exponent_difference;
//    }
//  }
//
//  /** 4. encode ... */
//  this_floating = 
//    (reinterpret_cast<const Bitset&>(this_significand)&significand_mask)
//    | this_exponent
//    | ( (this_floating&sign_mask)^(right_floating&sign_mask) );
//}

//template<typename Bitset>
//void floating_divide(
//  /*INOUT*/Bitset& this_floating,
//  /*Reg1*/ Bitset& dividend,
//  /*Reg2*/ Bitset& this_significand,
//  /*Reg3*/ Bitset& this_exponent,
//    const  Bitset& right_floating, 
//  /*Reg4*/ Bitset& divisor,
//  /*Reg5*/ Bitset& right_exponent,
//  /*Reg6*/ Bitset& exponent_difference,
//  /*Reg7*/ Bitset& divbit,
//  /*Reg8*/ Bitset& sticky_mask,
//  const Bitset zero_bitset,
//  const Bitset inf_bitset,
//  const Bitset nan_bitset,
//  const Bitset significand_mask,
//  const Bitset hidden_significand_mask,
//  const Bitset exponent_mask,
//  const Bitset sign_mask,
//  const Bitset exponent_bias_bitset,
//  const size_t exponent_offset,
//  const size_t hidden_significant_offset,
//  const size_t lastbit_offset)
//{
//  /** 0. check special case 'zero' */
//  if ( (this_floating & (~sign_mask)) == zero_bitset && (right_floating & (~sign_mask)) == zero_bitset ) {
//    this_floating = nan_bitset;
//    return;
//  } else if ( (this_floating & (~sign_mask)) == zero_bitset ) {
//    // not change
//    return;
//  } else if ( (right_floating & (~sign_mask)) == zero_bitset ) {
//    this_floating = inf_bitset | ((this_floating&sign_mask)^(right_floating&sign_mask));
//    return;
//  }
//
//  /** 0. check special case 'inf' and 'nan' */
//  if ( (this_floating & inf_bitset) == inf_bitset || (right_floating & inf_bitset) == inf_bitset ) {
//    this_floating = nan_bitset;
//    return;
//  }
//
//  /** 1. decode into divide-operator, exponent-bitset */
//  dividend = (this_floating&significand_mask)|hidden_significand_mask;
//  this_exponent = this_floating&exponent_mask;
//  divisor = (right_floating&significand_mask)|hidden_significand_mask;
//  right_exponent = right_floating&exponent_mask;
//  
//  /** 2. extend more accuracy, significand-largest-bit
//   * [1<<(significand_bits-1), 1<<significand_bits] 
//   * to 
//   * [1<<(significand_bits+1), 1<<(significand_bits+2)] */
//  assert( hidden_significant_offset <= lastbit_offset );
//  size_t shift = 0;
//  if ( hidden_significant_offset < lastbit_offset ) {
//    /**
//     *    11111111
//     * /  10000000
//     * ------------
//     * =  11111111
//     * 
//     * so we only reserved '0' + '1'(division) bit
//    */
//    const size_t reserved_bit = 1;
//    shift = lastbit_offset - hidden_significant_offset;
//    shift -= reserved_bit;
//    dividend <<= shift;
//    divisor <<= shift;
//  }
//
//  /** 3. significand division */
//  this_significand = 0;
//  size_t offset = hidden_significant_offset + shift;// contains hidden-significand
//  do {
//    if ( dividend >= divisor ) {
//      divbit = 1;  divbit <<= offset;
//      this_significand |= divbit;
//      dividend -= divisor;
//    }
//
//    dividend <<= 1;
//
//  } while ( --offset != size_t(-1) );
//
//  /** 4. this_exponent = this_exponent - right_exponent */
//  if (right_exponent > exponent_bias_bitset) { 
//    this_exponent -= (right_exponent - exponent_bias_bitset);
//  } else if (right_exponent < exponent_bias_bitset) { 
//    this_exponent += (exponent_bias_bitset - right_exponent);
//  }
//
//  /** 5. normalize, assert( this_significand&(~significand_mask) != 0 ); see(step 2.) */
//  {
//    // normalize for select 'shift'
//    size_t smaller_shift = 0;
//    while ( !this_significand.test(lastbit_offset - smaller_shift) ) {
//      ++smaller_shift;
//    }
//    smaller_shift = lastbit_offset - hidden_significant_offset - smaller_shift;
//
//    // normalize for 'significand' with 'round-up'
//    if ( need_roundup_for_nearest_even(this_significand, smaller_shift, sticky_mask) ) {
//      this_significand >>= smaller_shift;
//      this_significand += 1;
//      if ( this_significand.test(hidden_significant_offset + 1) ) {
//        this_significand >>= 1;
//        smaller_shift += 1;
//        // assert( round-bit == 0 );
//      }
//    } else {
//      this_significand >>= smaller_shift;
//    }
//
//    // normalize for 'exponent'
//    assert( smaller_shift <= shift );
//    exponent_difference = shift - smaller_shift;
//    exponent_difference <<= exponent_offset;
//    if ( inf_bitset - exponent_difference <= this_exponent ) {
//      /** special case 'infinite small' */
//      this_floating = zero_bitset;
//      return;
//    }
//    this_exponent -= exponent_difference;
//  }
//
//  /** 6. encode ... */
//  this_floating = 
//    (this_significand&significand_mask)
//    | this_exponent
//    | ( (this_floating&sign_mask)^(right_floating&sign_mask) );
//}


//#define floating_operation_special_case(PROC_ZERO, PROC_INF, PROC_NAN) \
//if ( (floating & ~sign_mask) == zero_bitset ) {       \
//  /** special case 'zero' */                         \
//  PROC_ZERO;                                         \
//}                                                    \
//if ( (floating & inf_bitset) == inf_bitset ) {       \
//  if ( (floating & significand_mask) != 0 ) {        \
//    /** special case 'NaN' */                        \
//    PROC_NAN;                                        \
//  }                                                  \
//                                                    \
//  /** special case 'inf' */                       \
//  PROC_INF;                              \
//}                               \
//
//
//
//template<typename Bitset> 
//void floating_significand(
//  const Bitset& floating, 
//        Bitset& significand,
//  const Bitset zero_bitset,
//  const Bitset inf_bitset,
//  const Bitset nan_bitset,
//  const Bitset significand_mask,
//  const Bitset sign_mask,
//  const Bitset exponent_bias_bitset,
//  const Bitset possign_bitset)
//{
//  floating_operation_special_case(
//    significand = zero_bitset;
//    return,
//    significand = nan_bitset;
//    return,
//    significand = nan_bitset;
//    return  
//    );
//
//  significand = floating & significand_mask;
//  significand |= exponent_bias_bitset;
//  significand |= possign_bitset;
//}
///**
// * We want to cover static_floating and dynamic_floating by function once.
// * But this idea should not work, in dynamic_floating operation, two floating has different floating_constants and different bitwise_op.
// * Therefore, we only write the static_floating opeartions in a class, because other methods are meaningless.
//*/
//template<typename Bitset>
//void floating_exponent(
//   const  Bitset& floating, 
//  /*OUT*/ Bitset& exponent, 
//  /*Reg1*/Bitset& exponent_part, 
//  /*Reg2*/Bitset& sign_part,
//  const Bitset zero_bitset,
//  const Bitset inf_bitset,
//  const Bitset nan_bitset,
//  const Bitset significand_mask,
//  const Bitset exponent_mask, 
//  const Bitset sign_mask,
//  const Bitset exponent_bias_bitset, 
//  const size_t last_bit_offset,
//  const size_t hidden_significant_offset,
//  const size_t exponent_offset,
//  const Bitset possign_bitset,
//  const Bitset negsign_bitset)
//{
//  floating_operation_special_case(
//    exponent = zero_bitset;
//    return,
//    exponent = inf_bitset;
//    return,
//    exponent = nan_bitset;
//    return  
//    );
//
//  // get unnormalized-exponent and exponent-sign
//  exponent = floating & exponent_mask;
//  if ( exponent == exponent_bias_bitset ) {
//    exponent = zero_bitset;
//    return;
//  } else if ( exponent > exponent_bias_bitset ) {
//    exponent -= exponent_bias_bitset;
//    sign_part = possign_bitset;
//  } else   /* exponent < exponent_bias_bitset */{
//    exponent = exponent_bias_bitset - exponent;
//    sign_part = negsign_bitset;
//  }
//        
//  // normalize select shift
//  size_t smaller_shift = 0;
//  while ( !exponent.test(last_bit_offset - smaller_shift) ) {
//    ++smaller_shift;
//  }
//  smaller_shift = last_bit_offset - hidden_significant_offset - smaller_shift;
//  // normalize significand
//  exponent >>= smaller_shift;
//  exponent &= significand_mask;
//  // normalize exponent
//  exponent_part = smaller_shift;
//  exponent_part <<= exponent_offset;
//  exponent_part += exponent_bias_bitset;
//
//  // combine
//  exponent |= exponent_part;
//  exponent |= sign_part;
//}
//    
//template<typename Bitset>
//void floating_sign(
//  const Bitset& floating, 
//        Bitset& sign,
//  const Bitset zero_bitset,
//  const Bitset exponent_bias_bitset,
//  const Bitset sign_mask)
//{
//  sign = zero_bitset;
//  sign |= exponent_bias_bitset;
//  sign |= (floating & sign_mask);
//}