/**
 * @license
 *   Please identify Author, 2019 - 2021
 * @author 
 *   LongJiangnan, Jiang1998Nan@outlook.com
 * @brief
 *   infinite precision calculation for rational number
 * @readme 
 *   next work is dynamic infinite precision rationalnumber
*/
#pragma once

#include <numeric>// std::gcd
#include <math.h>// _CSTD abs
#include <iosfwd>// std::basic_ostream
#include <exception>// std::underflow_error, std::overflow_error ...
namespace calculation 
  {
using _CSTD abs;

template<typename Integer>
Integer gcd(Integer _Ax, Integer _Bx) {
  Integer a = abs(_Ax);
  Integer b = abs(_Bx);
  if (b > a) {
    std::swap(a, b);
  }

  while (b != 0) {
    Integer tmp = b;
    b = a % b;
    a = tmp;
  }

  return a != 0 ? a : 1;
}

template<typename IntegerX, bool is_signed_X = std::is_signed_v<IntegerX>>
struct RationalX {
  using Rational = RationalX<IntegerX, is_signed_X>;
  using Integer = IntegerX;

  Integer numerator;
  Integer denominator;

  constexpr RationalX() : numerator(static_cast<Integer>(0)), denominator(static_cast<Integer>(1)) {}

  constexpr RationalX(Integer num, Integer den = static_cast<Integer>(1)) : numerator(num), denominator(den) {
    Integer common_divisor = gcd(this->numerator, this->denominator);
    this->numerator /= common_divisor;
    this->denominator /= common_divisor;
  }

  RationalX(int lll) {
    this->numerator = static_cast<Integer>(lll);
    this->denominator = static_cast<Integer>(1);
  }

  RationalX(size_t lll) {
    this->numerator = static_cast<Integer>(lll);
    this->denominator = static_cast<Integer>(1);
  }

  bool operator==(const Rational& right) const {
    // assert( simplest(this) && simplest(right) )
    return this->numerator == right.numerator
      && right.denominator == this->denominator;
  }

  bool operator!=(const Rational& right) const {
    return !((*this) == right);
  }

  bool operator<(const Rational& right) const {
    return this->numerator * right.denominator
      < right.numerator* this->denominator;
  }

  bool operator>(const Rational& right) const {
    return right < (*this);
  }

  bool operator>=(const Rational& right) const {
    return !(*this < right);
  }

  bool operator<=(const Rational& right) const {
    return !(*this > right);
  }

  RationalX operator-() const {
    return RationalX(-numerator, denominator);
  }

  Rational& operator=(const Rational&) = default;
  Rational& operator=(const Integer& number) {
    numerator = number;
    denominator = static_cast<Integer>(1);
    return *this;
  }

  Rational& operator+=(const Rational& right) {
    if (right.numerator == 0) {// (0 + right)
      return *this;
    }
    if (this->numerator == 0) {// (this + 0)
      this->numerator = right.numerator;
      this->denominator = right.denominator;
      return *this;
    }

    /**
     *  Lnum     Rnum     Lnum*Rden + Rnum*Lden
     * ------ + ------ = -----------------------
     *  Lden     Rden           Lden*Rden
     *      Lnum*(RdenC*denC) + Rnum*(LdenC*denC)     (Lnum*RdenC + Rnum*LdenC)*denC
     *   = --------------------------------------- = --------------------------------
     *            (LdenC*denC)*(RdenC*denC)             (Lden*RdenC)*denC
     *      Lnum*(Rden/denC) + Rnum*(Lden/denC)     denC
     *   = ------------------------------------- * ------
     *      Lden * Rden/denC                        denC
    */
    Integer denominator_gcd = gcd(this->denominator, right.denominator);
    this->numerator = this->numerator * (right.denominator / denominator_gcd) + right.numerator * (this->denominator / denominator_gcd);
    this->denominator = this->denominator * (right.denominator / denominator_gcd);

    Integer common_divisor = gcd(this->numerator, this->denominator);
    this->numerator /= common_divisor;
    this->denominator /= common_divisor;
    return *this;
  }
  Rational& operator+=(const Integer& num) {
    return *this += Rational(num);
  }

  Rational& operator-=(const Rational& right) {
    return *this += (-right);
  }
  Rational& operator-=(const Integer& num) {
    return *this -= Rational(num);
  }

  Rational& operator*=(const Rational& right) {
    if (this->numerator == 0 || right.numerator == 0) {
      this->numerator = 0;
      this->denominator = 1;
      return *this;
    }

    /**
     *  Lnum     Rnum     Lnum*Rnum
     * ------ * ------ = -----------
     *  Lden     Rden     Lden*Rden
     *      LnumC*gcdC * RnumD*gcdD     LnumC * RnumD *gcdC * gcdD
     *   = ------------------------- = -----------------------------
     *      LdenD*gcdD * RdenC*gcdC     LdenD * RdenC * gcdC * gcdD
     *      Lnum/gcdC * Rnum/gcdD     gcdC     gcdD
     *   = ----------------------- * ------ * ------
     *      Lden/gcdD * Rden/gcdC     gcdC     gcdD
    */
    Integer common_divisor1 = gcd(this->numerator, right.denominator);
    Integer common_divisor2 = gcd(right.numerator, this->denominator);
    this->numerator = (this->numerator / common_divisor1) * (right.numerator / common_divisor2);
    this->denominator = (this->denominator / common_divisor2) * (right.denominator / common_divisor1);

    Integer common_divisor = gcd(this->numerator, this->denominator);
    this->numerator /= common_divisor;
    this->denominator /= common_divisor;
    return *this;
  }
  Rational& operator*=(const Integer& num) {
    return *this *= Rational(num);
  }

  Rational& operator/=(const Rational& right) {
    return *this *= Rational(right.denominator, right.numerator);
  }
  Rational& operator/=(const Integer& num) {
    return *this *= Rational(static_cast<Integer>(1), num);
  }

  Rational operator+(const Rational& right) const {
    return Rational(*this) += right;
  }
  Rational operator+(const Integer& num) const {
    return Rational(*this) += Rational(num);
  }

  Rational operator-(const Rational& right) const {
    return Rational(*this) -= right;
  }
  Rational operator-(const Integer& num) const {
    return Rational(*this) -= Rational(num);
  }

  Rational operator*(const Rational& right) const {
    return Rational(*this) *= right;
  }
  Rational operator*(const Integer& num) const {
    return Rational(*this) *= Rational(num);
  }

  Rational operator/(const Rational& right) const {
    return RationalX(*this) /= right;
  }
  Rational operator/(const Integer& num) const {
    return RationalX(*this) /= Rational(num);
  }

  Rational& operator++() {
    numerator += denominator;
    return *this;
  }

  Rational operator++(int) {
    Rational copied = *this;
    ++(*this);
    return copied;
  }

  explicit
  RationalX(float value) {
    /**
     * <source>
     *   0 00000001 10101101011010110101000
     * </source>
     * <first>
     *   <tips> IEEE754-floating-formula: (-1)^S * (1+0.Fraction) * 2^(Exponent-Bias) </tips>
     *
     *   (-1)^0 * (1 + 0.10101101011010110101000) * 2^(00000001 - Bias)
     *     = 1 * 1.10101101011010110101000 * pow(2, _Exp)
     *     = 1 * 0.110101101011010110101000 * pow(2, _Exp)
     *     = 1 * 110101101011010110101000/pow(2,_Mn) * pow(2, _Exp)
     * </first>
     * <second>
     *   <tips> pow(2, X) = (1 << X) </tips>
     *
     *   _Nx     110101101011010110101000
     *   ---- = -------------------------- * ( 1 << _Exp )
     *   _Dx           1 << _Mn
     *
     *           110101101011010110101000 << 1
     *          ------------------------------- * ( 1 << (_Exp - 1) )
     *                 1 << _Mn
     *
     *           110101101011010110101000
     *          ------------------------------- * ( 1 << (_Exp - 1) )
     *                 1 << (_Mn-1)
     *
     *           110101101011010110101000
     *          ------------------------------- * (1 << 0)
     *                 1 << (_Mn - _Exp)
     * </second>
    */
    constexpr unsigned int sign_mask
      = 0b10000000000000000000000000000000;
    constexpr unsigned int exponent_mask
      = 0b01111111100000000000000000000000;
    constexpr unsigned int mantissa_mask
      = 0b00000000011111111111111111111111;
    constexpr unsigned int hidden_significant
      = 0b00000000100000000000000000000000;
    constexpr char exp2_bias = 127;

    unsigned int value_bits = reinterpret_cast<uint32_t&>(value);
    unsigned int exp2_bits = (value_bits & exponent_mask) >> 23;
    unsigned int signifi_bits = value_bits & mantissa_mask | hidden_significant;
    char exp2 = reinterpret_cast<char&>(exp2_bits) - exp2_bias;
    exp2 -= 23;

    // *this = significant * pow(2,exp2)
    if (exp2 > 0) {
      this->numerator = static_cast<Integer>(signifi_bits << exp2);
      this->denominator = static_cast<Integer>(1);
    }
    else if (exp2 < 0) {
      static_assert(sizeof(Integer) >= sizeof(unsigned int), "ratianal(float)");
      this->numerator = static_cast<Integer>(signifi_bits);
      this->denominator = static_cast<Integer>(1) << (-exp2);
    }
    else {
      this->numerator = static_cast<Integer>(signifi_bits);
      this->denominator = static_cast<Integer>(1);
    }

    // *this *= ~sign
    if ((value_bits & sign_mask) != 0) {
      this->numerator = -this->numerator;
    }

    // divide greater_common_divisor
    Integer common_divisor = gcd(numerator, denominator);
    numerator /= common_divisor;
    denominator /= common_divisor;
  }

  explicit
  RationalX(double value) {
    constexpr unsigned long long sign_mask
      = 0b1000000000000000000000000000000000000000000000000000000000000000;
    constexpr unsigned long long exponent_mask
      = 0b0111111111110000000000000000000000000000000000000000000000000000;
    constexpr unsigned long long mantissa_mask
      = 0b0000000000001111111111111111111111111111111111111111111111111111;
    constexpr unsigned long long hidden_significant
      = 0b0000000000010000000000000000000000000000000000000000000000000000;
    constexpr short exp2_bias = 1023;

    // seperate bits
    unsigned long long value_bits = reinterpret_cast<unsigned long long&>(value);
    unsigned long long exp2_bits = (value_bits & exponent_mask) >> 52;
    unsigned long long signifi_bits = value_bits & mantissa_mask | hidden_significant;
    short exp2 = reinterpret_cast<short&>(exp2_bits) - exp2_bias;
    exp2 -= 52;

    // *this = significant * pow(2,exp2)
    if (exp2 > 0) {
      this->numerator = static_cast<Integer>(signifi_bits << exp2);
      this->denominator = static_cast<Integer>(1);
    }
    else if (exp2 < 0) {
      static_assert(sizeof(Integer) >= sizeof(unsigned long long), "ratianal(double)");
      this->numerator = static_cast<Integer>(signifi_bits);
      this->denominator = static_cast<Integer>(1) << (-exp2);
    }
    else {
      this->numerator = static_cast<Integer>(signifi_bits);
      this->denominator = static_cast<Integer>(1);
    }

    // *this *= ~sign
    if ((value_bits & sign_mask) != 0) {
      this->numerator = -this->numerator;
    }

    // divide greater_common_divisor
    Integer common_divisor = gcd(numerator, denominator);
    numerator /= common_divisor;
    denominator /= common_divisor;
  }
};

template<typename __number_t>
struct RationalX<__number_t, false> {
  using Rational = RationalX<__number_t, false>;
  using Integer = __number_t;

  Integer numerator;
  Integer denominator;

  constexpr RationalX() : numerator(0), denominator(0) {}

  constexpr RationalX(Integer num, Integer den = static_cast<Integer>(1)) : numerator(num), denominator(den) {
    Integer common_divisor = gcd(this->numerator, this->denominator);
    this->numerator /= common_divisor;
    this->denominator /= common_divisor;
  }

  bool operator==(const Rational& right) const {
    // assert( simplest(this) && simplest(right) )
    return this->numerator == right.numerator
      && right.denominator == this->denominator;
  }

  bool operator!=(const Rational& right) const {
    return !((*this) == right);
  }

  bool operator<(const Rational& right) const {
    return this->numerator * right.denominator
      < right.numerator* this->denominator;
  }

  bool operator>(const Rational& right) const {
    return right < (*this);
  }

  bool operator>=(const Rational& right) const {
    return !(*this < right);
  }

  bool operator<=(const Rational& right) const {
    return !(*this > right);
  }

  Rational& operator=(const Rational&) = default;
  Rational& operator=(Integer number) {
    numerator = number;
    denominator = static_cast<Integer>(1);
    return *this;
  }

  Rational& operator+=(const Rational& right) {
    if (right.numerator == 0) {// (0 + right)
      return *this;
    }
    if (this->numerator == 0) {// (this + 0)
      this->numerator = right.numerator;
      this->denominator = right.denominator;
      return *this;
    }

    /**
     *  Lnum     Rnum     Lnum*Rden + Rnum*Lden
     * ------ + ------ = -----------------------
     *  Lden     Rden           Lden*Rden
     *      Lnum*(RdenC*denC) + Rnum*(LdenC*denC)     (Lnum*RdenC + Rnum*LdenC)*denC
     *   = --------------------------------------- = --------------------------------
     *            (LdenC*denC)*(RdenC*denC)             (Lden*RdenC)*denC
     *      Lnum*(Rden/denC) + Rnum*(Lden/denC)     denC
     *   = ------------------------------------- * ------
     *      Lden * Rden/denC                        denC
    */
    Integer denominator_gcd = gcd(this->denominator, right.denominator);
    this->numerator = this->numerator * (right.denominator / denominator_gcd) + right.numerator * (this->denominator / denominator_gcd);
    this->denominator = this->denominator * (right.denominator / denominator_gcd);

    Integer common_divisor = gcd(this->numerator, this->denominator);
    this->numerator /= common_divisor;
    this->denominator /= common_divisor;
    return *this;
  }
  Rational& operator+=(const Integer& num) {
    return *this += Rational(num);
  }

  Rational& operator-=(const Rational& right) {
    if (this->numerator == 0) {// (this - 0)
      return *this;
    }

    Integer denominator_gcd = gcd(this->denominator, right.denominator);
    if (this->numerator * (right.denominator / denominator_gcd) < right.numerator * (this->denominator / denominator_gcd)) {
      throw std::underflow_error("calculation::RationalX<number_t,unsigned>::operator-=( ... )");
    }
    this->numerator = this->numerator * (right.denominator / denominator_gcd) - right.numerator * (this->denominator / denominator_gcd);
    this->denominator = this->denominator * (right.denominator / denominator_gcd);

    Integer common_divisor = gcd(this->numerator, this->denominator);
    this->numerator /= common_divisor;
    this->denominator /= common_divisor;
    return *this;
  }
  Rational& operator-=(const Integer& num) {
    return *this -= Rational(num);
  }

  Rational& operator*=(const Rational& right) {
    if (this->numerator == 0 || right.numerator == 0) {
      this->numerator = 0;
      this->denominator = 1;
      return *this;
    }

    /**
     *  Lnum     Rnum     Lnum*Rnum
     * ------ * ------ = -----------
     *  Lden     Rden     Lden*Rden
     *      LnumC*gcdC * RnumD*gcdD     LnumC * RnumD *gcdC * gcdD
     *   = ------------------------- = -----------------------------
     *      LdenD*gcdD * RdenC*gcdC     LdenD * RdenC * gcdC * gcdD
     *      Lnum/gcdC * Rnum/gcdD     gcdC     gcdD
     *   = ----------------------- * ------ * ------
     *      Lden/gcdD * Rden/gcdC     gcdC     gcdD
    */
    Integer common_divisor1 = gcd(this->numerator, right.denominator);
    Integer common_divisor2 = gcd(right.numerator, this->denominator);
    this->numerator = (this->numerator / common_divisor1) * (right.numerator / common_divisor2);
    this->denominator = (this->denominator / common_divisor2) * (right.denominator / common_divisor1);

    Integer common_divisor = gcd(this->numerator, this->denominator);
    this->numerator /= common_divisor;
    this->denominator /= common_divisor;
    return *this;
  }
  Rational& operator*=(const Integer& num) {
    return *this *= Rational(num);
  }

  Rational& operator/=(const Rational& right) {
    return *this *= Rational(right.denominator, right.numerator);
  }
  Rational& operator/=(const Integer& num) {
    return *this *= Rational(static_cast<Integer>(1), num);
  }

  Rational operator+(const Rational& right) const {
    return Rational(*this) += right;
  }
  Rational operator+(const Integer& num) const {
    return Rational(*this) += Rational(num);
  }

  Rational operator-(const Rational& right) const {
    return Rational(*this) -= right;
  }
  Rational operator-(const Integer& num) const {
    return Rational(*this) -= Rational(num);
  }

  Rational operator*(const Rational& right) const {
    return Rational(*this) *= right;
  }
  Rational operator*(const Integer& num) const {
    return Rational(*this) *= Rational(num);
  }

  Rational operator/(const Rational& right) const {
    return Rational(*this) /= right;
  }
  Rational operator/(const Integer& num) const {
    return Rational(*this) /= Rational(num);
  }

  Rational& operator++() {
    numerator += denominator;
    return *this;
  }

  Rational operator++(int) {
    Rational copied = *this;
    ++(*this);
    return copied;
  }

  explicit
  RationalX(float value) {
    if (value < 0) {
      throw std::underflow_error("calculation::RationalX<number_t,unsigned>::RationalX(float)");
    }

    constexpr unsigned int exponent_mask
      = 0b01111111100000000000000000000000;
    constexpr unsigned int mantissa_mask
      = 0b00000000011111111111111111111111;
    constexpr unsigned int hidden_significant
      = 0b00000000100000000000000000000000;
    constexpr char exp2_bias = 127;

    unsigned int value_bits = reinterpret_cast<uint32_t&>(value);
    unsigned int exp2_bits = (value_bits & exponent_mask) >> 23;
    unsigned int signifi_bits = value_bits & mantissa_mask | hidden_significant;
    char exp2 = reinterpret_cast<char&>(exp2_bits) - exp2_bias;
    exp2 -= 23;

    // *this = significant * pow(2,exp2)
    if (exp2 > 0) {
      this->numerator = static_cast<Integer>(signifi_bits << exp2);
      this->denominator = static_cast<Integer>(1);
    }
    else if (exp2 < 0) {
      static_assert(sizeof(Integer) >= sizeof(unsigned int), "ratianal(float)");
      this->numerator = static_cast<Integer>(signifi_bits);
      this->denominator = static_cast<Integer>(1) << (-exp2);
    }
    else {
      this->numerator = static_cast<Integer>(signifi_bits);
      this->denominator = static_cast<Integer>(1);
    }

    // divide greater_common_divisor
    Integer common_divisor = gcd(numerator, denominator);
    numerator /= common_divisor;
    denominator /= common_divisor;
  }

  explicit
  RationalX(double value) {
    if (value < 0) {
      throw std::underflow_error("calculation::RationalX<number_t,unsigned>::RationalX(double)");
    }

    constexpr unsigned long long exponent_mask
      = 0b0111111111110000000000000000000000000000000000000000000000000000;
    constexpr unsigned long long mantissa_mask
      = 0b0000000000001111111111111111111111111111111111111111111111111111;
    constexpr unsigned long long hidden_significant
      = 0b0000000000010000000000000000000000000000000000000000000000000000;
    constexpr short exp2_bias = 1023;

    // seperate bits
    unsigned long long value_bits = reinterpret_cast<unsigned long long&>(value);
    unsigned long long exp2_bits = (value_bits & exponent_mask) >> 52;
    unsigned long long signifi_bits = value_bits & mantissa_mask | hidden_significant;
    short exp2 = reinterpret_cast<short&>(exp2_bits) - exp2_bias;
    exp2 -= 52;

    // *this = significant * pow(2,exp2)
    if (exp2 > 0) {
      this->numerator = static_cast<Integer>(signifi_bits << exp2);
      this->denominator = static_cast<Integer>(1);
    }
    else if (exp2 < 0) {
      static_assert(sizeof(Integer) >= sizeof(unsigned long long), "ratianal(double)");
      this->numerator = static_cast<Integer>(signifi_bits);
      this->denominator = static_cast<Integer>(1) << (-exp2);
    }
    else {
      this->numerator = static_cast<Integer>(signifi_bits);
      this->denominator = static_cast<Integer>(1);
    }

    // divide greater_common_divisor
    Integer common_divisor = gcd(numerator, denominator);
    numerator /= common_divisor;
    denominator /= common_divisor;
  }
};

template<typename Integer> inline
RationalX<Integer> operator+(RationalX<Integer> a, int b) {
  return a + RationalX<Integer>(static_cast<Integer>(b));
}
template<typename Integer> inline
RationalX<Integer> operator+(RationalX<Integer> a, long long b) {
  return a + RationalX<Integer>(static_cast<Integer>(b));
}
template<typename Integer> inline
RationalX<Integer> operator+(RationalX<Integer> a, unsigned int b) {
  return a + RationalX<Integer>(static_cast<Integer>(b));
}
template<typename Integer> inline
RationalX<Integer> operator+(RationalX<Integer> a, unsigned long long b) {
  return a + RationalX<Integer>(static_cast<Integer>(b));
}

template<typename Integer> inline
RationalX<Integer> operator-(RationalX<Integer> a, int b) {
  return a - RationalX<Integer>(static_cast<Integer>(b));
}
template<typename Integer> inline
RationalX<Integer> operator-(RationalX<Integer> a, long long b) {
  return a - RationalX<Integer>(static_cast<Integer>(b));
}
template<typename Integer> inline
RationalX<Integer> operator-(RationalX<Integer> a, unsigned int b) {
  return a - RationalX<Integer>(static_cast<Integer>(b));
}
template<typename Integer> inline
RationalX<Integer> operator-(RationalX<Integer> a, unsigned long long b) {
  return a - RationalX<Integer>(static_cast<Integer>(b));
}

template<typename Integer> inline
RationalX<Integer> operator*(RationalX<Integer> a, int b) {
  return a * RationalX<Integer>(static_cast<Integer>(b));
}
template<typename Integer> inline
RationalX<Integer> operator*(RationalX<Integer> a, long long b) {
  return a * RationalX<Integer>(static_cast<Integer>(b));
}
template<typename Integer> inline
RationalX<Integer> operator*(RationalX<Integer> a, unsigned int b) {
  return a * RationalX<Integer>(static_cast<Integer>(b));
}
template<typename Integer> inline
RationalX<Integer> operator*(RationalX<Integer> a, unsigned long long b) {
  return a * RationalX<Integer>(static_cast<Integer>(b));
}

template<typename Integer> inline
RationalX<Integer> operator/(RationalX<Integer> a, int b) {
  return a / RationalX<Integer>(static_cast<Integer>(b));
}
template<typename Integer> inline
RationalX<Integer> operator/(RationalX<Integer> a, long long b) {
  return a / RationalX<Integer>(static_cast<Integer>(b));
}
template<typename Integer> inline
RationalX<Integer> operator/(RationalX<Integer> a, unsigned int b) {
  return a / RationalX<Integer>(static_cast<Integer>(b));
}
template<typename Integer> inline
RationalX<Integer> operator/(RationalX<Integer> a, unsigned long long b) {
  return a / RationalX<Integer>(static_cast<Integer>(b));
}

template<typename Integer> inline
RationalX<Integer> abs(const RationalX<Integer>& x) {
  if constexpr ( std::is_signed_v<Integer> ) {
    return RationalX<Integer>{ abs(x.numerator), x.denominator };
  } else {
    return x;
  }
}

template<typename _Elem, typename _Traits, 
  typename Integer> inline
std::basic_ostream<_Elem,_Traits>& operator<<(std::basic_ostream<_Elem,_Traits>& _Ostr, const RationalX<Integer>& _R) {
    return _Ostr << '(' << _R.numerator << ',' << _R.denominator << ')';
}
 

/**
 * @author A.R.Barnett
 * 
 * @param num_it ArrayIterator or FunctionObject
 * 
 * @param den_it ArrayIterator or FunctionObject
 * 
 * @definition 
 *   f(x) = b[0]
 *     + a[1]/(b[1]
 *       + a[2]/(b[2]
 *         + a[3]/(b[3]
 *           + a[4]/(b[4]
 *             + ... 
 *                + a[n]/b[n] ))))
 * 
 * @alternative 
 *   f(x) = A[j]/B[j], when limit j->inf   :forward_recurrence
 *   
 *   f(x) = b[0] + sum<j=1,inf>( term[j-1]*(b[j]/(B[j]/B[j-1]) - 1) )  :sum_series
 * 
 *   f(x) = b[0]*product<j=1,inf>( (A[j]/A[j-1])/(B[j]/B[j-1]) )  :product_series
 * 
 * @reference
 *   tutorial = "http://www.maths.surrey.ac.uk/hosted-sites/R.Knott/Fibonacci/cfCALC.html"
*/
/*std::vector<double> a(100);
  std::vector<double> b(100);
  double x = 1.57;*/
  // cosine function
  /*b[0] = a[0] = 0;
  b[1] = 1; a[1] = 1;
  b[2] = 2-x*x; a[2] = x*x;
  b[3] = 3*4-x*x; a[3] = 2*x*x;
  for (size_t i = 4; i != a.size(); ++i) {
    b[i] = ((i-1)*2-1)*((i-1)*2)-x*x;
    a[i] = ((i-1)*2-3)*((i-1)*2-2)*x*x;
  }*/

  // tangent function
  /*b[0] = a[0] = 0;
  b[1] = 1/x; a[1] = 1;
  for (size_t i = 2; i != a.size(); ++i) {
    b[i] = (2*i-1)/x;
    a[i] = -1;
  }*/
template<typename Iterator>
auto evaluate_continued_fraction(Iterator num_it, Iterator den_it) {
/**
 * @continued_fraction
 *   f(x)
 *   = b[0]
 *     + a[1]/(b[1]
 *       + a[2]/(b[2]
 *         + a[3]/(b[3]
 *           + a[4]/(b[4]
 *             + ... 
 *                + a[n]/b[n] ))))
 * 
 * @continued_fraction forward recurrence
 *   limit j->inf
 *                  A[j]     A[j-1]*b[j] + A[j-2]*a[j]
 *   f(x) = Y[j] = ------ = ---------------------------
 *                  B[j]     B[j-1]*b[j] + B[j-2]*a[j]
 * 
 *                   1     A[-1]
 *   Y[-1] =  0  =  --- = -------
 *                   0     B[-1]  
 * 
 *                  b[0]     A[0]
 *   Y[0] = b[0] = ------ = ------
 *                   1       B[0]
 * 
 *           A[0]*b[1] + A[-1]*a[1]     b[0]*b[1] + a[1]
 *   Y[1] = ------------------------ = ------------------
 *           B[0]*b[1] + B[-1]*a[1]     b[1]
 * 
 * @continued_fraction to sum_series(Steed.) or product_series(Lentz.)
 * 
 *   First, we should find some properties of the variables 
 * 
 *                A[j]       A[j-1]*b[j] + A[j-2]*a[j]            A[j-2]
 *     Arel[j] = -------- = --------------------------- = b[j] + --------*a[j] = b[j] + a[j]/Arel[j-1]
 *                A[j-1]              A[j-1]                      A[j-1]
 * 
 *                B[j]       B[j-1]*b[j] + B[j-2]*a[j]            B[j-2]
 *     Brel[j] = -------- = --------------------------- = b[j] + --------*a[j] = b[j] + a[j]/Brel[j-1]
 *                B[j-1]              B[j-1]                      B[j-1]
 *            
 *   Then, we want to get the recurrence relation of Y[..]
 * 
 *      Y[j]       A[j]   A[j-1]     A[j]     B[j-1]
 *     -------- = ------/-------- = --------*--------  :This is derivation of Lentz's mathod
 *      Y[j-1]     B[j]   B[j-1]     A[j-1]   B[j]
 * 
 *     And
 * 
 *             A[j-1]*b[j] + A[j-2]*a[j]
 *     Y[j] = ---------------------------
 *             B[j-1]*b[j] + B[j-2]*a[j]
 * 
 *          Here has two kind division, The one                    Another
 * 
 *             ( A[j-1]*b[j] + A[j-2]*a[j] )/A[j-1]                   ( A[j-1]*b[j] + A[j-2]*a[j] )/B[j-1]
 *          = --------------------------------------               = --------------------------------------
 *             ( B[j-1]*b[j] + B[j-2]*a[j] )/A[j-1]                   ( B[j-1]*b[j] + B[j-2]*a[j] )/B[j-1]
 * 
 *                          Arel[j]                                   Y[j-1]*b[j] + A[j-2]/B[j-1]*a[j]
 *          = ------------------------------------                 = ----------------------------------
 *             1/Y[j-1]*b[j] + B[j-2]/A[j-1]*a[j]                                Brel[j]
 * 
 *                          Arel[j]                                   Y[j-1]*b[j] + A[j-2]/B[j-2]*B[j-2]/B[j-1]*a[j]
 *          = --------------------------------------------------   = ------------------------------------------------
 *             1/Y[j-1]*b[j] + B[j-2]/A[j-2]*A[j-2]/A[j-1]*a[j]                  Brel[j]
 * 
 *                          Arel[j]                                   Y[j-1]*b[j] + Y[j-2]/Brel[j-1]*a[j]
 *          = ---------------------------------------              = -------------------------------------
 *             b[j]/Y[j-1] + a[j]/(Y[j-2]*Arel[j-1])                             Brel[j]
 * 
 *   Now, we derivate the difference between Y[j] and Y[j-1]
 * 
 *                    Y[j-1]*b[j] + Y[j-2]/Brel[j-1]*a[j]
 *   Y[j] - Y[j-1] = ------------------------------------- - Y[j-1]
 *                             Brel[j]
 * 
 *                    Y[j-1]*b[j] - Y[j-1]*Brel[j] + Y[j-2]*(Brel[j]-b[j])
 *                 = ------------------------------------------------------   :a[j]/Brel[j-1] = Brel[j]-b[j]
 *                             Brel[j]
 * 
 *                    Y[j-1](b[j] - Brel[j]) + Y[j-2]*(Brel[j] - b[j])
 *                 = --------------------------------------------------
 *                             Brel[j]
 * 
 *                    (Y[j-1] - Y[j-2])*(b[j] - Brel[j])
 *                 = ------------------------------------                     :Y[j-2]*(Brel[j] - b[j]) = -Y[j-2]*(b[j] - Brel[j])
 *                             Brel[j]
 * 
 *                 = (Y[j-1] - Y[j-2])*(b[j]/Brel[j] - 1)                     :This is derivation of Steed's mathod
*/
	using Number = decltype(*num_it/(*den_it));
	const Number eps = std::numeric_limits<Number>::epsilon();
#if 0
	Number num = *num_it++;
	Number den = *den_it++;
	Number series = den;

	num = *num_it++;
	den = *den_it++;
	Number Brel = den;
	Number term = num/den;
	do {
		series += term;
		num = *num_it++;
		den = *den_it++;
		Brel = den + num/Brel;
		term = term*(den/Brel - 1);
	} while ( abs(term) >= eps*abs(series) );
	
	return series;
#else
	Number num = *num_it++;
	Number den = *den_it++;
	Number series = den; if (series == 0) series = 10e-50;
	
	num = *num_it++;
	den = *den_it++;
	Number Arel = den + num/series; if(Arel == 0) Arel = 1e-50;
	Number Brel = den;              if(Brel == 0) Brel = 1e-50;
	Number term = Arel/Brel;
	do {
		series *= term;
		num = *num_it++;
		den = *den_it++;
		Arel = den + num/Arel; if(Arel == 0) Arel = 1e-50;
		Brel = den + num/Brel; if(Brel == 0) Brel = 1e-50;
		term = Arel/Brel;
	} while ( abs(1 - term) >= eps );

	return series;
#endif
};
  }// namespace calculation