#pragma once
/*{ "calculation/fundamental/rational":{
  "Description": "infinite precision calculation"
  "License": "Please identity Author"
  "Author": "LongJiangnan",
  "Date": "2019-2021",
} }*/


#include <numeric>// std::gcd
#include <math.h>// _CSTD abs
#include <iosfwd>// std::basic_ostream
#include <exception>// std::underflow_error, std::overflow_error ...
namespace calculation 
{

using _CSTD abs;

using std::gcd;

using std::lcm;

template<typename UnsignedInteger>
void add_unsigned_rational(
  UnsignedInteger& this_numerator, 
  UnsignedInteger& this_denominator, 
  const UnsignedInteger& right_numerator, 
  const UnsignedInteger& right_denominator,
  const UnsignedInteger upper_boundary)
{
  if ( this_numerator == 0 ) {
    this_numerator = right_numerator;
    this_denominator = right_denominator;
    return;
  }
  if ( right_numerator == 0 ) {
    // do nothing
    return;
  }

  /**
   *  numL     numR     numL*denR     numR*denL
   * ------ + ------ = ----------- + -----------
   *  denL     denR     denL*denR     denR*denL
   * 
   *                    numL*denR + numR*denL
   *                 = -----------------------
   *                    denL*denR
   * 
   *                    (numL*denR + numR*denL) / divC
   *                 = --------------------------------  :denL%divC = 0, denR%divC = 0
   *                    denL*denR / divC
   * 
   * assert( (numL*denR + numR*denL) / divC < MAX )
   * assert( numL*denR/divC <= MAX 
   *      && numR*denL/divC <= MAX
   *      && numL*denR/divC + numR*denL/divC < MAX )
   * assert( denR/divC <= MAX/numL
   *      && denL/divC <= MAX/numR
   *      && numL*denR/divC + numR*denL/divC < MAX )
   * 
   * assert( denL*denR/divC <= MAX )
   * assert( denL/divC <= MAX/denR )
  */
  UnsignedInteger
  common_divisor = gcd(this_denominator, right_denominator);
    
  UnsignedInteger
  another_numerator = right_denominator / common_divisor;
  if ( another_numerator > upper_boundary / this_numerator ) {
    throw std::overflow_error("add_unsigned_rational(...)");
  }
  another_numerator *= this_numerator;

  this_numerator = this_denominator / common_divisor;
  if ( this_numerator > upper_boundary / right_numerator ) {
    throw std::overflow_error("add_unsigned_rational(...)");
  }
  this_numerator *= right_numerator;

  if ( this_numerator > upper_boundary - another_numerator ) {
    throw std::overflow_error("add_unsigned_rational(...)");
  }
  this_numerator += another_numerator;

  this_denominator = this_denominator / common_divisor;
  if ( this_denominator > upper_boundary / right_denominator ) {
    throw std::overflow_error("add_unsigned_rational(...)");
  }
  this_denominator *= right_denominator;

  common_divisor = gcd(this_numerator, this_denominator);
  this_numerator /= common_divisor;
  this_denominator /= common_divisor;
}

template<typename UnsignedInteger>
void sub_unsigned_rational(
  UnsignedInteger& this_numerator, 
  UnsignedInteger& this_denominator, 
  const UnsignedInteger& right_numerator, 
  const UnsignedInteger& right_denominator,
  const UnsignedInteger upper_boundary)
{
  if ( this_numerator == 0 ) {
    this_numerator = right_numerator;
    this_denominator = right_denominator;
    return;
  }
  if ( right_numerator == 0 ) {
    // do nothing
    return;
  }

  /**
   *  numL     numR     numL*denR     numR*denL
   * ------ - ------ = ----------- - -----------
   *  denL     denR     denL*denR     denR*denL
   * 
   *                    numL*denR - numR*denL
   *                 = -----------------------
   *                    denL*denR
   * 
   *                    (numL*denR - numR*denL) / divC
   *                 = --------------------------------  :denL%divC = 0, denR%divC = 0
   *                    denL*denR / divC
   * 
   * assert( (numL*denR - numR*denL) / divC >= 0)
   * assert( numL*denR/divC <= MAX
   *      && numR*denL/divC <= MAX
   *      && numL*denR/divC >= 0
   *      && numR*denL/divC >= 0
   *      && numL*denR/divC >= numR*denL/divC )
   * assert( denR/divC <= MAX/numL
   *      && denL/divC <= MAX/numR
   *      && numL*denR/divC >= numR*denL/divC )
   * 
   * assert( denL*denR/divC <= MAX )
   * assert( denL/divC <= MAX/denR )
  */
  UnsignedInteger
  common_divisor = gcd(this_denominator, right_denominator);
    
  UnsignedInteger
  another_numerator = right_denominator / common_divisor;
  if ( another_numerator > upper_boundary / this_numerator ) {
    throw std::overflow_error("sub_unsigned_rational(...)");
  }
  another_numerator *= this_numerator;

  this_numerator = this_denominator / common_divisor;
  if ( this_numerator > upper_boundary / right_numerator ) {
    throw std::overflow_error("sub_unsigned_rational(...)");
  }
  this_numerator *= right_numerator;

  if ( another_numerator < this_numerator ) {
    throw std::underflow_error("sub_unsigned_rational(...)");
  }
  this_numerator = another_numerator - this_numerator;

  this_denominator = this_denominator / common_divisor;
  if ( this_denominator > upper_boundary / right_denominator ) {
    throw std::overflow_error("sub_unsigned_rational(...)");
  }
  this_denominator *= right_denominator;

  common_divisor = gcd(this_numerator, this_denominator);
  this_numerator /= common_divisor;
  this_denominator /= common_divisor;
}

template<typename UnsignedInteger>
void mul_unsigned_rational(
  UnsignedInteger& this_numerator,
  UnsignedInteger& this_denominator,
  const UnsignedInteger& right_numerator,
  const UnsignedInteger& right_denominator,
  const UnsignedInteger upper_boundary) 
{
  if ( this_numerator == 0 || right_numerator == 0 ) {
    this_numerator = 0;
    this_denominator = 1;
    return;
  }
			
  /**
   *  numL     numR     numL*numR
   * ------ * ------ = -----------
   *  denL     denR     denL*denR
   * 
   * assert( numL*numR <= MAX )
   * assert( numL <= MAX/numR )
   * 
   * assert( denL*denR <= MAX )
   * assert( denL <= MAX/denR )
  */

  if ( this_numerator > upper_boundary/right_numerator ) {
    throw std::overflow_error("mul_unsigned_rational");
  }
  this_numerator *= right_numerator;

  if ( this_denominator > upper_boundary/right_denominator ) {
    throw std::overflow_error("mul_unsigned_rational");
  }
  this_denominator *= right_denominator;

  UnsignedInteger
  common_divisor = gcd(this_numerator, this_denominator);
  this_numerator /= common_divisor;
  this_denominator /= common_divisor;
}

template<typename UnsignedInteger>
bool less_unsigned_rational(
  const UnsignedInteger& this_numerator,
  const UnsignedInteger& this_denominator,
  const UnsignedInteger& right_numerator,
  const UnsignedInteger& right_denominator,
  const UnsignedInteger upper_boundary) {
  if (this_numerator == 0) {
    return right_numerator != 0;
  }
  if (right_numerator == 0) {
    return false;
  }

  UnsignedInteger
  common_divisor = gcd(this_denominator, right_denominator);
    
  UnsignedInteger
  the_numerator = right_denominator / common_divisor;
  if ( the_numerator > upper_boundary / this_numerator ) {
    throw std::overflow_error("sub_unsigned_rational(...)");
  }
  the_numerator *= this_numerator;

  UnsignedInteger
  another_numerator = this_denominator / common_divisor;
  if ( another_numerator > upper_boundary / right_numerator ) {
    throw std::overflow_error("sub_unsigned_rational(...)");
  }
  another_numerator *= right_numerator;

  return the_numerator < another_numerator;
}


/** TEMPLATE rational number */
template<size_t __b, bool __sb = true, bool __opt = false>
struct rationalX {
  // no implementation
};

/** TEMPLATE SPECIAL 32bits unsigned rational */
template<>
struct rationalX<32,false, true> {
  typedef unsigned int unsigned_integer_type;
  unsigned_integer_type numerator;
  unsigned_integer_type denominator;

  rationalX() = default;
  rationalX(unsigned int _Nx, unsigned int _Dx = 1) noexcept {
    unsigned int _Cd = gcd(_Nx, _Dx);
    this->numerator = _Nx/_Cd;
    this->denominator = _Dx/_Cd;
  }
  rationalX(unsigned long long _Nx, unsigned long long _Dx = 1){
    // reduce
    unsigned long long _Cd = gcd(_Nx, _Dx);
    _Nx /= _Cd;
    _Dx /= _Cd;
    // check upper boundaty
    unsigned long long _Max = static_cast<unsigned long long>( std::numeric_limits<unsigned_integer_type>::max() );
    if ( _Nx > _Max || _Dx > _Max ) {
      throw std::overflow_error("calculation::rationalX<32,false, true>"
        "::rationalX(unsigned long long, unsigned long long)");
    }
    // encode
    this->numerator = static_cast<unsigned_integer_type>(_Nx);
    this->denominator = static_cast<unsigned_integer_type>(_Dx);
  }
  rationalX(int _Nx, int _Dx = 1) {
    if ( _Nx < 0 || _Dx < 0 ) {
      throw std::overflow_error("calculation::rationalX<32,false, true>"
        "::rationalX(int, int)");
    }
    int _Cd = gcd(_Nx, _Dx);
    this->numerator = _Nx / _Cd;
    this->denominator = _Dx / _Cd;
  }
  rationalX(long long _Nx, long long _Dx = 1) {
    // check lower boundary
    if ( _Nx < 0 || _Dx < 0 ) {
      throw std::overflow_error("calculation::rationalX<32,false, true>"
        "::rationalX(long long, long long)");
    }
    // reduce
    long long _Cd = gcd(_Nx, _Dx);
    _Nx /= _Cd;
    _Dx /= _Cd;
    // check upper boundary
    long long _Max = static_cast<long long>( std::numeric_limits<unsigned_integer_type>::max() );
    if ( _Nx > _Max || _Dx > _Max ) {
      throw std::overflow_error("calculation::rationalX<32,false, true>"
        "::rationalX(long long, long long)");
    }
    // encode
    this->numerator = static_cast<unsigned_integer_type>(_Nx);
    this->denominator = static_cast<unsigned_integer_type>(_Dx);
  }
  template<typename Number> explicit operator Number() const {
    return static_cast<Number>(numerator) / static_cast<Number>(denominator);
  }

  rationalX& operator=(const rationalX&) = default;
  rationalX& operator=(unsigned_integer_type _Nx) noexcept {
    numerator = _Nx;
    denominator = 1;
    return *this;
  }

  bool operator==(rationalX right) const {
    return this->numerator == right.numerator
      && this->denominator == right.denominator;
  }
  bool operator!=(rationalX right) const {
    return !((*this) == right);
  }
  bool operator<(rationalX right) const {
    return less_unsigned_rational(
      this->numerator, this->denominator,
      right.numerator, right.denominator,
      std::numeric_limits<unsigned_integer_type>::max());
  }
  bool operator>(rationalX right) const {
    return right < (*this);
  }
  bool operator>=(rationalX right) const {
    return !(*this < right);
  }
  bool operator<=(rationalX right) const {
    return !(*this > right);
  }

  rationalX& operator+=(rationalX right) {
    add_unsigned_rational(
      this->numerator, this->denominator,
      right.numerator, right.denominator,
      std::numeric_limits<unsigned_integer_type>::max()
      );
    return *this;
  }
  rationalX& operator-=(rationalX right) {
    sub_unsigned_rational(
      this->numerator, this->denominator,
      right.numerator, right.denominator,
      std::numeric_limits<unsigned_integer_type>::max()
      );
    return *this;
  }
  rationalX& operator*=(rationalX right) {
    mul_unsigned_rational(
      this->numerator, this->denominator,
      right.numerator, right.denominator,
      std::numeric_limits<unsigned_integer_type>::max()
      );
    return *this;
  }
  rationalX& operator/=(rationalX right) {
    return *this *= rationalX(right.denominator,right.numerator);
  }
  rationalX& operator++() {
    return *this += rationalX(denominator, unsigned_integer_type(1));
  }
  rationalX operator++(int) {
    rationalX copied = *this;
    ++(*this);
    return copied;
  }

  rationalX operator+(rationalX right) const {
    return rationalX(*this) += right;
  }
  rationalX operator-(rationalX right) const {
    return rationalX(*this) -= right;
  }
  rationalX operator*(rationalX right) const {
    return rationalX(*this) *= right;
  }
  rationalX operator/(rationalX right) const {
    return rationalX(*this) /= right;
  }
};

/** TEMPLATE SPECIAL 64bits unsigned rational */
template<>
struct rationalX<64,false, true> {
  typedef unsigned long long unsigned_integer_type;
  unsigned_integer_type numerator;
  unsigned_integer_type denominator;

  rationalX() = default;
  rationalX(unsigned int _Nx, unsigned int _Dx = 1) noexcept {
    unsigned int _Cd = gcd(_Nx, _Dx);
    this->numerator = _Nx/_Cd;
    this->denominator = _Dx/_Cd;
  }
  rationalX(unsigned long long _Nx, unsigned long long _Dx = 1){
    unsigned long long _Cd = gcd(_Nx, _Dx);
    this->numerator = _Nx/_Cd;
    this->denominator = _Dx/_Cd;
  }
  rationalX(int _Nx, int _Dx = 1) {
    if ( _Nx < 0 || _Dx < 0 ) {
      throw std::overflow_error("calculation::rationalX<64,false, true>"
        "::rationalX(int, int)");
    }
    int _Cd = gcd(_Nx, _Dx);
    this->numerator = _Nx / _Cd;
    this->denominator = _Dx / _Cd;
  }
  rationalX(long long _Nx, long long _Dx = 1) {
    if ( _Nx < 0 || _Dx < 0 ) {
      throw std::overflow_error("calculation::rationalX<64,false, true>"
        "::rationalX(long long, long long)");
    }
    long long _Cd = gcd(_Nx, _Dx);
    this->numerator = static_cast<unsigned_integer_type>( _Nx / _Cd );
    this->denominator = static_cast<unsigned_integer_type>( _Dx / _Cd );
  }
  template<typename Number> explicit operator Number() const {
    return static_cast<Number>(numerator) / static_cast<Number>(denominator);
  }

  rationalX& operator=(const rationalX&) = default;
  rationalX& operator=(unsigned_integer_type _Nx) noexcept {
    numerator = _Nx;
    denominator = 1;
    return *this;
  }

  bool operator==(rationalX right) const {
    return this->numerator == right.numerator
      && this->denominator == right.denominator;
  }
  bool operator!=(rationalX right) const {
    return !((*this) == right);
  }
  bool operator<(rationalX right) const {
    return less_unsigned_rational(
      this->numerator, this->denominator,
      right.numerator, right.denominator,
      std::numeric_limits<unsigned_integer_type>::max());
  }
  bool operator>(rationalX right) const {
    return right < (*this);
  }
  bool operator>=(rationalX right) const {
    return !(*this < right);
  }
  bool operator<=(rationalX right) const {
    return !(*this > right);
  }

  rationalX& operator+=(rationalX right) {
    add_unsigned_rational(
      this->numerator, this->denominator,
      right.numerator, right.denominator,
      std::numeric_limits<unsigned_integer_type>::max()
      );
    return *this;
  }
  rationalX& operator-=(rationalX right) {
    sub_unsigned_rational(
      this->numerator, this->denominator,
      right.numerator, right.denominator,
      std::numeric_limits<unsigned_integer_type>::max()
      );
    return *this;
  }
  rationalX& operator*=(rationalX right) {
    mul_unsigned_rational(
      this->numerator, this->denominator,
      right.numerator, right.denominator,
      std::numeric_limits<unsigned_integer_type>::max()
      );
    return *this;
  }
  rationalX& operator/=(rationalX right) {
    return *this *= rationalX(right.denominator,right.numerator);
  }
  rationalX& operator++() {
    return *this += rationalX(denominator, unsigned_integer_type(1));
  }
  rationalX operator++(int) {
    rationalX copied = *this;
    ++(*this);
    return copied;
  }

  rationalX operator+(rationalX right) const {
    return rationalX(*this) += right;
  }
  rationalX operator-(rationalX right) const {
    return rationalX(*this) -= right;
  }
  rationalX operator*(rationalX right) const {
    return rationalX(*this) *= right;
  }
  rationalX operator/(rationalX right) const {
    return rationalX(*this) /= right;
  }
};

/** TEMPLATE SPECIAL 32bits signed rational */
template<>
struct rationalX<32,true, true> {
  typedef int integer_type;
  typedef unsigned int unsigned_integer_type;
  integer_type numerator;
  integer_type denominator;

  rationalX() = default;
  rationalX(int _Nx, int _Dx = 1) noexcept {
    integer_type _Cd = gcd(_Nx, _Dx);
    this->numerator = _Nx/_Cd;
    this->denominator = _Dx/_Cd;
  }
  rationalX(long long _Nx, long long _Dx = 1){
    // reduce
    long long _Cd = gcd(_Nx, _Dx);
    _Nx /= _Cd;
    _Dx /= _Cd;
    //check lower boundary
    long long _Lowest = static_cast<long long>( std::numeric_limits<integer_type>::lowest() );
    if ( _Nx < _Lowest || _Dx < _Lowest ) {
      throw std::underflow_error("calculation::rationalX<32,true, true>"
        "::rationalX(long long, long long)");
    }
    // check upper boundary
    long long _Max = static_cast<long long>( std::numeric_limits<integer_type>::max() );
    if ( _Nx > _Max || _Dx > _Max ) {
      throw std::overflow_error("calculation::rationalX<32,true, true>"
        "::rationalX(long long, long long)");
    }
    // encode
    this->numerator = static_cast<integer_type>(_Nx);
    this->denominator = static_cast<integer_type>(_Dx);
  }
  rationalX(unsigned int _Nx, unsigned int _Dx = 1) {
    // reduce
    unsigned int _Cd = gcd(_Nx, _Dx);
    _Nx /= _Cd;
    _Dx /= _Cd;
    // check upper boundary
    unsigned int _Max = static_cast<unsigned int>( std::numeric_limits<integer_type>::max() );
    if ( _Nx > _Max || _Dx > _Max ) {
      throw std::overflow_error("calculation::rationalX<32,true, true>"
        "::rationalX(unsigned int, unsigned int)");
    }
    // encode
    this->numerator = static_cast<integer_type>(_Nx);
    this->denominator = static_cast<integer_type>(_Dx);
  }
  rationalX(unsigned long long _Nx, unsigned long long _Dx = 1) {
    // reduce
    unsigned long long _Cd = gcd(_Nx, _Dx);
    _Nx /= _Cd;
    _Dx /= _Cd;
    // check upper boundary
    unsigned long long _Max = static_cast<unsigned long long>( std::numeric_limits<integer_type>::max() );
    if ( _Nx > _Max || _Dx > _Max ) {
      throw std::overflow_error("calculation::rationalX<32,true, true>"
        "::rationalX(unsigned long long, unsigned long long)");
    }
    // encode
    this->numerator = static_cast<integer_type>(_Nx);
    this->denominator = static_cast<integer_type>(_Dx);
  }
  template<typename Number> explicit operator Number() const {
    return static_cast<Number>(numerator) / static_cast<Number>(denominator);
  }

  rationalX& operator=(const rationalX&) = default;
  rationalX& operator=(integer_type _Nx) {
    numerator = _Nx;
    denominator = 1; 
    return *this;
  }
		
  bool operator==(rationalX right) const {
    return this->numerator == right.numerator
      && this->denominator == right.denominator;
  }
  bool operator!=(rationalX right) const {
    return !((*this) == right);
  }
  bool operator<(rationalX right) const {
    constexpr bool _Pos = 0;
    constexpr bool _Neg = 1;
    bool this_sign = this->numerator < 0;
    bool right_sign = right.numerator < 0;
    if ( this_sign == right_sign ) {
      if ( this_sign == _Pos ) {
        // *this < right
        return less_unsigned_rational(
          static_cast<unsigned_integer_type>(abs(this->numerator)), static_cast<unsigned_integer_type>(this->denominator),
          static_cast<unsigned_integer_type>(abs(right.numerator)), static_cast<unsigned_integer_type>(right.denominator),
          std::numeric_limits<unsigned_integer_type>::max() );
      } else /* this_sign == _neg */ {
        // *this * -1 < right * -1 => *this > right
         return less_unsigned_rational(
          static_cast<unsigned_integer_type>(abs(right.numerator)), static_cast<unsigned_integer_type>(right.denominator),
          static_cast<unsigned_integer_type>(abs(this->numerator)), static_cast<unsigned_integer_type>(this->denominator),
          std::numeric_limits<unsigned_integer_type>::max() );
      }
    } else {

      return this_sign > right_sign;

    }
  }
  bool operator>(rationalX right) const {
    return right < (*this);
  }
  bool operator>=(rationalX right) const {
    return !(*this < right);
  }
  bool operator<=(rationalX right) const {
    return !(*this > right);
  }
		
  rationalX operator-() const noexcept {
    return rationalX(-numerator, denominator);
  }
  rationalX& operator+=(rationalX right) {
    auto this_numerator = static_cast<unsigned_integer_type>( abs(this->numerator) );
    auto this_denominator = static_cast<unsigned_integer_type>( this->denominator );
    auto right_numerator = static_cast<unsigned_integer_type>( abs(right.numerator) );
    auto right_denominator = static_cast<unsigned_integer_type>( right.denominator );
    bool this_sign = this->numerator < 0;
    bool right_sign = right.numerator < 0;

    constexpr bool _Neg = true;
    constexpr bool _Pos = false;
    if ( this_sign == right_sign ) {

      // this + right | -this + -right
      add_unsigned_rational(
        this_numerator, this_denominator,
        right_numerator, right_denominator,
        std::numeric_limits<unsigned_integer_type>::max() );
      // not change sign

    } else
    if ( this_sign == _Pos && right_sign == _Neg ) {

      if ( less_unsigned_rational(right_numerator,right_denominator,this_numerator,this_denominator,std::numeric_limits<unsigned_integer_type>::max()) ) {
        // this - right
        sub_unsigned_rational(
          this_numerator, this_denominator,
          right_numerator, right_denominator,
          std::numeric_limits<unsigned_integer_type>::max() );
        // not change sign
      } else {
        // -right + this
        sub_unsigned_rational(
          right_numerator, right_denominator,
          this_numerator, this_denominator,
          std::numeric_limits<unsigned_integer_type>::max() );
        this_sign = _Neg;
        this_numerator = right_numerator;
        this_denominator = right_denominator;
      }

    } 
    else /* this_sign == _Neg && right_sign == _Pos */ {

      if ( less_unsigned_rational(this_numerator,this_denominator,right_numerator,right_denominator,std::numeric_limits<unsigned_integer_type>::max()) ) {
        // right - this
        sub_unsigned_rational(
          right_numerator, right_denominator,
          this_numerator, this_denominator,
          std::numeric_limits<unsigned_integer_type>::max() );
        this_sign = _Pos;
        this_numerator = right_numerator;
        this_denominator = right_denominator;
      } else {
        // -this + right
        sub_unsigned_rational(
          this_numerator, this_denominator,
          right_numerator, right_denominator,
          std::numeric_limits<unsigned_integer_type>::max() );
        // not change sign
      }

    }

    const auto _Max_abs = static_cast<unsigned_integer_type>( std::numeric_limits<integer_type>::max() );
    const auto _Lowest_abs = static_cast<unsigned_integer_type>( abs(std::numeric_limits<integer_type>::lowest()) );
    if ( this_sign == _Pos ) {
      if ( this_numerator > _Max_abs || this_denominator > _Max_abs ) {
        throw std::overflow_error("ratioanX::operator+");
      }
      this->numerator = static_cast<integer_type>(this_numerator);
      this->denominator = static_cast<integer_type>(this_denominator);
    }
    else /* this_sign == _Neg */ {
      if ( this_numerator >= _Lowest_abs || this_denominator > _Lowest_abs ) {
        throw std::underflow_error("ratioanX::operator+");
      }
      this->numerator = - static_cast<integer_type>(this_numerator);
      this->denominator = static_cast<integer_type>(this_denominator);
    }

    return *this;
  }
  rationalX& operator-=(rationalX right) {
    return *this += (-right);
  }
  rationalX& operator*=(rationalX right) {
    auto this_numerator = static_cast<unsigned_integer_type>( abs(this->numerator) );
    auto this_denominator = static_cast<unsigned_integer_type>( this->denominator );
    auto right_numerator = static_cast<unsigned_integer_type>( abs(right.numerator) );
    auto right_denominator = static_cast<unsigned_integer_type>( right.denominator );

    mul_unsigned_rational(
      this_numerator, this_denominator,
      right_numerator, right_denominator,
      std::numeric_limits<unsigned_integer_type>::max() );
    bool this_sign = (this->numerator < 0) ^ (right.numerator < 0);

    constexpr bool _Neg = true;
    constexpr bool _Pos = false;
    const auto _Max_abs = static_cast<unsigned_integer_type>( std::numeric_limits<integer_type>::max() );
    const auto _Lowest_abs = static_cast<unsigned_integer_type>( abs(std::numeric_limits<integer_type>::lowest()) );
    if ( this_sign == _Pos ) {
      if ( this_numerator > _Max_abs || this_denominator > _Max_abs ) {
        throw std::overflow_error("ratioanX::operator*");
      }
      this->numerator = static_cast<integer_type>(this_numerator);
      this->denominator = static_cast<integer_type>(this_denominator);
    }
    else /* this_sign == _Neg */ {
      if ( this_numerator >= _Lowest_abs || this_denominator > _Lowest_abs ) {
        throw std::underflow_error("ratioanX::operator*");
      }
      this->numerator = - static_cast<integer_type>(this_numerator);
      this->denominator = static_cast<integer_type>(this_denominator);
    }

    return *this;
  }
  rationalX& operator/=(rationalX right) {
    return *this *= rationalX(right.denominator,right.numerator);
  }
  rationalX& operator++() {
    return *this += rationalX(denominator, integer_type(1));
  }
  rationalX operator++(int) {
    rationalX copied = *this;
    ++(*this);
    return copied;
  }

  rationalX operator+(rationalX right) const {
    return rationalX(*this) += right;
  }
  rationalX operator-(rationalX right) const {
    return rationalX(*this) -= right;
  }
  rationalX operator*(rationalX right) const {
    return rationalX(*this) *= right;
  }
  rationalX operator/(rationalX right) const {
    return rationalX(*this) /= right;
  }
};

/** TEMPLATE SPECIAL 64bits signed rational */
template<>
struct rationalX<64,true, true> {
  typedef long long integer_type;
  typedef unsigned long long unsigned_integer_type;
  integer_type numerator;
  integer_type denominator;

  rationalX() = default;
  rationalX(int _Nx, int _Dx = 1) noexcept {
    integer_type _Cd = gcd(_Nx, _Dx);
    this->numerator = _Nx/_Cd;
    this->denominator = _Dx/_Cd;
  }
  rationalX(long long _Nx, long long _Dx = 1){
    integer_type _Cd = gcd(_Nx, _Dx);
    this->numerator = _Nx / _Cd;
    this->denominator = _Dx / _Cd;
  }
  rationalX(unsigned int _Nx, unsigned int _Dx = 1) {
    unsigned int _Cd = gcd(_Nx, _Dx);
    this->numerator = static_cast<integer_type>(_Nx / _Cd);
    this->denominator = static_cast<integer_type>(_Dx / _Cd);
  }
  rationalX(unsigned long long _Nx, unsigned long long _Dx = 1) {
    // reduce
    unsigned long long _Cd = gcd(_Nx, _Dx);
    _Nx /= _Cd;
    _Dx /= _Cd;
    // check upper boundary
    unsigned long long _Max = static_cast<unsigned long long>( std::numeric_limits<integer_type>::max() );
    if ( _Nx > _Max || _Dx > _Max ) {
      throw std::overflow_error("calculation::rationalX<64,true, true>"
        "::rationalX(unsigned long long, unsigned long long)");
    }
    // encode
    this->numerator = static_cast<integer_type>(_Nx);
    this->denominator = static_cast<integer_type>(_Dx);
  }
  template<typename Number> explicit operator Number() const {
    return static_cast<Number>(numerator) / static_cast<Number>(denominator);
  }

  rationalX& operator=(const rationalX&) = default;
  rationalX& operator=(integer_type _Nx) {
    numerator = _Nx;
    denominator = 1; 
    return *this;
  }
		
  bool operator==(rationalX right) const {
    return this->numerator == right.numerator
      && this->denominator == right.denominator;
  }
  bool operator!=(rationalX right) const {
    return !((*this) == right);
  }
  bool operator<(rationalX right) const {
    constexpr bool _Pos = 0;
    constexpr bool _Neg = 1;
    bool this_sign = this->numerator < 0;
    bool right_sign = right.numerator < 0;
    if ( this_sign == right_sign ) {
      if ( this_sign == _Pos ) {
        // *this < right
        return less_unsigned_rational(
          static_cast<unsigned_integer_type>(abs(this->numerator)), static_cast<unsigned_integer_type>(this->denominator),
          static_cast<unsigned_integer_type>(abs(right.numerator)), static_cast<unsigned_integer_type>(right.denominator),
          std::numeric_limits<unsigned_integer_type>::max() );
      } else /* this_sign == _neg */ {
        // *this * -1 < right * -1 => *this > right
         return less_unsigned_rational(
          static_cast<unsigned_integer_type>(abs(right.numerator)), static_cast<unsigned_integer_type>(right.denominator),
          static_cast<unsigned_integer_type>(abs(this->numerator)), static_cast<unsigned_integer_type>(this->denominator),
          std::numeric_limits<unsigned_integer_type>::max() );
      }
    } else {

      return this_sign > right_sign;

    }
  }
  bool operator>(rationalX right) const {
    return right < (*this);
  }
  bool operator>=(rationalX right) const {
    return !(*this < right);
  }
  bool operator<=(rationalX right) const {
    return !(*this > right);
  }
		
  rationalX operator-() const noexcept {
    return rationalX(-numerator, denominator);
  }
  rationalX& operator+=(rationalX right) {
    auto this_numerator = static_cast<unsigned_integer_type>( abs(this->numerator) );
    auto this_denominator = static_cast<unsigned_integer_type>( this->denominator );
    auto right_numerator = static_cast<unsigned_integer_type>( abs(right.numerator) );
    auto right_denominator = static_cast<unsigned_integer_type>( right.denominator );
    bool this_sign = this->numerator < 0;
    bool right_sign = right.numerator < 0;

    constexpr bool _Neg = true;
    constexpr bool _Pos = false;
    if ( this_sign == right_sign ) {

      // this + right | -this + -right
      add_unsigned_rational(
        this_numerator, this_denominator,
        right_numerator, right_denominator,
        std::numeric_limits<unsigned_integer_type>::max() );
      // not change sign

    } else
    if ( this_sign == _Pos && right_sign == _Neg ) {

      if ( less_unsigned_rational(right_numerator,right_denominator,this_numerator,this_denominator,std::numeric_limits<unsigned_integer_type>::max()) ) {
        // this - right
        sub_unsigned_rational(
          this_numerator, this_denominator,
          right_numerator, right_denominator,
          std::numeric_limits<unsigned_integer_type>::max() );
        // not change sign
      } else {
        // -right + this
        sub_unsigned_rational(
          right_numerator, right_denominator,
          this_numerator, this_denominator,
          std::numeric_limits<unsigned_integer_type>::max() );
        this_sign = _Neg;
        this_numerator = right_numerator;
        this_denominator = right_denominator;
      }

    } 
    else /* this_sign == _Neg && right_sign == _Pos */ {

      if ( less_unsigned_rational(this_numerator,this_denominator,right_numerator,right_denominator,std::numeric_limits<unsigned_integer_type>::max()) ) {
        // right - this
        sub_unsigned_rational(
          right_numerator, right_denominator,
          this_numerator, this_denominator,
          std::numeric_limits<unsigned_integer_type>::max() );
        this_sign = _Pos;
        this_numerator = right_numerator;
        this_denominator = right_denominator;
      } else {
        // -this + right
        sub_unsigned_rational(
          this_numerator, this_denominator,
          right_numerator, right_denominator,
          std::numeric_limits<unsigned_integer_type>::max() );
        // not change sign
      }

    }

    const auto _Max_abs = static_cast<unsigned_integer_type>( std::numeric_limits<integer_type>::max() );
    const auto _Lowest_abs = static_cast<unsigned_integer_type>( abs(std::numeric_limits<integer_type>::lowest()) );
    if ( this_sign == _Pos ) {
      if ( this_numerator > _Max_abs || this_denominator > _Max_abs ) {
        throw std::overflow_error("ratioanX::operator+");
      }
      this->numerator = static_cast<integer_type>(this_numerator);
      this->denominator = static_cast<integer_type>(this_denominator);
    }
    else /* this_sign == _Neg */ {
      if ( this_numerator >= _Lowest_abs || this_denominator > _Lowest_abs ) {
        throw std::underflow_error("ratioanX::operator+");
      }
      this->numerator = - static_cast<integer_type>(this_numerator);
      this->denominator = static_cast<integer_type>(this_denominator);
    }

    return *this;
  }
  rationalX& operator-=(rationalX right) {
    return *this += (-right);
  }
  rationalX& operator*=(rationalX right) {
    auto this_numerator = static_cast<unsigned_integer_type>( abs(this->numerator) );
    auto this_denominator = static_cast<unsigned_integer_type>( this->denominator );
    auto right_numerator = static_cast<unsigned_integer_type>( abs(right.numerator) );
    auto right_denominator = static_cast<unsigned_integer_type>( right.denominator );

    mul_unsigned_rational(
      this_numerator, this_denominator,
      right_numerator, right_denominator,
      std::numeric_limits<unsigned_integer_type>::max() );
    bool this_sign = (this->numerator < 0) ^ (right.numerator < 0);

    constexpr bool _Neg = true;
    constexpr bool _Pos = false;
    const auto _Max_abs = static_cast<unsigned_integer_type>( std::numeric_limits<integer_type>::max() );
    const auto _Lowest_abs = static_cast<unsigned_integer_type>( abs(std::numeric_limits<integer_type>::lowest()) );
    if ( this_sign == _Pos ) {
      if ( this_numerator > _Max_abs || this_denominator > _Max_abs ) {
        throw std::overflow_error("ratioanX::operator*");
      }
      this->numerator = static_cast<integer_type>(this_numerator);
      this->denominator = static_cast<integer_type>(this_denominator);
    }
    else /* this_sign == _Neg */ {
      if ( this_numerator >= _Lowest_abs || this_denominator > _Lowest_abs ) {
        throw std::underflow_error("ratioanX::operator*");
      }
      this->numerator = - static_cast<integer_type>(this_numerator);
      this->denominator = static_cast<integer_type>(this_denominator);
    }

    return *this;
  }
  rationalX& operator/=(rationalX right) {
    return *this *= rationalX(right.denominator,right.numerator);
  }
  rationalX& operator++() {
    return *this += rationalX(denominator, integer_type(1));
  }
  rationalX operator++(int) {
    rationalX copied = *this;
    ++(*this);
    return copied;
  }

  rationalX operator+(rationalX right) const {
    return rationalX(*this) += right;
  }
  rationalX operator-(rationalX right) const {
    return rationalX(*this) -= right;
  }
  rationalX operator*(rationalX right) const {
    return rationalX(*this) *= right;
  }
  rationalX operator/(rationalX right) const {
    return rationalX(*this) /= right;
  }

  //  explicit 
  rationalX(float value) {
    /** 
     * significant-bitset move to numerator, 
     * exponential-bitset move to denominator, 
     * then reduce.
     *
     * 0 01000001 10101101011010110101000
     * 
     * = (1 + 0.10101101011010110101000) * 2^(01000001 - Bias) * (-1)^0  :IEEE754-floating
     * 
     * =  1.10101101011010110101000 * pow(2, exp2)
     *     
     * =  11010110101101011010100 * pow(2, exp2-significant_bits)
     *
     *    11010110101101011010100 << (exp2-significant_bits)
     * = ----------------------------------------------------    :Number*pow(2,X) = Number<<X
     *                       1
     *
     *    11010110101101011010100
     * = ---------------------------------------------------
     *                       1 >> (exp2-significant_bits)
    */
    constexpr unsigned int significant_mask
      = 0b00000000011111111111111111111111;
    constexpr unsigned int exponent_mask
      = 0b01111111100000000000000000000000;
    constexpr unsigned int sign_mask 
      = 0b10000000000000000000000000000000;
    constexpr unsigned int hidden_significant_mask
      = 0b00000000100000000000000000000000;
    constexpr unsigned int zero
      = 0b00000000000000000000000000000000;
    constexpr char exp2_bias = 127;
    constexpr char significant_bits = 23;

    unsigned int value_bitset = reinterpret_cast<uint32_t&>(value);
    if ( (value_bitset & exponent_mask) == exponent_mask && (value_bitset & significant_mask) != zero ) {
      throw std::exception("nan");
    }
    if ( value_bitset == zero ) {
      this->numerator = this->denominator = 0;
      return;
    }

    unsigned int exp2_bitset = (value_bitset & exponent_mask) >> significant_bits;
    unsigned int signifi_bitset = value_bitset & significant_mask | hidden_significant_mask;
    char exp2 = reinterpret_cast<char&>(exp2_bitset) - exp2_bias;
	
    // significant-bitset move to numerator, exponential-bitset move to denominator
    this->numerator = static_cast<integer_type>(signifi_bitset);
    this->denominator = static_cast<integer_type>(1);
    char dn_rshift = exp2 - significant_bits;
    if (dn_rshift > 0) {
      this->numerator <<= dn_rshift;
    } else if (dn_rshift < 0) {
      this->denominator <<= -dn_rshift;
    }

    // then reduce
    integer_type common_divisor = std::gcd(numerator, denominator);
    numerator /= common_divisor;
    denominator /= common_divisor;

    // sign
    if ((value_bitset & sign_mask) != 0) {
      this->numerator = -this->numerator;
    }
  }

  explicit 
  rationalX(double value) {
    static_assert(sizeof(integer_type) >= sizeof(unsigned long long), "ratianal(double)");
    constexpr unsigned long long significant_mask
      = 0b0000000000001111111111111111111111111111111111111111111111111111;
    constexpr unsigned long long exponent_mask
      = 0b0111111111110000000000000000000000000000000000000000000000000000;
    constexpr unsigned long long sign_mask     
      = 0b1000000000000000000000000000000000000000000000000000000000000000;
    constexpr unsigned long long hidden_significant_mask
      = 0b0000000000010000000000000000000000000000000000000000000000000000;
    constexpr unsigned int zero
      = 0b0000000000000000000000000000000000000000000000000000000000000000;
    constexpr short exp2_bias = 1023;
    constexpr char significant_bits = 52;

    unsigned long long value_bitset = reinterpret_cast<unsigned long long&>(value);
    if ( (value_bitset & exponent_mask) == exponent_mask && (value_bitset & significant_mask) != zero ) {
      throw std::exception("nan");
    }
    if ( value_bitset == zero ) {
      this->numerator = this->denominator = 0;
      return;
    }

    unsigned long long exp2_bitset = (value_bitset & exponent_mask) >> significant_bits;
    unsigned long long signifi_bitset = value_bitset & significant_mask | hidden_significant_mask;
    short exp2 = reinterpret_cast<short&>(exp2_bitset) - exp2_bias;
			
    // significant-bitset move to numerator, exponential-bitset move to denominator
    this->numerator = static_cast<integer_type>(signifi_bitset);
    this->denominator = static_cast<integer_type>(1);
    short dn_rshift = exp2 - significant_bits;
    if (dn_rshift > 0) {
      this->numerator <<= dn_rshift;
    } else if (dn_rshift < 0) {
      this->denominator <<= -dn_rshift;
    }

    // then reduce
    integer_type common_divisor = std::gcd(numerator, denominator);
    numerator /= common_divisor;
    denominator /= common_divisor;

    // sign
    if ((value_bitset & sign_mask) != 0) {
      this->numerator = -this->numerator;
    }
}
};


#define _calculation_rational_operator_with_literal(_OP_, _LITERAL_TYPE_)                   \
template<size_t b, bool s, bool opt> inline                                                 \
rationalX<b,s,opt> operator##_OP_##(const rationalX<b,s,opt>& left, _LITERAL_TYPE_ right) { \
  return left _OP_ rationalX<b,s,opt>(right, static_cast<_LITERAL_TYPE_>(1));               \
}

_calculation_rational_operator_with_literal(+, int)
_calculation_rational_operator_with_literal(+, long long)
_calculation_rational_operator_with_literal(+, unsigned int)
_calculation_rational_operator_with_literal(+, unsigned long long)

#undef _calculation_rational_operator_with_literal


using rational32_t = rationalX<32,true, true>;
using urational32_t = rationalX<32,false, true>;
using rational64_t = rationalX<64,true, true>;
using urational64_t = rationalX<64,false, true>;

template<size_t b, bool sb, bool opt> inline
rationalX<b,sb,opt> abs(const rationalX<b,sb,opt>& x) {
  return rationalX<b,sb,opt>( abs(x.numerator), x.denominator );
}

template<typename _Elem, typename _Traits, 
  size_t b, bool sb, bool opt> inline
std::basic_ostream<_Elem,_Traits>& operator<<(std::basic_ostream<_Elem,_Traits>& _Ostr, rationalX<b,sb,opt> _R) {
    return _Ostr << '(' << _R.numerator << ',' << _R.denominator << ')';
}

}// namespace calculation