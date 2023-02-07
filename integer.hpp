///
/// Infinite Precision Calculation for Dynamic Integer number
/// 
///@license Free, 2021-2022
///@author LongJiangnan, Jiang1998Nan@outlook.com
///@readme next work is calculate 'gamma(x)' by math::rational<math::dynamic_unsinged<>>
#pragma once
#ifndef _MATH_BEGIN
#define _MATH_BEGIN namespace math {
#define _MATH_END }
#endif

#include <cassert>
#include <vector>
#include <string>
#include <iosfwd>
#include <algorithm>// std::minmax()

_MATH_BEGIN
template<typename Word = unsigned char, 
	bool ReverseShift/* usually depends on platform|compiler */ = true>
class dynamic_unsigned {
public:
	//using Word = unsigned long long;
	static_assert((sizeof(Word) & (sizeof(Word) - 1)) == 0 
		&& sizeof(Word) <= sizeof(unsigned long long));

	/* Why list these quantities, because programmer have only 'bitwise operation' and 'incr|decr' before 'integer operation' . */
	static constexpr size_t bitsperword = ReverseShift ? (sizeof(Word) << 3) : (sizeof(Word) >> 3);/* sizeof(Word)*8 */
	static constexpr size_t log2_bitsperword = 
		bitsperword == 1 ? 0
		: bitsperword == 2 ? 1
			: bitsperword == 4 ? 2
				: bitsperword == 8 ? 3
					: bitsperword == 16 ? 4
						: bitsperword == 32 ? 5
							: /*bitsperword == 64 ? */6;
	static constexpr size_t modmask_bitsperword =
		bitsperword == 1 ? 0
		: bitsperword == 2 ? 1
			: bitsperword == 4 ? 3
				: bitsperword == 8 ? 7
					: bitsperword == 16 ? 15
						: bitsperword == 32 ? 31
							: /*bitsperword == 64 ? */63;
	static constexpr size_t bitsperword_mask = modmask_bitsperword;

	static size_t _Mulbitsperword(size_t dividend) {
		if constexpr (ReverseShift) {
			return dividend << log2_bitsperword;
		} else {
			return dividend >> log2_bitsperword;
		}
	}
	static size_t _Divbitsperword(size_t dividend) {
		if constexpr (ReverseShift) {
			return dividend >> log2_bitsperword;
		} else {
			return dividend << log2_bitsperword;
		}
	}
	static size_t _Modbitsperword(size_t dividend) {
		return dividend & modmask_bitsperword;
	}
	static Word _Leftshift(Word w, size_t bits) {
		if constexpr (ReverseShift) {
			return w >> bits;
		} else {
			return w << bits;
		}
	}
	static Word _Rightshift(Word w, size_t bits) {
		if constexpr (ReverseShift) {
			return w << bits;
		} else {
			return w >> bits;
		}
	}

	std::vector<Word> _My_wordset;
	size_t _My_bits;
	/// In order to correct operator~(), we introduce a state:
	/// 
	///		bool _I_am_inversed;
	/// 
	/// Very successful!
	/// But it also makes the state chaotic, and cannot but mark reasions in many places.
	/// So, we turn a idea, introduce a variable:
	Word _My_defaultword;
	/// Natural successful.

	static size_t _Correct_bits(typename std::vector<Word>::const_iterator _first, typename std::vector<Word>::const_iterator _last) {
		auto first = std::reverse_iterator<decltype(_last)>(_last);
		auto last = std::reverse_iterator<decltype(_first)>(_first);
		for ( ; first != last; ++first) {
			if (*first != static_cast<Word>(0)) {
				size_t bits = _Mulbitsperword(std::distance(first, last));
				for (size_t i = bitsperword - 1; i != -1; --i, --bits) {
					if ((*first) & _Rightshift(Word(1),i)) {
						break;
					}
				}
				/// first                              last
				///  |                                  |
				/// [11111111 11111111 11111111 11100000]
				///  |<-- -- -- -- 27 -- -- -- -->|
				return bits;
			}
		}

		return 0;
	}

	dynamic_unsigned()
		: _My_wordset(), _My_bits(0), _My_defaultword(static_cast<Word>(0)) {}

	explicit dynamic_unsigned(size_t bits, Word defaultword = static_cast<Word>(0))
		: _My_wordset(_Divbitsperword(bits) + (_Modbitsperword(bits)!=0?1:0), defaultword), _My_bits(bits), _My_defaultword(defaultword) {}

	dynamic_unsigned(std::vector<Word>&& wordset, size_t bits, Word defaultword = static_cast<Word>(0))
		: _My_wordset(std::move(wordset)), _My_bits(bits), _My_defaultword(defaultword) { 
		assert(_Mulbitsperword(_My_wordset.size()) >= _My_bits); 
		_My_bits = _Correct_bits(_My_wordset.begin(), 
									 std::next(_My_wordset.begin(), _Divbitsperword(_My_bits)+(_Modbitsperword(_My_bits)!=0?1:0)));
	}

	dynamic_unsigned(const std::string& str, bool reversed = true)
		: _My_wordset(), _My_bits(0), _My_defaultword(static_cast<Word>(0)) {
		if (!str.empty()) {
			if (str.size() >= 2 && str[0] == '0' && str[1] == 'b') {
				if (!reversed) {
					this->from_binstring(std::next(str.begin(),2), str.end());
				} else {
					this->from_binstring(str.rbegin(), std::prev(str.rend(),2));
				}
			} else if (str.size() >= 2 && str[0] == '0' && str[1] == 'x') {
				// ...
			} else {
				if (!reversed) {
					this->from_decstring(str.begin(), str.end());
				} else {
					this->from_decstring(str.rbegin(), str.rend());
				}
			}
		}
	}

	dynamic_unsigned(const char* c_str) : dynamic_unsigned(std::string(c_str), true) {}

	dynamic_unsigned(int i)
		: _My_wordset(sizeof(int)/sizeof(Word) + (sizeof(int)%sizeof(Word)!=0?1:0)), _My_bits(32), _My_defaultword(static_cast<Word>(0)) {
		std::memcpy(_My_wordset.data(), &i, sizeof(int));
		_My_bits = _Correct_bits(_My_wordset.begin(), _My_wordset.end());
	}

	size_t size() const {
		return _My_bits;
	}

	bool test(size_t pos) const {
		return _My_wordset[_Divbitsperword(pos)] & _Rightshift(static_cast<Word>(1),_Modbitsperword(pos));
	}
	
	void set(size_t pos) {
		_My_wordset[_Divbitsperword(pos)] |= _Rightshift(static_cast<Word>(1),_Modbitsperword(pos));
	}

	void clear() {
		_My_wordset = std::vector<Word>(0);
		_My_bits = 0;
		_My_defaultword = static_cast<Word>(0);
	}

	bool any() const {
		const Word None = static_cast<Word>(0);
		for (Word word : _My_wordset)
			if (word != None)
				return true;

		return false;
	}

	bool none() const {
		return !any();
	}

	bool all() const {
		// Test integral_part == All
		const Word All = (~static_cast<Word>(0));
		const size_t _My_intwords = _Divbitsperword(_My_bits);
		auto wordfirst = _My_wordset.begin();
		auto wordlast = std::next(_My_wordset.begin(), _My_intwords);
		for ( ; wordfirst != wordlast; ++wordfirst)
			if ((*wordfirst) != All)
				return false;

		// Test remaind_part[i] != None
		const Word None = static_cast<Word>(0);
		const size_t _My_rembits = _Modbitsperword(_My_bits);
		if (_My_rembits != 0) {
			for (size_t i = 0; i != _My_rembits; ++i)
				if (((*wordfirst) & _Rightshift(static_cast<Word>(1),i)) == None)
					return false;
		}

		return true;
	}

	explicit operator bool() const {
		return this->any() ? true : false;
	}

	bool operator!() const {
		return !static_cast<bool>(*this);
	}

	bool operator==(const dynamic_unsigned& right) const {
		// Compare tail part
		if (this->_My_defaultword != right._My_defaultword)
			return false;

		// Compare common part
		const size_t _Comm_words = std::min(_My_wordset.size(), right._My_wordset.size());
		auto first1 = this->_My_wordset.begin();
		auto last1  = std::next(this->_My_wordset.begin(), _Comm_words);
		auto first2 = right._My_wordset.begin();
		for ( ; first1 != last1; ++first1, ++first2)
			if ((*first1) != (*first2))
				return false;

		// Compare individual part are 'defaultword'
		if (this->_My_wordset.size() > _Comm_words) 
		{
			const Word _Comm_defaultword = _My_defaultword;
			last1 = this->_My_wordset.end();
			for ( ; first1 != last1; ++first1)
				if ((*first1) != _Comm_defaultword)
					return false;
		}
		else if (right._My_wordset.size() != _Comm_words) 
		{
			const Word _Comm_defaultword = _My_defaultword;
			auto last2 = right._My_wordset.end();
			for ( ; first2 != last2; ++first2)
				if ((*first2) != _Comm_defaultword)
					return false;
		}

		return true;
	}

	bool operator!=(const dynamic_unsigned& right) const {
		return !((*this) == right);
	}
	
	dynamic_unsigned& inc() {
		const size_t _My_rembits = _Modbitsperword(_My_bits);
		const size_t _My_words = _Divbitsperword(_My_bits) + (_My_rembits!=0?1:0);
		const Word None = static_cast<Word>(0);
		const Word All = ~static_cast<Word>(0);
		for (size_t i = 0; i != _My_words; ++i) {
			if (_My_wordset[i] == All) {
				_My_wordset[i] = None;
			} else {
				++_My_wordset[i];
				if (i+1 == _My_words) {
					if (_My_rembits != 0 && (_My_wordset[i] & _Rightshift(Word(1),_My_rembits))) {
						++_My_bits;
					}
				}
				return *this;
			}
		}

		assert(_My_rembits == 0);
		const size_t _My_wordunused = _My_wordset.size() - _My_words;
		if (_My_wordunused == 0) {
			_My_wordset.push_back(static_cast<Word>(1));
		} else {
			_My_wordset[_My_words+1] = static_cast<Word>(1);
		}
		++_My_bits;
		return *this;
	}

	dynamic_unsigned& dec() {
		const size_t _My_words = _Divbitsperword(_My_bits) + (_Modbitsperword(_My_bits)!=0?1:0);
		const Word None = static_cast<Word>(0);
		const Word All = ~static_cast<Word>(0);
		for (size_t i = 0; i != _My_words; ++i) {
			if (_My_wordset[i] == None) {
				_My_wordset[i] = All;
			} else {
				--_My_wordset[i];
				if (i+1 == _My_words) {
					if (_My_wordset[i] == None) {
						--_My_bits;
					} else if (!( _My_wordset[i] & _Rightshift(Word(1),_Modbitsperword(_My_bits - 1)) )) {
						--_My_bits;
					}
				}
				return *this;
			}
		}

		assert(_My_bits == 0);
		return *this;
	}

	dynamic_unsigned& operator<<=(size_t shift_bits) {
		if (shift_bits == 0) 
		{
			return *this;
		} 
		else if (shift_bits >= _My_bits) 
		{
			_My_wordset.resize(0);
			_My_bits = 0;
			return *this;
		} 
		else 
		{
			///          [11111111 11111111 11111111 11100000] smaller_shift 10
			/// [11111111 11111111 11111111 10000000]
			///           |--------+------| |---+--|
			///                    |            |
			/// [0,final)th from ( 'remainder pa|rt' of [10/8,final)th of source | 'shifted part' of [1,last)th of source )
			///                                 |
			///                          final word is 'remainder part' of final of source
			const size_t _My_words = _Divbitsperword(_My_bits) + (_Modbitsperword(_My_bits)!=0?1:0);
			size_t remove_words = _Divbitsperword(shift_bits);
			size_t remainder_bits = _Modbitsperword(shift_bits);
	#if 0
			size_t complement_bits = remainder_bits != 0 ? (bitsperword_mask&(~remainder_bits))+1 : 0;
	#else
			if (remainder_bits == 0) {/* Optimize for integral */
				std::move(
					std::next(_My_wordset.begin(), remove_words), 
					std::next(_My_wordset.begin(), _My_words),
					_My_wordset.begin());
				_My_bits -= shift_bits;
				return *this;
			}
			size_t complement_bits = (bitsperword_mask&(~remainder_bits))+1;
	#endif
			auto first = std::next(_My_wordset.begin(), remove_words);
			auto next = std::next(first);
			auto last = std::next(_My_wordset.begin(), _My_words);
			auto dest = _My_wordset.begin();
			assert(first != last);/* @see Case:boundary, first == last if and only if (_Divbitsperword(shift_bits) >= _My_words) => (shift_bits >= ceil(_My_bits,bitsperword)). */
			for ( ; next != last; ++first, ++next, ++dest)
				(*dest) = _Leftshift((*first),remainder_bits) | _Rightshift((*next),complement_bits);
			(*dest) = _Leftshift((*first),remainder_bits) | _My_defaultword;
			_My_bits -= shift_bits;
			return *this;
		}
	}

	dynamic_unsigned& operator>>=(size_t shift_bits) {
		if (shift_bits == 0 || _My_bits == 0) 
		{
			return *this;
		} 
		else 
		{
			/// [11111111 11111111 11111111 11100000] larger_shift 10
			/// [00000000 00111111 11111111 11111111 11111000]
			///  |---+--| |---+--| |-------------+----------|
			///      |        |                  |
			///  insert (10/8)|words             |
			///               |                  |
			///   (10/8)th word is 'shifted part'| of 0th of source 
			///                                  |
			///       [10/8+1,last)th from ( 'remainder pa|rt' of [0,final)th of source | 'shifted part' of [1,last)th of source )
			/// 
			/// For inplace and memory, first we insert to end"not begin", then reverse shift"inplace", final clear new words.
			const size_t _My_words = _Divbitsperword(_My_bits) + (_Modbitsperword(_My_bits)!=0?1:0);
			size_t _My_wordunused = _My_wordset.size() - _My_words;
			size_t insert_words = _Divbitsperword(shift_bits);
			size_t remainder_bits = _Modbitsperword(shift_bits);
	#if 0
			size_t complement_bits = remainder_bits != 0 ? (bitsperword_mask&(~remainder_bits))+1 : 0;
	#else
			if ( remainder_bits == 0 ) {/* Optimize:integral */
				if (_My_wordunused < insert_words) {
					_My_wordset.insert(_My_wordset.end(), insert_words - _My_wordunused, _My_defaultword);
					_My_wordunused = insert_words;
				}
				auto rdest = std::move( 
					std::next(_My_wordset.rbegin(),_My_wordunused), 
					_My_wordset.rend(), 
					std::next(_My_wordset.rbegin(),_My_wordunused - insert_words) 
				);
				std::fill(rdest, _My_wordset.rend(), _My_defaultword);
				_My_bits += shift_bits;
				return *this;
			}
			size_t complement_bits = (bitsperword_mask&(~remainder_bits))+1;
	#endif
			/* _Modbitsperword(_My_bits) + _Modbitsperword(shift_bits) > 8 => _Modbitsperword(_My_bits) > 8 - _Modbitsperword(shift_bits) */
			if (_Modbitsperword(_My_bits) > complement_bits || (_Modbitsperword(_My_bits) == 0 && remainder_bits != 0)) {
				if (_My_wordunused < insert_words+1) {
					_My_wordset.insert(_My_wordset.end(), insert_words+1 - _My_wordunused, _My_defaultword);
					_My_wordunused = insert_words;/* not is insert_words + 1, first 'shifted part' is 'None' if carry */
				} else {
					_My_wordunused -= 1;
				}
			} else {
				if (_My_wordunused < insert_words) {
					_My_wordset.insert(_My_wordset.end(), insert_words - _My_wordunused, _My_defaultword);
					_My_wordunused = insert_words;
				}
			}
			auto rfirst = std::next(_My_wordset.rbegin(), _My_wordunused);
			auto rlast = _My_wordset.rend();
			auto rnext = std::next(rfirst);
			auto rdest = std::next(_My_wordset.rbegin(), _My_wordunused - insert_words);
			for ( ; rnext != rlast; ++rfirst, ++rnext, ++rdest)
				(*rdest) = _Rightshift((*rfirst),remainder_bits) | _Leftshift((*rnext),complement_bits);
			(*rdest) = _Rightshift((*rfirst),remainder_bits) | _My_defaultword;
			std::fill(std::next(rdest), rlast, _My_defaultword);
			_My_bits += shift_bits;
			return *this;
		}
	}
	
	dynamic_unsigned operator<<(size_t shift_bits) const {
		return dynamic_unsigned(*this) <<= shift_bits;
	}

	dynamic_unsigned operator>>(size_t shift_bits) const {
		return dynamic_unsigned(*this) >>= shift_bits;
	}

	template<typename TransformOp> inline
	static std::vector<Word> transform(const std::vector<Word>& source, TransformOp func) {
		std::vector<Word> destination = std::vector<Word>(source.size(), static_cast<Word>(0));
		typename std::vector<Word>::const_iterator first = source.begin();
		typename std::vector<Word>::const_iterator last  = source.end();
		typename std::vector<Word>::iterator       dest  = destination.begin();
		for ( ; first != last; ++first, ++dest) {
			*dest = func(*first);
		}

		return destination;
	}

	template<typename TransformOp> inline
	static std::vector<Word>& transform(std::vector<Word>& source, TransformOp assignfunc) {
		typename std::vector<Word>::iterator first = source.begin();
		typename std::vector<Word>::iterator last  = source.end();
		for ( ; first != last; ++first) {
			assignfunc(*first);
		}

		return source;
	}

	dynamic_unsigned& inverse() {
		transform(_My_wordset, [](Word& a){ a = ~a; });
		_My_defaultword = ~_My_defaultword;
		return *this;
	}

	dynamic_unsigned operator~() const {
		return dynamic_unsigned(
			transform(_My_wordset, [](Word a) { return ~a; }),
			_My_bits,
			~_My_defaultword
		);
	}

	template<typename TransformOp> inline
	static std::vector<Word> transform(const std::vector<Word>& left, Word left_default, const std::vector<Word>& right, Word right_default, TransformOp func) {
		const std::pair<size_t, size_t> temp = std::minmax(left.size(), right.size());
		const size_t common_words = temp.first;
		const size_t required_words = temp.second;

		std::vector<Word> result = std::vector<Word>(required_words);
		typename std::vector<Word>::const_iterator first1 = left.begin();
		typename std::vector<Word>::const_iterator last1  = std::next(left.begin(), common_words);
		typename std::vector<Word>::const_iterator first2 = right.begin();
		typename std::vector<Word>::iterator       dest   = result.begin();
		for ( ; first1 != last1; ++first1, ++first2, ++dest) {
			(*dest) = func((*first1),(*first2));
		}

		if (left.size() > common_words) {
			last1 = left.end();
			for ( ; first1 != last1; ++first1, ++dest) {
				(*dest) = func((*first1),right_default);
			}
		} else if (right.size() > common_words) {
			auto last2 = right.end();
			for ( ; first2 != last2; ++first2, ++dest) {
				(*dest) = func(left_default,(*first2));
			}
		}

		return result;
	}

	template<typename TransformOp> inline
	static std::vector<Word>& transform(std::vector<Word>& left, Word left_default, const std::vector<Word>& right, Word right_default, TransformOp assignfunc) {
		const std::pair<size_t, size_t> temp = std::minmax(left.size(), right.size());
		const size_t common_words = temp.first;
		const size_t required_words = temp.second;

		if (left.size() != required_words) {
			/* left.append(required_words - left.size(), left_default); */
			left.resize(required_words, left_default);
		}
		typename std::vector<Word>::iterator       first1 = left.begin();
		typename std::vector<Word>::iterator       last1  = std::next(left.begin(), right.size());
		typename std::vector<Word>::const_iterator first2 = right.begin();
		for ( ; first1 != last1; ++first1, ++first2) {
			assignfunc((*first1),(*first2));
		}

		if (left.size() > right.size()) {
			last1 = left.end();
			for (; first1 != last1; ++first1) {
				assignfunc((*first1),right_default);
			}
		}

		return left;
	}

	dynamic_unsigned& operator&=(const dynamic_unsigned& right) {
		transform(this->_My_wordset, this->_My_defaultword, right._My_wordset, right._My_defaultword, [](Word& a, Word b){ a &= b; });
		size_t max_bits = std::max(this->_My_bits, right._My_bits);
		this->_My_bits = _Correct_bits(this->_My_wordset.begin(), 
												 std::next(this->_My_wordset.begin(), _Divbitsperword(max_bits)+(_Modbitsperword(max_bits)!=0?1:0)));
		this->_My_defaultword &= right._My_defaultword;
		return *this;
	}
	
	dynamic_unsigned operator&(const dynamic_unsigned& right) const {
		return dynamic_unsigned(
			transform(this->_My_wordset, this->_My_defaultword, right._My_wordset, right._My_defaultword, [](Word a, Word b){ return a&b; }),
			std::max(this->_My_bits, right._My_bits),
			this->_My_defaultword & right._My_defaultword
		);
	}

	dynamic_unsigned& operator|=(const dynamic_unsigned& right) {
		transform(this->_My_wordset, this->_My_defaultword, right._My_wordset, right._My_defaultword, [](Word& a, Word b){ a |= b; });
		size_t max_bits = std::max(this->_My_bits, right._My_bits);
		this->_My_bits = _Correct_bits(this->_My_wordset.begin(), 
												 std::next(this->_My_wordset.begin(), _Divbitsperword(max_bits)+(_Modbitsperword(max_bits)!=0?1:0)));
		this->_My_defaultword |= right._My_defaultword;
		return *this;
	}
	
	dynamic_unsigned operator|(const dynamic_unsigned& right) const {
		return dynamic_unsigned(
			transform(this->_My_wordset, this->_My_defaultword, right._My_wordset, right._My_defaultword, [](Word a, Word b){ return a|b; }),
			std::max(this->_My_bits, right._My_bits),
			this->_My_defaultword | right._My_defaultword
		);
	}

	dynamic_unsigned& operator^=(const dynamic_unsigned& right) {
		transform(this->_My_wordset, this->_My_defaultword, right._My_wordset, right._My_defaultword, [](Word& a, Word b){ a ^= b; });
		size_t max_bits = std::max(this->_My_bits, right._My_bits);
		this->_My_bits = _Correct_bits(this->_My_wordset.begin(), 
												 std::next(this->_My_wordset.begin(), _Divbitsperword(max_bits)+(_Modbitsperword(max_bits)!=0?1:0)));
		this->_My_defaultword ^= right._My_defaultword;
		return *this;
	}
	
	dynamic_unsigned operator^(const dynamic_unsigned& right) const {
		return dynamic_unsigned(
			transform(this->_My_wordset, this->_My_defaultword, right._My_wordset, right._My_defaultword, [](Word a, Word b){ return a^b; }),
			std::max(this->_My_bits, right._My_bits),
			this->_My_defaultword ^ right._My_defaultword
		);
	}

public:
	bool operator<(const dynamic_unsigned& right) const {
		if (this->_My_defaultword != right._My_defaultword) {
			return this->_My_defaultword < right._My_defaultword;
		} else {
			if (this->_My_bits != right._My_bits) {
				return this->_My_bits < right._My_bits;
			} else {
				const size_t common_words = _Divbitsperword(_My_bits) + (_Modbitsperword(_My_bits)!=0?1:0);
				const size_t this_wordunused = this->_My_wordset.size() - common_words;
				auto first1 = std::next(this->_My_wordset.rbegin(), this_wordunused);
				auto last1  = this->_My_wordset.rend();
				const size_t right_wordunused = right._My_wordset.size() - common_words;
				auto first2 = std::next(right._My_wordset.rbegin(), right_wordunused);
				for ( ; first1 != last1; ++first1, ++first2) {
					if ((*first1) != (*first2)) {
						return (*first1) < (*first2);
					}
				}

				return false;
			}
		}
	}

	bool operator>(const dynamic_unsigned& right) const {
		return right < (*this);
	}

	bool operator<=(const dynamic_unsigned& right) const {
		return !((*this) > right);
	}

	bool operator>=(const dynamic_unsigned& right) const {
		return !((*this) < right);
	}

	dynamic_unsigned& operator+=(dynamic_unsigned right) {
		dynamic_unsigned carry = (*this) & right;
		for ( ; carry.any(); carry = (*this) & right) {
			(*this) ^= right;
			right = carry >> static_cast<size_t>(1);
		}

		return ((*this) ^= right);
	}

	dynamic_unsigned& operator*=(dynamic_unsigned right) {
		dynamic_unsigned temp = *this;
		(*this) = dynamic_unsigned();
		
		while (right.any()) {
			if (right.test(0)) {
				(*this) += temp;
			}
			/* No overflow */
			temp >>= 1;
			right <<= 1;
		}

		return *this;
	}

	dynamic_unsigned& operator-=(dynamic_unsigned right) {
		for ( ; ; ) {
			(*this) ^= right;
			right &= (*this);
			if ( right.none() ) {
				break;
			}
			//#if _DEBUG
			//if ( right.test(_Bits-1) ) {// _Right > _Left
			//  throw std::underflow_error("std::operator-=(std::bitset<...>&, std::bitset<...>)");
			//}
			//#endif
			right >>= 1;
		}

		return *this;
	}

	dynamic_unsigned divide(dynamic_unsigned divisor) {
		const Word None = static_cast<Word>(0);
		if (this->_My_defaultword != None || divisor._My_defaultword != None) {
			throw std::overflow_error("dynamic_unsigned::divide(...)");
		}

		/**
		 * @figure 
		 *    001110          max bits is '6' = 8 - 3 + 1
		 *   ---------+
		 *    10110001| 101   1. 10110001<<5 = 001 < 101,  2. 10110001<<4 = 0001 >= 101, 10110001-(101>>4) = 111101
		 *        101    
		 *   ---------
		 *    101111          1. 101111<<3 = 111 >= 101, 101111-(101>>3) = 10101
		 *       101
		 *   --------
		 *    10101           1. 10101<<2 = 101 >= 101, 10101-(101>>2) = 1
		 *      101
		 *   -------
		 *    1               1. 1 < 101, so break
		*/
		if (this->_My_bits < divisor._My_bits) {/* special case '0' */
			return dynamic_unsigned();
		}
		dynamic_unsigned quotient = dynamic_unsigned(this->_My_bits - divisor._My_bits + 1);
		size_t i = quotient._My_bits - 1;
		divisor >>= i;
		for ( ; ; ) {
			if ( (*this) >= divisor ) {
				quotient.set(i);
				(*this) -= divisor;
			}
			
			if (i == 0) {
				break;
			}

			--i;
			divisor <<= 1;
		}

		if (!(quotient.test(quotient.size() - 1))) {
			--quotient._My_bits;
		}

		return quotient;
	}

	dynamic_unsigned& operator/=(dynamic_unsigned right) {
		dynamic_unsigned quotient = this->divide(right);
		std::swap(*this, quotient);
		return *this;
	}

	dynamic_unsigned& operator%=(dynamic_unsigned right) {
		this->divide(right);
		return *this;
	}

	std::string to_binstring(bool reversed = true) const {
		std::string str = std::string(this->size(), '\0');
		for (size_t i = 0; i != this->size(); ++i) {
			str[i] = this->test(i) ? '1' : '0';
		}

		if (reversed) {
			std::reverse(str.begin(), str.end());
		}

		if (str.empty()) {
			str = "0";
		}

		return str;
	}

	template<typename Iter>
	void from_binstring(Iter first, Iter last) {
		this->clear();
		
		_My_defaultword = static_cast<Word>(0);
		_My_bits = std::distance(first, last);
		const size_t _My_rembits = _Modbitsperword(_My_bits);
		const size_t _My_intwords = _Divbitsperword(_My_bits);
		const size_t _My_words = _My_intwords + (_My_rembits != 0 ? 1 : 0);
		_My_wordset.resize(_My_words, Word(0));

		auto wordfirst = _My_wordset.begin();
		auto wordlast = std::next(_My_wordset.begin(), _My_intwords);
		for ( ; wordfirst != wordlast; ++wordfirst) {
			Word& word = *wordfirst;
			/*Word eachbit = 1;
			while (eachbit != 0) {
				if (*bitfirst++ == '1') {
					word |= eachbit;
				}
				eachbit <<= 1;
			}*/
			for (size_t i = 0; i != bitsperword; ++i) {
				if (*first++ == '1') {
					word |= (Word(1) << i);
				}
			}
		}

		if (_My_rembits != 0) {
			Word& word = *wordfirst;
			for (size_t i = 0; i != _My_rembits; ++i) {
				if (*first++ == '1') {
					word |= (Word(1) << i);
				}
			}
		}

	}

	std::string to_octstring(bool reversed = true) const {
		std::string str = std::string(this->size()/3 + (this->size()%3!=0?1:0), '\0');
		dynamic_unsigned number = (*this);
		for (size_t i = 0; i != str.size(); ++i) {
			str[i] = ('0' + (number._My_wordset[0] & 8));
			number <<= 3;
		}

		if (reversed) {
			std::reverse(str.begin(), str.end());
		}

		return str;
	}

	std::string to_decstring(bool reversed = true) const {
		std::string str;
		dynamic_unsigned number = (*this);
		dynamic_unsigned base = "0b1010";
		while (number.size() != 0) {
			dynamic_unsigned temp = number.divide(base);
			std::swap(number, temp);
			str += ('0' + (temp._My_wordset.empty() ? 0 : static_cast<char>(temp._My_wordset[0])));
		}

		if (reversed) {
			std::reverse(str.begin(), str.end());
		}

		if (str.empty()) {
			str = "0";
		}

		return str;
	}
	
	template<typename Iter>
	void from_decstring(Iter _first, Iter _last) {
		this->clear();
		if (_first != _last) {
			dynamic_unsigned base = "0b1010";
			dynamic_unsigned remainder = dynamic_unsigned(4);
			auto first = std::reverse_iterator<Iter>(_last);
			auto last = std::reverse_iterator<Iter>(_first);
			for ( ; ; ) {
				remainder._My_wordset[0] = ((*first) - '0');
				(*this) += remainder;
				if (++first == last) {
					break;
				}
				(*this) *= base;
			}
		}
	}

	/// introduction
	///   A expression is form of 
	///     ans = a + (b + (c + d)) or ans = add(a, add(b, add(c,d)))
	///   Since operation order, we must first calculate 
	///     ans0 = c + d, then calculate 
	///     ans1 = b + ans0, final calculate
	///     ans  = a + ans1 .
	/// 
	///   Require three variables, little expensive, we hope a better method.
	///     ans1 = b + ans0 = ans0 + b
	///     ans  = a + ans1 = ans1 + a
	///    ----------------------------
	///   =>ans0 = c + d
	///     ans1 = ans0 + b
	///     ans  = ans1 + a
	///    ----------------------------
	///   =>ans = c + d
	///     ans += b
	///     ans += a
	///   Only one variable, very beauty.
	/// 
	///   But another expression
	///     ans = a + ((b * c) + (d * e))
	///    -------------------------------
	///   =>ans0 = b * c
	///     ans1 = d * e
	///     ans2 = ans0 + ans1
	///     ans  = a + ans2
	///    -------------------------------
	///   =>ans = b * c
	///     ans += ans1, Here we can't write the form of one variable, why?
	///     ans += a 
	/// 
	///   To illustrate why, first we introduce anothor form
	///                          ans           ans = a+
	///                          / \                  \
	///                         a + ans1              b+
	///     a + (b + (c + d)) =     /  \     =          \        What did we do ?
	///                            b + ans0             c+
	///                                / \                \
	///                               c + d                d
	///                          ans           ans = a+
	///                          / \                  \
	///                         a + ans2             ans2         
	///      a + (b*c + d*e)  =     /  \     =       /  \         What is difference?
	///                         ans0 + ans1         b*   d*
	///                          / \   / \           \    \
	///                         b * c d * e           c    e
	/// 
	/// expression-tree 
	///		We shrink left-leaf node to its parent, tree of 'a + (b + (c + d))' became a path(each degree of vertex is 1).
	///		Graph Thery, we can direct through a path(each vertex passes only once),
	///		So we need only one variable if and only if the expression is a path.

	//struct Expr {
	//  virtual size_t required_words() const { return 0; }
	//  virtual size_t required_bits() const { return 0; }
	//  virtual void get_result(std::vector<Word>& result) const {}
	//};

	//struct ExprOp : public std::vector<Expr*> {
	//  using _Mybase = std::vector<Expr*>;
	//  using _Mybase::_Mybase;
	//  ExprOp(ExprOp&& right) noexcept : _Mybase(std::move(right)) {}
	//  ~ExprOp() { for(Expr* expr : *this){ delete expr; } }
	//};

	//struct AndExpr0 : public Expr {
	//  const std::vector<Word>& left;
	//  size_t left_bits;
	//  const std::vector<Word>& right;
	//  size_t right_bits;

	//  AndExpr0(const std::vector<Word>& arg0, size_t arg1, const std::vector<Word>& arg2, size_t arg3) : left(arg0), left_bits(arg1), right(arg2), right_bits(arg3) {}

	//  virtual size_t required_words() const override {
	//    return std::max(left.size(), right.size());
	//  }
	//  
	//  virtual size_t required_bits() const override {
	//    return std::max(left_bits, right_bits);
	//  }

	//  virtual void get_result(std::vector<Word>& result) const override {
	//    // Compute this expr(left &= right)
	//    const size_t common_size = std::min(left.size(), right.size());
	//    auto first1 = left.begin();
	//    auto first2 = right.begin();
	//    auto dest   = result.begin();
	//    auto dlast  = std::next(result.begin(), common_size);
	//    for ( ; dest != dlast; ++first1, ++first2, ++dest)
	//      (*dest) = (*first1) & (*first2);

	//    // Compute this expr(left & 0)|expr(0 & right)
	//    /* assert( [dest,std::next(dlast,X.size-common_size)) == Word(0) ) */
	//  }
	//};

	//struct AndExpr1 : public Expr {
	//  const std::vector<Word>& left;
	//  size_t left_bits;
	//  const Expr& right;

	//  AndExpr1(const std::vector<Word>& arg0, size_t arg1, const Expr& arg2) : left(arg0), left_bits(arg1), right(arg2) {}

	//  virtual size_t required_words() const override {
	//    return std::max(left.size(), right.required_words());
	//  }

	//  virtual size_t required_bits() const override {
	//    return std::max(left_bits, right.required_bits());
	//  }

	//  virtual void get_result(std::vector<Word>& result) const override {
	//    const size_t right_size = right.required_words();

	//    // Compute right expr
	//    right.get_result(result);

	//    // Compute this expr(left &= right)
	//    const size_t common_size = std::min(left.size(), right_size);
	//    auto first1 = result.begin();
	//    auto last1  = std::next(result.begin(), common_size);
	//    auto first2 = left.begin();
	//    for ( ; first1 != last1; ++first1, ++first2)
	//      (*first1) &= (*first2);

	//    if (left.size() > common_size) {
	//      // Compute this expr(left &= 0)
	//      /* assert( [first1,std::next(last1,left.size()-common_size)) == Word(0) ) */
	//    } else if (right_size > common_size) {
	//      // Compute this expr(left = 0 & right)
	//      std::advance(last1, (right_size - common_size));
	//      for ( ; first1 != last1; ++first1)
	//        (*first1) = 0;
	//    }
	//  }
	//};

	//struct AndExpr2 : public Expr {
	//  Expr left;
	//  Expr right;

	//  virtual size_t required_words() const override {
	//    return std::max(left.required_words(), right.required_words());
	//  }
	//  
	//  virtual size_t required_bits() const override {
	//    return std::max(left.required_bits(), right.required_bits());
	//  }

	//  virtual void get_result(std::vector<Word>& result) const override {
	//    const size_t left_size = left.required_words();
	//    const size_t right_size = right.required_words();

	//    // Compute left expr
	//    left.get_result(result);

	//    // Compute right expr
	//    auto temp = std::vector<Word>(right_size, Word(0));
	//    right.get_result(temp);
	//    
	//    // Compute this expr(left &= right)
	//    const size_t common_size = std::min(left_size, right_size);
	//    auto first1 = result.begin();
	//    auto last1  = std::next(result.begin(), common_size);
	//    auto first2 = temp.begin();
	//    for ( ; first1 != last1; ++first1, ++first2)
	//      (*first1) &= (*first2);

	//    if (left_size > common_size) {
	//      // Compute this expr(left &= 0)
	//      std::advance(last1, (left_size - common_size));
	//      for (; first1 != last1; ++first1)
	//        (*first1) = 0;
	//    } else if (right_size > common_size) {
	//      // Compute this expr(left = 0 & right)
	//      std::advance(last1, (right_size - common_size));
	//      for ( ; first1 != last1; ++first1)
	//        (*first1) = 0;
	//    }
	//  }
	//};

	//dynamic_unsigned(ExprOp&& expr) {
	//  this->_My_word_set.resize(expr.back()->required_words(), static_cast<Word>(0));
	//  this->bits = expr.back()->required_bits();
	//  expr.back()->get_result(this->_My_word_set);
	//}

	/*ExprOp operator&(const dynamic_unsigned& right) const {
		return ExprOp{ new AndExpr0(this->_My_word_set, this->bits, right._My_word_set, right.bits) };
	}

	ExprOp operator&(ExprOp&& right) const {
		right.push_back(new AndExpr1(this->_My_word_set, this->bits, *right.back()));
		return std::move(right);
	}

	friend ExprOp operator&(ExprOp&& left, const dynamic_unsigned& right) {
		left.push_back(new AndExpr1(right._My_word_set, right.bits, *left.back()));
		return std::move(left);
	}*/
};// end of class dynamic_unsigned<..>

template<typename Word, bool ReverseShift>
std::string to_string(const dynamic_unsigned<Word, ReverseShift>& bitset) {
	return bitset.to_decstring();
}

template<typename Word, bool ReverseShift>
std::ostream& operator<<(std::ostream& ostr, const dynamic_unsigned<Word, ReverseShift>& bitset) {
	if (ostr.flags() & std::ios::oct) {
		return ostr << bitset.to_octstring();
	} else if(ostr.flags() & std::ios::binary) {
		return ostr << bitset.to_binstring();
	} else {
		return ostr << bitset.to_decstring();
	}
}
_MATH_END

/*
#include <iostream>
#include <clmagic/calculation/number/dynamicint.h>

int main() {
	srand(time(0));

	using namespace::calculation;

	std::cout.setf(std::ios::binary);
	{
		dynamic_unsigned<> a = "1101010111101110101111111111111111111100000000";
		std::cout << "source              is:" << a << std::endl;
		std::cout << "(~source)           is:" << (~a) << std::endl;
		std::cout << "(~source)&111111    is:" << ((~a)&(dynamic_unsigned<>("111111"))) << std::endl;
		std::cout << "(~source)|111111    is:" << ((~a)|(dynamic_unsigned<>("111111"))) << std::endl;
		std::cout << "(source)&(~111111)  is:" << ((a)&(~dynamic_unsigned<>("111111"))) << std::endl;
		std::cout << "(~source)&(~111111) is:" << ((~a)&(~dynamic_unsigned<>("111111"))) << std::endl;
		std::cout << "(source)|(~111111)  is:" << ((a)|(~dynamic_unsigned<>("111111"))) << std::endl;
		std::cout << "(~source)|(~111111) is:" << ((~a)|(~dynamic_unsigned<>("111111"))) << std::endl;
		std::cout << "(source)^(~111111)  is:" << ((a)^(~dynamic_unsigned<>("111111"))) << std::endl;
		std::cout << "(~source)^(~111111) is:" << ((~a)^(~dynamic_unsigned<>("111111"))) << std::endl;
	}
	std::cout << std::endl;
	std::cout << std::endl;
	std::cout << std::endl;
	{
		dynamic_unsigned<> a = "101010111101110101111111111111111111100000000";
		std::cout << "source is:\n\t";
		for (size_t i = 0; i != a.size(); ++i)
			std::cout << a.test(i);
		std::cout << std::endl;

		std::cout << "left shift start\n";
		for (size_t i = 0; i != a.size(); ++i) {
			std::cout << '\t';
			dynamic_unsigned<> b = a;
			std::cout << (b <<= i) << std::endl;
		}
	}
	{
		dynamic_unsigned<> a = "1111101110101111111111111111111100000000";
		a._My_wordset.resize(1000);
		std::cout << "source is:\n\t";
		for (size_t i = 0; i != a.size(); ++i)
			std::cout << a.test(i);
		std::cout << std::endl;

		std::cout << "right shift start\n";
		for (size_t i = 0; i != a.size(); ++i) {
			std::cout << '\t';
			dynamic_unsigned<> b = a;
			std::cout << (b >>= i) << std::endl;
		}
	}
	std::cout << std::endl;
	std::cout << std::endl;
	std::cout << std::endl;
	{
		dynamic_unsigned<> a = "0";
		std::cout << "inc:\n";
		for (size_t i = 0; i != 10000; ++i) {
			std::cout << a.inc() << std::endl;
		}
		std::cout << "dec:\n";
		for (size_t i = 0; i != 10000; ++i) {
			std::cout << a.dec() << std::endl;
		}
	}


	return 0;
}
*/