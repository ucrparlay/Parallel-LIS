#ifndef _UTIL_HPP_
#define _UTIL_HPP_

#include <cstdint>


static inline uint32_t log2u_floor(const uint32_t x)
{
	return sizeof(x)*8-1-(uint32_t)__builtin_clz(x);
}

static inline uint32_t log2u_ceil(const uint32_t x)
{
	const auto floor = log2u_floor(x);
	return (1u<<floor)==x? floor: floor+1;
}

static inline uint32_t div_ceil(const uint32_t x, const uint32_t y)
{
	return (x+y-1)/y;
}

static inline uint32_t bit_lowest(const uint32_t x)
{
	return (uint32_t)__builtin_ctz(x);
}

static inline uint32_t cnt_bit(const uint32_t x)
{
	return (uint32_t)__builtin_popcount(x);
}

struct completeBST
{
	completeBST(uint32_t n_) :
		n(n_),
		bitwise(log2u_floor(n)),
		base((n-(1u<<bitwise)+1)*2)
	{
	}

	completeBST(const completeBST &other) = default;
	completeBST& operator=(const completeBST &other) = default;

	uint32_t index2rank(uint32_t index) const
	{
		const auto depth = log2u_floor(index);
		const auto offset = index-(1u<<depth);
		auto rank = ((offset<<1)|1)<<(bitwise-depth);
		return rank>base? (rank+base)/2: rank;
	}

	uint32_t rank2index(uint32_t rank) const
	{
		if(rank>base)
			rank = rank*2 - base;

		const auto lowest_bit = bit_lowest(rank);
		const auto depth = bitwise-lowest_bit;
		const auto offset = (rank>>(lowest_bit+1));
		const auto index = offset+(1u<<depth);
		return index;
	}

	uint32_t get_size(uint32_t index) const
	{
		const auto depth = log2u_floor(index);
		const auto diff = bitwise-depth;
		uint32_t size = (1u<<diff)-1;
		const auto lindex = index<<diff;
		const auto rindex = lindex+size;
		if(lindex<=n)
		{
			if(rindex<=n) size += size+1;
			else size += n-lindex+1;
		}
		return size;
	}

	uint32_t size() const { return n; }

private:
	uint32_t n, bitwise, base;
};


#endif // _UTIL_HPP_