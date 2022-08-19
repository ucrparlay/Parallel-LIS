#ifndef __PRIMITIVE_HPP_
#define __PRIMITIVE_HPP_

#include <cstdint>
#include <iterator>
#include <type_traits>
#include <utility>
#include <memory>
#include <typeinfo>
#include "parallelize.h"

// extern uint32_t threshold_scan;

static inline int log2i(const uint32_t x)
{
	return sizeof(x)*8-1-__builtin_clz(x);
}

template<typename T>
void scan_exclusive_nonrecursive(T *a, uint32_t n)
{
	const int depth = log2i(n);
	const int threshold = 1<<20;
	int d;
	// Up-sweep
	for(d=0; d<depth; ++d)
	{
		const uint32_t stride = 1u<<d;
		const auto end = n&~(2*stride-1);
		a[(n&~(stride-1))-1] = a[end-1];
		if(end/stride<threshold)
		{
			T t=0, sum=0;
			for(uint32_t k=0; k<end; k+=stride)
			{
				t = a[k+stride-1];
				a[k+stride-1] = sum;
				sum += t;
			}
			break;
		}
		else par_for(uint32_t k=0; k<end; k+=stride*2)
			a[k+stride*2-1] += a[k+stride-1];
	}
	// Down-sweep
	for(d--; d>=0; --d)
	{
		const uint32_t stride = 1u<<d;
		const auto end = n&~(2*stride-1);

		par_for(uint32_t k=0; k<end; k+=stride*2)
		{
			const auto t = a[k+stride-1];
			a[k+stride-1] = a[k+stride*2-1];
			a[k+stride*2-1] += t;
		}
		a[(n&~(stride-1))-1] += a[end-1]*((n>>d)&1);
	}
}

template<typename T>
static void scan_exclusive_impl(T* In, T* Out, T* B, T* C, uint32_t n)
{

	if (n==0) return;
	if (n == 1) {
		Out[0] = 0;
		return;
	}


	if(n<(1<<7))
	{
		Out[0] = 0;
		for(uint32_t i=1; i<n; ++i)
			Out[i] = Out[i-1]+In[i-1];
		return;
	}

	par_for(uint32_t i = 0; i < n/2; i++) 
		B[i] = In[2*i] + In[2*i+1];

	scan_exclusive_impl(B, C, B+n/2, C+n/2, n/2);
	
	par_for(uint32_t i=0; i<n/2*2; i++) {
		if(i&1)
		{
			Out[i] = In[i-1]+C[i/2];
		}
		else
		{
			Out[i] = C[i/2];
		}
	}
	if(n&1) Out[n-1] = Out[n-2]+In[n-2];
/*
	fprintf(stderr, "In:\t");
	for(uint32_t i=0; i<n; ++i) fprintf(stderr, "%d ", In[i]);
	fprintf(stderr, "\n");
	fprintf(stderr, "C:\t");
	for(uint32_t i=0; i<n/2; ++i) fprintf(stderr, "%d ", C[i]);
	fprintf(stderr, "\n");
	fprintf(stderr, "Out:\t");
	for(uint32_t i=0; i<n; ++i) fprintf(stderr, "%d ", Out[i]);
	fprintf(stderr, "\n");
*/
}

template<typename T>
static void scan_exclusive_impl_lite(T* A, uint32_t n, uint32_t stride)
{
/*
	if(n==0) return;
	if(n==1)
	{
		A[stride-1] = 0;
		return;
	}
*/

	if(n<(1<<20))
	{
		T t=0, sum=0;
		for(uint32_t i=0; i<n*stride; i+=stride)
		{
			t = A[i+stride-1];
			A[i+stride-1] = sum;
			sum += t;
		}
		return;
	}

	const auto end = (n/2)*(2*stride);
	A[n*stride-1] = A[end-1];
	par_for(uint32_t i=0; i<end; i+=stride*2)
		A[i+stride*2-1] += A[i+stride-1];

	scan_exclusive_impl_lite(A, n/2, stride*2);

	par_for(uint32_t i=0; i<end; i+=stride*2)
	{
		const auto t = A[i+stride-1];
		A[i+stride-1] = A[i+stride*2-1];
		A[i+stride*2-1] += t;
	}

	A[n*stride-1] += A[end-1]*(n&1);
}

template<typename T>
inline void scan_exclusive(T *a, uint32_t n)
{

	auto t1 = new T[n];
	auto t2 = new T[n];
	auto output = new T[n];
	scan_exclusive_impl(a, output, t1, t2, n);
	par_for(uint32_t i=0; i<n; ++i) a[i] = output[i];
	delete []output;
	delete []t2;
	delete []t1;

	// scan_exclusive_impl_lite(a, n, 1);
	// scan_exclusive_nonrecursive(a, n);
}

namespace detail{

template<typename InputIter, typename OutputIter, typename T, typename Index, class BinOp>
void scan_inclusive_impl(InputIter i_first, OutputIter d_first, T *si, T *sd, Index n, const BinOp &op)
{
	if (n==0) return;
	if (n==1)
	{
		*d_first = *i_first;
		return;
	}

	par_for(int i=0; i<n/2; ++i)
	{
		auto it = i_first+2*i;
		si[i] = op(*it, *(it+1));
	}

	scan_inclusive_impl(si, sd, si+n/2, sd+n/2, n/2, op);
	*d_first = *i_first;
	
	par_for(int i=1; i<n; ++i)
	{
		if (i&1) *(d_first+i) = sd[i/2];
		else *(d_first+i) = sd[i/2-1] + *(i_first+i);
	}
}

}

template<class Allocator=void, typename InputIter, typename OutputIter, class BinOp>
void scan_inclusive(InputIter first, InputIter last, OutputIter d_first, BinOp op);

template<class Allocator, typename InputIter, typename OutputIter, class BinOp>
void scan_inclusive(InputIter first, InputIter last, OutputIter d_first, BinOp op)
{
	static_assert(std::is_same_v<
		typename std::iterator_traits<InputIter>::iterator_category,
		std::random_access_iterator_tag
	>);
	static_assert(std::is_same_v<
		typename std::iterator_traits<OutputIter>::iterator_category,
		std::random_access_iterator_tag
	>);

	const auto n = std::distance(first, last);

	std::conditional_t<
		std::is_void_v<Allocator>,
		std::allocator<typename std::iterator_traits<OutputIter>::value_type>,
		Allocator
	> allocator;
	auto si_first = allocator.allocate(n);
	auto sd_first = allocator.allocate(n);
 
	detail::scan_inclusive_impl(first, d_first, si_first, sd_first, n, op);

	allocator.deallocate(sd_first, n);
	allocator.deallocate(si_first, n);
}

template<typename T, class UnaryPredicate>
uint32_t partition(T *a, uint32_t n, UnaryPredicate p)
{
	auto flag = new int[n];
	auto output = new T[n];

	par_for(uint32_t i=0; i<n; ++i)
		flag[i] = p(a[i]);
	scan_exclusive(flag, n);
	const uint32_t cnt_less = flag[n-1]+p(a[n-1]);
	par_for(uint32_t i=0; i<n; ++i)
		if(p(a[i])) output[flag[i]] = std::move(a[i]);
		else output[i-flag[i]+cnt_less] = std::move(a[i]);
	par_for(uint32_t i=0; i<n; ++i)
		a[i] = std::move(output[i]);

	delete []output;
	delete []flag;

	return cnt_less;
}

#endif // __PRIMITIVE_HPP_