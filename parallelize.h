#ifndef __PARALLELIZE_H_
#define __PARALLELIZE_H_

#ifndef ENABLE_PARALLELISM
	#define ENABLE_PARALLELISM 1
#endif 

#if ENABLE_PARALLELISM
#include <cilk/cilk.h>
#include <cilk/cilk_api.h>
//#include <parlay/slice.h>
//#include <parlay/primitives.h>
#include <functional>
#include <sample_sort.h>
#include <merge.h>
#include <sequence_ops.h>
#else
// #include <algorithm>
#endif

#include <utility>

#if ENABLE_PARALLELISM
	#define _spawn cilk_spawn
	#define _sync cilk_sync
	#define par_for cilk_for
#else
	#define _spawn
	#define _sync
	#define par_for for
#endif

template<class Iter, class Comp=std::less<>>
void parallel_sort(Iter begin, Iter end, Comp comp=Comp())
{
#if ENABLE_PARALLELISM
	pbbs::sample_sort(begin, std::distance(begin,end), comp);
	// parlay::sort_inplace(parlay::make_slice(begin,end), std::move(comp));
#else
	std::sort(begin, end, comp);
#endif
}

#endif // __PARALLELIZE_H_