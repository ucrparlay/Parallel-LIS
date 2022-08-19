#ifndef __INTERVAL_HPP_
#define __INTERVAL_HPP_

#include <iostream>
#include <cstdint>
#include <cassert>
#include <algorithm>
#include <utility>
#include <memory>
#include <random>
#include <tuple>
#include <set>
#include <chrono>
#include <type_traits>
// #include <random_shuffle.h>
#include <sample_sort.h>
#include <parlay/utilities.h>
#include "primitive.hpp"
#include "util.hpp"
#include "interval_common.hpp"
// #include "lib/DappurCodebase/function.hpp"
using std::cerr;
using std::cout;
using std::endl;
using namespace std::chrono;
using namespace parlay;

template<typename Iter, typename T, class Comp>
auto merge_one(Iter begin, Iter end, T x, Comp comp)
{
	const auto n = std::distance(begin, end);
	const auto p = std::lower_bound(begin, end, x, comp)-begin;
	auto res = pbbs::sequence<T>::no_init(n+1);
	par_for(uint32_t i=0; i<p; ++i)
		res[i] = *(begin+i);
	res[p] = x;
	par_for(uint32_t i=p; i<n; ++i)
		res[i+1] = *(begin+i);
	return res;
}

void calc_usage(size_t *const usage, size_t begin, size_t end)
{
	size_t n = end-begin;
	if(n==0) return;
	if(n==1)
	{
		usage[begin] = n;
		return;
	}
	auto p = n/2;
	usage[begin+p] = n;
	_spawn calc_usage(usage, begin, begin+p);
	calc_usage(usage, begin+p+1, end);
}

template<typename T, template<typename> class U, template<typename> class Alloc=std::allocator>
struct interval
{
	typedef U<T> IU;

	struct treap_node{
		int priority;
		treap_node *lch, *rch;
		T key;
		typename IU::type_data userdata;
		uint32_t size;

		void update(){
			/*
			const auto userdata_tmp = IU::g(key, true);
			if constexpr(IU::is_nested_alloc)
				assert(userdata.allocator.buffer_end);
			
			if constexpr(has_replace_one<typename IU::type_data, decltype(userdata_tmp)>::value)
				userdata.replace_one(userdata_tmp);
			else userdata = userdata_tmp;
			*/
			userdata = IU::g(key);

			/*
			if constexpr(std::is_integral_v<decltype(userdata)>){
				assert(userdata!=0);
			}
			*/
			size = 1;
			if(lch) userdata = IU::f(lch->userdata, userdata), size += lch->size;
			if(rch) userdata = IU::f(userdata, rch->userdata), size += rch->size;
		}

		void update_size(){
			size = 1;
			if(lch) size += lch->size;
			if(rch) size += rch->size;
		}

		void update_recur()
		{
			if(lch) lch->update_recur();
			if(rch) rch->update_recur();
			update();
		}
	};
	typedef std::tuple<treap_node*,treap_node*,treap_node*> treap_component;
	treap_node *root;
	Alloc<treap_node> allocator;

	// typedef void (*func_update_t)(const T &element, const typename IU::type_data &userdata);

	// treap_node* build_treap(treap_node *a, uint32_t n);

	// typename IU::type_data query_impl(const treap_node *u, const T &begin, const T &end);
	// typename IU::type_data query_impl(const treap_node *u, const T &qb, const T &qe, const T &bb, const T &be);
	template<typename G, typename F>
	auto query_left(const treap_node *u, const T &end, const G &retrieve_info, const F &merge_info) const;
	template<typename G, typename F>
	auto query_right(const treap_node *u, const T &begin, const G &retrieve_info, const F &merge_info) const;

	template<typename G, typename F>
	auto query_left(const T &end, const G &retrieve_info, const F &merge_info) const
	{
		return query_left(root, end, retrieve_info, merge_info);
	}

	auto query_left(const T &end) const
	{
		return query_left(root, end, [](const typename IU::type_data &x){return x;}, IU::f);
	}

	template<typename G>
	auto query_pivot(const treap_node *u, const G &retrieve_info) const
	{
		while(true)
		{
			if(retrieve_info(IU::g(u->key,true))<0)
				return -long(u->key);
			if(u->rch && retrieve_info(u->rch->userdata)<0)
				u = u->rch;
			else u = u->lch;
		}
	}

	template<typename G, typename F>
	auto query_special(const treap_node *u, const T &end, const G &retrieve_info, const F &merge_info) const
	{
		if(u==nullptr) return retrieve_info(typename IU::type_data());

		const auto &ue = u->key;
		if(!IU::compare(ue,end)) // ue >= end
			return query_special(u->lch, end, retrieve_info, merge_info);

		const auto info_mid = retrieve_info(IU::g(u->key,true));
		auto res = merge_info(
			info_mid<0? -long(ue): info_mid,
			query_special(u->rch, end, retrieve_info, merge_info)
		);
		if(res<0) return res;
		if(u->lch)
		{
			const auto info_left = retrieve_info(u->lch->userdata);
			res = merge_info(info_left<0? query_pivot(u->lch, retrieve_info): info_left, res);
			if(res<0) return res;
		}
		return res;
	}

	template<typename G, typename F>
	auto query_special(const T &end, const G &retrieve_info, const F &merge_info) const
	{
		return query_special(root, end, retrieve_info, merge_info);
	}

	auto query_special(const T &end) const
	{
		return query_special(root, end, [](const typename IU::type_data &x){return x;}, IU::f);
	}

	template<typename G, typename F>
	auto query_impl(const treap_node *u, const T &begin, const T &end, const G &retrieve_info, const F &merge_info) const
	{
		if(u==nullptr) return retrieve_info(typename IU::type_data());

		const auto &ue = u->key;
		if(!IU::compare(begin,ue)) return query_impl(u->rch, begin, end, retrieve_info, merge_info);
		if(!IU::compare(ue,end)) return query_impl(u->lch, begin, end, retrieve_info, merge_info);
		const auto res_left = _spawn query_right(u->lch, begin, retrieve_info, merge_info);
		const auto res_right = query_left(u->rch, end, retrieve_info, merge_info);
		_sync;
		return merge_info(
			merge_info(
				res_left,
				retrieve_info(IU::g(u->key,true))
			),
			res_right
		);
	}
	template<typename G, typename F>
	auto query(const T &begin, const T &end, const G &retrieve_info, const F &merge_info) const
	{
		/*
		using trait_G = DPCB::function_trait<G>;
		using trait_F = DPCB::function_trait<F>;
		static_assert(trait_G::type_argument::size==1, "`merge_info` requires one argument");
		static_assert(trait_F::type_argument::size==2, "`merge_info` requires two arguments");
		static_assert(std::is_convertible_v<typename IU::type_data, typename trait_G::type_argument::template type<0>>);
		static_assert(std::is_same_v<typename trait_F::type_return, typename trait_G::type_return>);
		static_assert(std::is_convertible_v<typename trait_G::type_return, typename trait_F::type_argument::template type<0>>);
		static_assert(std::is_convertible_v<typename trait_G::type_return, typename trait_F::type_argument::template type<1>>);
		*/
		return query_impl(root, begin, end, retrieve_info, merge_info);
	}
	auto query(const T &begin, const T &end) const
	{
		return query_impl(root, begin, end, [](const typename IU::type_data &x){return x;}, IU::f);
	}

	typename IU::type_data query_all(){
		return root->userdata;
	}

	bool update_impl(treap_node *curr, const T &element, void (*const assign)(T&, const T&))
	{
		if(!curr) return false;

		bool retval = true;
		if(IU::compare(element, curr->key))
			retval = update_impl(curr->lch, element, assign);
		else if(IU::compare(curr->key, element))
			retval = update_impl(curr->rch, element, assign);
		else assign(curr->key, element);

		curr->update();
		return retval;
	}

	template<typename Iter>
	void update_impl(treap_node *u, Iter begin, Iter end, void (*const assign)(T&, const T&))
	{
		if(u==nullptr || begin==end) return;

		const auto n = std::distance(begin, end);
		auto &ue = u->key;
		// const auto lb = partition(begin, n, [&](const auto &x){return IU::compare(x,ue);});
		const auto lb = std::lower_bound(begin, end, ue, IU::compare)-begin;
		// auto less_val = [&](const T& x){return IU::compare(x,ue);};
		// pbbs::sequence<T> seq(begin, n);
		// const auto lb = pbbs::binary_search(seq, less_val);
/*
		auto seq = pbbs::sequence<uint64_t>(begin, n);
		auto fl = pbbs::delayed_seq<bool>(n, [&](const size_t i){return !IU::compare(seq[i],ue);});
		const auto [out, nlb] = pbbs::split_two(seq, fl);
		par_for(size_t i=0; i<nlb; ++i)
			assert(!IU::compare(out[i],ue)==false);
		par_for(size_t i=nlb; i<out.size(); ++i)
			assert(!IU::compare(out[i],ue)==true);
		assert(out.size()==n);
		par_for(size_t i=0; i<out.size(); ++i)
			*(begin+i) = out[i];
		
		if(nlb!=lb)
		{
			cerr << "nlb: " << nlb << '\n';
			cerr << "lb: " << lb << '\n';
		}
		assert(nlb==lb);
*/		
		const bool match_u = lb!=n&&!IU::compare(ue,*(begin+lb));
		const auto ub = match_u? lb+1: lb;

		if((n<(1<<5)&&std::is_integral_v<typename IU::type_data>)||n<(1<<2))
		{
			update_impl(u->lch, begin, begin+lb, assign);
			update_impl(u->rch, begin+ub, end, assign);
		}
		else
		{
			_spawn update_impl(u->lch, begin, begin+lb, assign);
			update_impl(u->rch, begin+ub, end, assign);
			_sync;
		}

		assign(ue, *(begin+lb));
		if constexpr(has_update<IU,typename IU::type_data&,Iter,Iter>::value)
		{
			const pbbs::range ra(begin,begin+lb);
			if(match_u)
			{
				/*
				auto f = [&](){
					_spawn IU::update(u->userdata, begin, begin+lb);
					IU::update(u->userdata, begin+lb, begin+ub);
					IU::update(u->userdata, begin+ub, end);
					_sync;
				};
				_spawn f();
				*/
				const auto rb = merge_one(begin+ub,end,ue,IU::type_data::IU::compare);
				const auto res = pbbs::merge(ra, rb, IU::type_data::IU::compare);
				// _sync;
				par_for(uint32_t i=0; i<res.size(); ++i)
					*(begin+i) = res[i];
				// std::inplace_merge(begin, begin+lb, begin+ub, IU::type_data::IU::compare);
				// std::inplace_merge(begin, begin+ub, end, IU::type_data::IU::compare);
			}
			else
			{
				/*
				auto f = [&](){
					_spawn IU::update(u->userdata, begin, begin+lb);
					IU::update(u->userdata, begin+ub, end);
					_sync;
				};
				_spawn f();
				*/
				const pbbs::range rb(begin+ub,end);
				const auto res = pbbs::merge(ra, rb, IU::type_data::IU::compare);
				// _sync;
				par_for(uint32_t i=0; i<res.size(); ++i)
					*(begin+i) = res[i];
				// std::inplace_merge(begin, begin+lb, end, IU::type_data::IU::compare);
			}
			IU::update(u->userdata, begin, end);
		}
		else u->update();
	}

	// return whether element is found and updated
	bool update(const T &element, void (*const assign)(T&, const T&)) // TODO: may support rvalue ref T_&&
	{
		return update_impl(root, element, assign);
	}

	template<typename Iter>
	void update(Iter begin, Iter end, void (*const assign)(T&, const T&))
	{
		// parallel_sort(begin, end, IU::compare);
		return update_impl(root, begin, end, assign);
	}

/*
	template<typename Iter>
	void update_impl(treap_node *curr, Iter begin, Iter end, void (*const assign)(T&, const T&))
	{
		if(!curr) return;
		if(begin==end) return;

		// TODO: eliminate the branches
		if(IU::compare(curr->key, *begin))
			update_impl(curr->rch, begin, end, assign);
		else if(IU::compare(*(end-1), curr->key))
			update_impl(curr->lch, begin, end, assign);
		else
		{
			assign(curr->key, *begin);
			auto it = std::lower_bound(begin, end, curr->key, IU::compare);
			if(!IU::compare()&&!IU::compare())
			{
				it--;
			}
			_spawn update_impl(curr, begin, it, assign);
			update_impl(curr, it, end, assign);
			_sync;
		}
		curr->update();
	}

	template<typename Iter>
	void update(Iter begin, Iter end, void (*const assign)(T&, const T&)){
		update_impl(root, begin, end, assign);
	}
*/

	template<typename Iter>
	treap_node* make_treap(Iter begin, Iter end, bool do_update=true, bool cnt_time=false)
	{
		const auto time_begin = high_resolution_clock::now();

		static_assert(std::is_base_of_v<
			std::random_access_iterator_tag,
			typename std::iterator_traits<Iter>::iterator_category
		>);
		const auto n = std::distance(begin, end);

		if constexpr(!std::is_integral_v<typename IU::type_data>)
			cerr << "Start sorting\n";
//		pbbs::sample_sort(begin, n, IU::compare);
		if constexpr(!std::is_integral_v<typename IU::type_data>)
			cerr << "Finished\n";
		// std::sort(begin, end, IU::compare);
//		auto t = allocator.allocate(n);
//		auto priority = new int[n];

		// std::mt19937 g(6);
		// std::random_device rd;
		// std::mt19937 g(rd());
		// std::shuffle(priority, priority+n, g);
		// auto seq_pri = pbbs::sequence<int>(priority,n);
		// pbbs::random_shuffle(seq_pri);
/*
		par_for(uint32_t i=0; i<n; ++i)
		{
			int priority = hash64(i+(uint64_t)this);
			if(priority<0) priority *= -1;
			t[i].priority = priority;
			t[i].key = *(begin+i);
		}
*/
		if constexpr(!std::is_integral_v<typename IU::type_data>)
			cerr << "Start Building Treap [" << n << "]\n";

		if constexpr(!std::is_integral_v<typename IU::type_data>)
		{
			auto copyied_input = new typename std::iterator_traits<Iter>::value_type[n];
			par_for(uint32_t i=0; i<n; ++i)
				copyied_input[i] = *(begin+i);
			begin = copyied_input;
		}

		treap_node *treap_root = build_treap(nullptr, n, do_update, begin, cnt_time);
		if constexpr(!std::is_integral_v<typename IU::type_data>)
			delete[] begin;

		const auto time_end = high_resolution_clock::now();
		const auto duration = time_end - time_begin;
		const auto msec = duration_cast<milliseconds>(duration).count();
		if(cnt_time)
			cerr << "Time of make_treap: " << msec << " ms\n";
		return treap_root;
	}
/*
	template<typename Iter>
	void insert_impl(treap_node *u, Iter begin, Iter end, void (*const assign)(T&, const T&))
	{
		const auto range = std::equal_range(begin, end, u->key, IU::compare);
	}
*/
	template<typename Iter>
	void insert(Iter begin, Iter end, void (*const assign)(T&, const T&))
	{
		const auto time_begin = high_resolution_clock::now();
		pbbs::sample_sort(begin, std::distance(begin,end), IU::compare);
		auto t = make_treap(begin, end, /*false*/true, false);
		root = treap_union(root, t, assign/*, true*/);
		const auto time_end = high_resolution_clock::now();
		const auto duration = time_end - time_begin;
		const auto msec = duration_cast<milliseconds>(duration).count();
		(void)msec;
		// cerr << "Time of insertion: " << msec << " ms\n";
	}

	interval() : root(nullptr)
	{
	}

	template<typename Iter>
	interval(Iter begin, Iter end) :
		root(make_treap(begin, end, true, !std::is_integral_v<typename IU::type_data>))
	{
	}

	interval(const T &element);

	interval(interval &&other)
		: root(other.root)
	{
		other.root = nullptr;
	}

	void destroy(treap_node *u)
	{
	//	u->key.~T();
	//	u->userdata.~typename IU::type_data();
		u->~treap_node();
		allocator.deallocate(u,1);
	}

	void destroy_sequential(treap_node *u)
	{
		if(u==nullptr) return;
		destroy_sequential(u->lch);
		destroy_sequential(u->rch);
		destroy(u);
	}

	void destroy_parallel(treap_node *u)
	{
		if(u==nullptr) return;
		if(u->size<(1<<9))
		{
			destroy_sequential(u);
			return;
		}
		_spawn destroy_parallel(u->lch);
		destroy_parallel(u->rch);
		destroy(u);
		_sync;
	}

	~interval()
	{
		destroy_parallel(root);
	}

	void copy(treap_node *&dst, const treap_node *src, treap_node *buffer)
	{
		if(src==nullptr) return;

		const auto lsize = src->lch? src->lch->size: 0;
		dst = &buffer[lsize];
		*dst = *src;
		if(lsize>(1<<9))
		{
			_spawn copy(dst->lch, src->lch, buffer);
		}
		else copy(dst->lch, src->lch, buffer);
		copy(dst->rch, src->rch, dst+1);
	}

	interval(const interval &other)
	{
		auto buffer = allocator.allocate(other.root->size);
		copy(root, other.root, buffer);
		/*
		const auto buffer = std::make_unique<T[]>(other.size());
		other.output(&buffer[0]);
		root = interval{&buffer[0], &buffer[other.size()]}.root;
		*/
	}

	void merge(interval &&other)
	{
		// std::set<uintptr_t> pset;
		// check_uniqueness(pset, root, 0);
		// check_uniqueness(pset, other.root, 0);
		root = treap_union(root, other.root, [](T&, const T&){});
		other.root = nullptr;
	}

	void check_uniqueness(std::set<uintptr_t> &pset, const treap_node *u, int d) const
	{
		if(u==nullptr) return;
		check_uniqueness(pset, u->lch, d+1);
		check_uniqueness(pset, u->rch, d+1);
		const auto res = pset.insert((uintptr_t)u);
		assert(res.second);
		// assert(d<30);
	}

	interval& operator=(interval &&other)
	{
		root = other.root;
		allocator = std::move(other.allocator);
		other.root = nullptr;
		return *this;
	}

	template<typename Interval>
	void replace_one(const Interval &other)
	{// TODO
		if(root==nullptr)
			root = allocator.allocate(1);
		assert(other.size()==1);
		unsigned int priority = hash64((uint64_t)root);
		root->priority = priority>>1;
		root->lch = root->rch = nullptr;
		root->key = other.root->key;
		root->userdata = other.root->userdata;
		root->size = other.root->size;
	}

	// interval& operator=(const interval&) = default;

	size_t size() const{
		return root? root->size: 0;
	}

	treap_component treap_extract(treap_node *root)
	{
		auto lch = root->lch;
		auto rch = root->rch;
		root->lch = root->rch = nullptr;
		return {lch, root, rch};
	}

	static treap_node* get_max_priority(treap_node *a, uint32_t n)
	{
		// if(n==1) return a;

		static const uint32_t threshold_scan = 1u<<6;
		if(n<threshold_scan)
		{
			uint32_t max_index = 0;
			auto max_priority = a[0].priority;
			for(uint32_t i=1; i<n; ++i)
				if(a[i].priority>max_priority)
				{
					max_index = i;
					max_priority = a[i].priority;
				}
			return &a[max_index];
		}

		treap_node *lmax=nullptr, *rmax=nullptr;
		_spawn [&](){lmax = get_max_priority(a, n/2);}();
		rmax = get_max_priority(a+n/2, n-n/2);
		_sync;

		return lmax->priority>rmax->priority?lmax:rmax;
	}

	template<typename Iter>
	treap_node* build_treap(treap_node *a, uint32_t n, bool do_update, const Iter begin, bool cnt_time=false)
	{
		if(n==0) return nullptr;

		a = allocator.allocate(1);

		int prio = hash64((uint64_t)a^(uint64_t)this);//g();
		if(prio<0) prio *= -1;
		a->priority = prio;

		const uint32_t p = n/2;
		a->key = *(begin+p);

		if(n==1)
		{
			a->lch = nullptr;
			a->rch = nullptr;
			if(do_update)
			{
				//memset(&a->userdata, 0, sizeof(a->userdata));
				// a->userdata = IU::g(a->key);
				a->update();
				//new(&a->userdata) typename IU::type_data(IU::g(a->key));
			}
			else a->update_size();
			return a;
		}

		const auto time_begin = high_resolution_clock::now();

		// uint32_t p = get_max_priority(a, n)-a;
		auto &ap = *a;
		/*
		std::swap(a[p], a[n-1]);
		auto pivot = a[n-1].key;
		p = partition(a, n-1, [=](const treap_node &k){
			return IU::compare(k.key, pivot);
		});
		std::swap(a[p], a[n-1]);
		*/

		if(p>(1<<9))
		{
			_spawn [&](){ap.lch = build_treap(a, p, do_update, begin);}();
			ap.rch = build_treap(a+p+1, n-p-1, do_update, begin+p+1);
			_sync;
		}
		else
		{
			ap.lch = build_treap(a, p, do_update, begin);
			ap.rch = build_treap(a+p+1, n-p-1, do_update, begin+p+1);
		}

		// auto ldata = a[p].lch?a[p].lch->userdata: typename IU::type_data();
		// auto rdata = a[p].rch?a[p].rch->userdata: typename IU::type_data();
		if(do_update)
		{
			if constexpr(!IU::is_nested)
			{
				// ap.userdata = IU::g(ap.key);
				ap.update();
			}
			else
			{
				// a[p].userdata = IU::g(a[p].key);
				// a[p].update();
				const auto sorted_lm = merge_one(begin,begin+p,*(begin+p),IU::type_data::IU::compare);
				const pbbs::range sorted_r(begin+p+1,begin+n);
				const auto sorted_whole = pbbs::merge(sorted_lm, sorted_r, IU::type_data::IU::compare);
				// assert(sorted_whole.size()==n);
				par_for(uint32_t i=0; i<n; ++i)
					*(begin+i) = sorted_whole[i];
				// a[p].userdata.copy_out(&sorted_whole[0], &sorted_whole[n], 0, UINT64_MAX, [](const auto u){return u;});
				// pbbs::sample_sort(&sorted_whole[0], n, IU::type_data::IU::compare);
				// for(uint32_t i=0; i<n; ++i)
					// assert(sorted_whole[i]==*(begin+i));
				ap.userdata = typename IU::type_data(begin, begin+n);
				ap.update_size();
			}
		}
		else ap.update_size();

		const auto time_end = high_resolution_clock::now();
		const auto duration = time_end - time_begin;
		const auto msec = duration_cast<milliseconds>(duration).count();
		if(cnt_time)
			cerr << "Time of build_treap: " << msec << " ms\n";
		return &ap;
	}

	treap_node* treap_join_direct(treap_node *tl, treap_node *k, treap_node *tr)
	{
		// TODO: correct the priority
		// if(tl) assert(k->priority>tl->priority);
		// if(tr) assert(k->priority>tr->priority);
		k->lch = tl, k->rch = tr;
		k->update();
		return k;
	}

	void rotate_right(treap_node *&k)
	{
		auto p = k;
		k = p->lch;
		p->lch = k->rch;
		k->rch = p;

		k->userdata = std::move(p->userdata);
		k->size = p->size;
		p->update();
	}

	void rotate_left(treap_node *&k)
	{
		auto p = k;
		k = p->rch;
		p->rch = k->lch;
		k->lch = p;

		k->userdata = std::move(p->userdata);
		k->size = p->size;
		p->update();
	}

	void treap_rebalance(treap_node *&k)
	{
		assert(k!=nullptr);

		auto kp = k->priority;
		auto lp = k->lch?k->lch->priority:-1;
		auto rp = k->rch?k->rch->priority:-1;

		if(kp>lp && kp>rp) return;

		if(lp>rp)
		{
			rotate_right(k);
			treap_rebalance(k->rch);
		}
		else
		{
			rotate_left(k);
			treap_rebalance(k->lch);
		}
	}

	treap_node* treap_join(treap_node *tl, treap_node *tm, treap_node *tr)
	{
		tm->lch = tl, tm->rch = tr;
		tm->update();
		if constexpr(std::is_integral_v<typename IU::type_data>)
			treap_rebalance(tm);
		return tm;
	}

	treap_component treap_split(treap_node *root, treap_node *k)
	{
		if(!root) return {nullptr, nullptr, nullptr};

		if(IU::compare(k->key,root->key))
		{
			auto [tl,tm,tr] = treap_split(root->lch, k);
			return {tl, tm, treap_join_direct(tr,root,root->rch)};
		}

		if(IU::compare(root->key,k->key))
		{
			auto [tl,tm,tr] = treap_split(root->rch, k);
			return {treap_join_direct(root->lch,root,tl), tm, tr};
		}

		//k->key==root->key
		return treap_extract(root);
	}

	// TODO: change the return type
	treap_component treap_split_last(treap_node *root)
	{
		auto [tl,tm,tr] = treap_extract(root);

		if(tr==nullptr) return {tl, tm, nullptr};

		auto [rtl,rtm,rtr] = treap_split_last(tr);
		(void)rtr; // rtr==nullptr
		return {treap_join(tl,tm,rtl), rtm, nullptr};
	}

	treap_node* treap_join(treap_node *tl, treap_node *tr)
	{
		if(tl==nullptr) return tr;

		auto [ltl,ltm,ltr] = treap_split_last(tl);
		(void)ltr; // ltr==nullptr
		return treap_join(ltl, ltm, tr);
	}

	// merge `t2` to `t1`
	[[nodiscard]] treap_node* treap_union(treap_node *t1, treap_node *t2, void (*const assign)(T&, const T&), bool update_t2=false)
	{
		if(!t2) return t1;
		if(!t1)
		{
			// if(update_t2) t2->update_recur();
			// cerr << t2->size << '\n';
			return t2;
		}
		// if(t1->priority>t2->priority) std::swap(t1, t2);

		auto [l2,m2,r2] = treap_extract(t2);
		auto [l1,m1,r1] = treap_split(t1, m2);
		// assert(!m1);

		// const auto lsize = (l1?l1->size:0) + (l2?l2->size:0);
		treap_node *tl, *tr;
		_spawn [&](){tl = treap_union(l1, l2, assign, update_t2);}();
		tr = treap_union(r1, r2, assign, update_t2);
		_sync;

		// return treap_join_direct(tl, m2, tr);

		// TODO: consider how to delete `m2`
		if(m1) // the same key exists in the original tree
		{
			assign(m1->key, m2->key);
		}
		else m1 = m2;
		return treap_join(tl, m1, tr);
	}

	void debug_output_impl(treap_node *curr) const
	{
		if(curr==nullptr) return;
		debug_output_impl(curr->lch);
		cerr << '(';
		if constexpr(std::is_same_v<T, uint64_t>)
			cerr << uint32_t(curr->key>>32) << ',' << uint32_t(curr->key);
		else
			cerr << curr->key;
		cerr << ',';
		cerr << curr->userdata << ") ";
		// cerr << curr->size << ") ";
		debug_output_impl(curr->rch);
	}

	void debug_output() const
	{
		cerr << '[';
		debug_output_impl(root);
		cerr << "]\n";
	}

	uint32_t check_height(const treap_node *u, const uint32_t d) const
	{
		if(u==nullptr) return d;
		const auto ld = _spawn check_height(u->lch, d+1);
		const auto rd = check_height(u->rch, d+1);
		_sync;
		return std::max(ld, rd);
	}
	uint32_t check_height() const
	{
		return check_height(root, 0);
	}

	template<typename Iter, typename Invocable>
	void output_impl(treap_node *curr, Iter it, Invocable f) const
	{
		if(curr==nullptr) return;
		const auto lsize = curr->lch?curr->lch->size:0;
		if(lsize>(1<<9))
		{
			_spawn output_impl(curr->lch, it, f);
			*(it+lsize) = f(curr->key);
			output_impl(curr->rch, it+(lsize+1), f);
			_sync;
		}
		else
		{
			output_impl(curr->lch, it, f);
			*(it+lsize) = f(curr->key);
			output_impl(curr->rch, it+(lsize+1), f);
		}
	}

	template<typename Iter>
	void output(Iter begin) const
	{
		output_impl(root, begin, [](const T &key){return key;});
	}

	template<typename Iter, class Invocable>
	void output(Iter begin, Invocable f) const
	{
		output_impl(root, begin, f);
	}

	template<typename Iter>
	uint32_t copy_useful_impl(treap_node *u, Iter iter)
	{
		if(u==nullptr) return 0;
		if(u->key>0) return 0; // u has NOT done
	}

	template<typename Iter>
	uint32_t copy_useful(Iter iter)
	{
		return copy_useful_impl(root, iter);
	}

	template<typename Iter, typename Invocable>
	uint32_t copy_out_left(treap_node *u, Iter it_s, const T &end, Invocable f) const
	{
		if(u==nullptr) return 0;

		const auto &ue = u->key;
		if(!IU::compare(ue,end))
		{
			const auto lsize = copy_out_left(u->lch, it_s, end, f);
			// u->update();
			return lsize;
		}

		const uint32_t lsize = u->lch?u->lch->size:0;
		if(lsize>(1<<9))
			_spawn output_impl(u->lch, it_s, f);
		else
			output_impl(u->lch, it_s, f);
		*(it_s+lsize) = f(u->key);
		const auto rsize = copy_out_left(u->rch, it_s+lsize+1, end, f);
		_sync;

		// u = u->rch;
		return lsize+1+rsize;
	}
	template<typename Iter, typename Invocable>
	uint32_t copy_out_right(treap_node *u, Iter it_t, const T &begin, Invocable f) const
	{
		if(u==nullptr) return 0;

		const auto &ue = u->key;
		if(!IU::compare(begin,ue))
		{
			const auto rsize = copy_out_right(u->rch, it_t, begin, f);
			// u->update();
			return rsize;
		}

		const uint32_t rsize = u->rch?u->rch->size:0;
		if(rsize>(1<<9))
			_spawn output_impl(u->rch, it_t-rsize, f);
		else
			output_impl(u->rch, it_t-rsize, f);
		*(it_t-(rsize+1)) = f(u->key);
		const auto lsize = copy_out_right(u->lch, it_t-(rsize+1), begin, f);
		_sync;

		// u = u->lch;
		return lsize+1+rsize;
	}
	template<typename Iter, typename Invocable>
	void copy_out_impl(treap_node *u, Iter it_s, Iter it_t, const T &begin, const T &end, Invocable f) const
	{
		if(u==nullptr) return;

		const auto &ue = u->key;
		// const uint32_t lsize = u->lch?u->lch->size:0;
		if(!IU::compare(begin,ue))
		{
			copy_out_impl(u->rch, it_s, it_t, begin, end, f);
			// u->update();
			return;
		}
		if(!IU::compare(ue,end))
		{
			copy_out_impl(u->lch, it_s, it_t, begin, end, f);
			// u->update();
			return;
		}

		const auto rsize = _spawn copy_out_right(u->lch, it_t, begin, f);
		const auto lsize = copy_out_left(u->rch, it_s, end, f);
		_sync;
		assert(it_s+lsize==it_t-(rsize+1));
		*(it_s+lsize) = f(u->key);
		// u = treap_join(u->lch, u->rch);
	}

	template<typename Iter, typename Invocable>
	void copy_out(Iter it_s, Iter it_t, const T &begin, const T &end, Invocable f) const
	{
		copy_out_impl(root, it_s, it_t, begin, end, f);
	}

	template<typename Iter>
	void remove_impl(treap_node *&u, Iter begin, Iter end)
	{
		if(u==nullptr) return;
		if(begin==end) return;

		const auto range = std::equal_range(begin, end, u->key, IU::compare);
		assert(range.second-range.first<=1);
		if(range.first-begin>(1<<9))
			_spawn remove_impl(u->lch, begin, range.first);
		else
			remove_impl(u->lch, begin, range.first);
		remove_impl(u->rch, range.second, end);
		_sync;
		if(range.first==range.second) // current key is not deleted
		{
			if constexpr(has_remove<IU, typename IU::type_data&,Iter,Iter>::value)
			{
				IU::remove(u->userdata, begin, end);
				u->update_size();
			}
			else u->update();
		}
		else
		{
		//	auto old_u = u;
			u = treap_join(u->lch, u->rch);
		//	allocator.deallocate(old_u, 1);
		}
	}

	template<typename Iter>
	void remove(Iter begin, Iter end)
	{
		remove_impl(root, begin, end);
	}

	uint32_t count_left(const treap_node *u, const T &end) const
	{
		/*
		if(u==nullptr) return 0;

		const auto &ue = u->key;
		if(!IU::compare(ue,end))
			return count_left(u->lch, end);

		auto res = 1+count_left(u->rch, end);
		if(u->lch)
			res += u->lch->size;
		return res;
		*/
		uint32_t res = 0;
		while(u!=nullptr)
		{
			const auto &ue = u->key;
			if(!IU::compare(ue,end))
			{
				u = u->lch;
				continue;
			}
			if(u->lch) res += u->lch->size;
			res++;
			u = u->rch;
		}
		return res;
	}
	uint32_t count_right(const treap_node *u, const T &begin) const
	{
		if(u==nullptr) return 0;

		const auto &ue = u->key;
		if(!IU::compare(begin,ue))
			return count_right(u->rch, begin);

		auto res = count_right(u->lch, begin)+1;
		if(u->rch)
			res += u->rch->size;
		return res;
	}
	uint32_t count_impl(treap_node *u, const T &begin, const T &end) const
	{
		if(u==nullptr) return 0;

		const auto &ue = u->key;
		if(!IU::compare(begin,ue)) return count_impl(u->rch, begin, end);
		if(!IU::compare(ue,end)) return count_impl(u->lch, begin, end);
		const auto lsize = count_right(u->lch, begin);
		const auto rsize = count_left(u->rch, end);
		return lsize+1+rsize;
	}
	uint32_t count(const T &begin, const T &end) const
	{
		return count_impl(root, begin, end);
	}

	// k starts from 0
	template<typename Invocable>
	std::pair<T,unsigned> get_kth(uint32_t k, Invocable g) const// g: count size
	{
		auto curr = root;
		unsigned int cnt = 0;
		for(; curr; cnt++)
		{
			uint32_t lsize = curr->lch?g(curr->lch->userdata):0;
			if(k<lsize) curr = curr->lch;
			else if(g(IU::g(curr->key,true))>0)
			{
				if(lsize==k) return {curr->key,cnt};
				else k-=lsize+1, curr=curr->rch;
			}
			else k-=lsize, curr=curr->rch;
		}
		assert(false);
	}

	void check_correctness(const treap_node *u) const
	{
		if(u==nullptr) return;
		check_correctness(u->lch);
		check_correctness(u->rch);
		auto augval = IU::g(u->key,true);
		if(u->lch)
			augval = IU::f(u->lch->userdata, augval);
		if(u->rch)
			augval = IU::f(augval, u->rch->userdata);
		assert(u->userdata==augval);
	}

	void check_correctness() const
	{
		check_correctness(root);
	}
};

// typedef std::tuple<treap_node*,bool,treap_node*> treap_debris;

// treap_node* to_treap(uint32_t *a, int *priority, uint32_t n, treap_node **data=nullptr);
// treap_component treap_split(treap_node *root, treap_node *k);
// treap_node* treap_union(treap_node *t1, treap_node *t2);

template<typename T, template<typename> class U, template<typename> class Alloc>
inline interval<T,U,Alloc>::interval(const T &element)
{
	root = allocator.allocate(1); // TODO
	root->priority = (unsigned int)(hash64((uint64_t)root))>>1;
	root->key = element;
	root->lch = root->rch = nullptr;
	root->update();
}
/*
template<typename T, template<typename> class U, template<typename> class Alloc>
template<typename Iter>
interval<T,U,Alloc>::interval(Iter begin, Iter end)
	: n(std::distance(begin, end))
{
	static_assert(std::is_base_of_v<
		std::random_access_iterator_tag,
		typename std::iterator_traits<Iter>::iterator_category
	>);
	std::sort(begin, end, IU::compare);
	auto t = new treap_node[n];
	// auto priority = new int[n];

	for(uint32_t i=0; i<n; ++i)
	{
		t[i].priority = i;//priority[i];
		t[i].key = *(begin+i);
		// t[i].userdata = IU::g(t[i].key);
	}
	root = build_treap(t, n);

	// delete priority;
}
*/
/*
template<typename T, template<typename> class U, template<typename> class Alloc>
typename interval<T,U,Alloc>::IU::type_data interval<T,U,Alloc>::query_impl(const treap_node *u, const T &begin, const T &end)
{
	if(u==nullptr) return typename IU::type_data();
	if(begin==end) return typename IU::type_data();

	const auto &ue = u->key;
	if(IU::compare(ue,begin))
		return query_impl(u->rch, begin, end);
	if(!IU::compare(ue,end))
		return query_impl(u->lch, begin, end);
	return IU::f(
		IU::f(query_impl(u->lch, begin, ue), IU::g(u->key)),
		query_impl(u->rch, ue+1, end)
	);
}
*/

template<typename T, template<typename> class U, template<typename> class Alloc>
template<typename G, typename F>
inline auto interval<T,U,Alloc>::query_left(const treap_node *u, const T &end, const G &retrieve_info, const F &merge_info) const
{
	unsigned int cnt = 0;
	auto res = retrieve_info(typename IU::type_data());
	for(; u!=nullptr; cnt++)
	{
		const auto &ue = u->key;
		if(!IU::compare(ue,end))
		{
			u = u->lch;
			continue;
		}

		if(u->lch)
		{
			res = merge_info(res, retrieve_info(u->lch->userdata));
			cnt++;
		}
		res = merge_info(
			res,
			retrieve_info(IU::g(u->key,true))
		);

		u = u->rch;
	}

	// if(u==nullptr) return std::make_pair(retrieve_info(typename IU::type_data()), 1u);
	// return std::make_pair(res,cnt);
	return res;
}

template<typename T, template<typename> class U, template<typename> class Alloc>
template<typename G, typename F>
inline auto interval<T,U,Alloc>::query_right(const treap_node *u, const T &begin, const G &retrieve_info, const F &merge_info) const
{
	if(u==nullptr) return retrieve_info(typename IU::type_data());

	const auto &ue = u->key;
	if(!IU::compare(begin,ue))
		return query_right(u->rch, begin, retrieve_info, merge_info);

	auto res = merge_info(
		query_right(u->lch, begin, retrieve_info, merge_info),
		retrieve_info(IU::g(u->key,true))
	);
	if(u->rch)
		res = merge_info(res, retrieve_info(u->rch->userdata));
	return res;
}

#endif // __INTERVAL_HPP_