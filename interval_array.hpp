#ifndef __INTERVAL_ARRAY_HPP_
#define __INTERVAL_ARRAY_HPP_

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
#include <atomic>
#include <type_traits>
#include <sample_sort.h>
#include "primitive.hpp"
#include "util.hpp"
#include "interval_common.hpp"
using std::cerr;
using std::cout;
using std::endl;
using namespace std::chrono;

template<typename T, template<typename> class U, template<typename> class Alloc=std::allocator>
struct interval_array
{
	typedef U<T> IU;

	struct treap_node{
		T key;
		typename IU::type_data userdata;
		// uint32_t size;
	};
	Alloc<treap_node> allocator;
	completeBST layout;
	treap_node *data_arena;

	inline treap_node* get_lch(const treap_node *u) const
	{
		const uint32_t index = (u-data_arena)*2;
		return index<=size()? &data_arena[index]: nullptr;
	}

	inline treap_node* get_rch(const treap_node *u) const
	{
		const uint32_t index = (u-data_arena)*2+1;
		return index<=size()? &data_arena[index]: nullptr;
	}

	inline treap_node* get_parent(const treap_node *u) const
	{
		const uint32_t index = (u-data_arena)/2;
		return index>0? &data_arena[index]: nullptr;
	}

	auto get_userdata_update_(const treap_node *u) const
	{
		const auto lch = get_lch(u);
		const auto rch = get_rch(u);
		auto userdata = IU::g(u->key);
		if(lch) userdata = IU::f(lch->userdata, userdata);
		if(rch) userdata = IU::f(userdata, rch->userdata);
		return userdata;
	}

	void node_update(treap_node *u)
	{
		u->userdata = get_userdata_update_(u);
		// node_update_size(u);
	}

	/*
	inline void node_update_size(treap_node *u)
	{
		const auto lch = get_lch(u);
		const auto rch = get_rch(u);
		u->size = 1;
		if(lch) u->size += lch->size;
		if(rch) u->size += rch->size;
	}
	*/

	template<typename G, typename F>
	auto query_left(const treap_node *u, const T &end, const G &retrieve_info, const F &merge_info) const
	{
		unsigned int cnt = 0;
		auto res = retrieve_info(typename IU::type_data());
		for(; u!=nullptr; cnt++)
		{
			const auto &ue = u->key;
			if(!IU::compare(ue,end))
			{
				u = get_lch(u);
				continue;
			}

			if(get_lch(u))
			{
				res = merge_info(res, retrieve_info(get_lch(u)->userdata));
				cnt++;
			}
			res = merge_info(
				res,
				retrieve_info(IU::g(u->key,true))
			);

			u = get_rch(u);
		}

		// if(u==nullptr) return std::make_pair(retrieve_info(typename IU::type_data()), 1u);
		// return std::make_pair(res,cnt);
		return res;
	}

	template<typename G, typename F>
	auto query_left(const T &end, const G &retrieve_info, const F &merge_info) const
	{
		return size()? query_left(&data_arena[1], end, retrieve_info, merge_info): 0;
	}

	auto query_left(const T &end) const
	{
		return size()? query_left(&data_arena[1], end, [](const typename IU::type_data &x){return x;}, IU::f): 0;
	}

	template<typename G, typename F>
	auto query_right(const treap_node *u, const T &begin, const G &retrieve_info, const F &merge_info) const
	{
		if(u==nullptr) return retrieve_info(typename IU::type_data());

		const auto &ue = u->key;
		auto *rch = get_rch(u);
		if(!IU::compare(begin,ue))
			return query_right(rch, begin, retrieve_info, merge_info);

		auto res = merge_info(
			query_right(get_lch(u), begin, retrieve_info, merge_info),
			retrieve_info(IU::g(u->key,true))
		);
		if(rch)
			res = merge_info(res, retrieve_info(rch->userdata));
		return res;
	}

	template<typename G, typename F>
	auto query_right(const T &begin, const G &retrieve_info, const F &merge_info) const
	{
		return query_right(&data_arena[1], begin, retrieve_info, merge_info);
	}

	auto query_right(const T &begin) const
	{
		return query_right(&data_arena[1], begin, [](const typename IU::type_data &x){return x;}, IU::f);
	}

	template<typename G>
	auto query_pivot(const treap_node *u, const G &retrieve_info) const
	{
		while(true)
		{
			if(retrieve_info(IU::g(u->key,true))<0)
				return -long(u->key);
			if(get_rch(u) && retrieve_info(get_rch(u)->userdata)<0)
				u = get_rch(u);
			else u = get_lch(u);
		}
	}

	template<typename G, typename F>
	auto query_special(const treap_node *u, const T &end, const G &retrieve_info, const F &merge_info) const
	{
		if(u==nullptr) return retrieve_info(typename IU::type_data());

		const auto &ue = u->key;
		if(!IU::compare(ue,end)) // ue >= end
			return query_special(get_lch(u), end, retrieve_info, merge_info);

		const auto info_mid = retrieve_info(IU::g(u->key,true));
		auto res = merge_info(
			info_mid<0? -long(ue): info_mid,
			query_special(get_rch(u), end, retrieve_info, merge_info)
		);
		if(res<0) return res;
		if(get_lch(u))
		{
			const auto info_left = retrieve_info(get_lch(u)->userdata);
			res = merge_info(info_left<0? query_pivot(get_lch(u), retrieve_info): info_left, res);
			if(res<0) return res;
		}
		return res;
	}

	template<typename G, typename F>
	auto query_special(const T &end, const G &retrieve_info, const F &merge_info) const
	{
		return query_special(&data_arena[1], end, retrieve_info, merge_info);
	}

	template<typename Iter>
	void update_impl(treap_node *u, Iter begin, Iter end, void (*const assign)(T&, const T&))
	{
		if(u==nullptr || begin==end) return;

		const auto n = std::distance(begin, end);
		auto &ue = u->key;
		// const auto lb = partition(begin, n, [&](const auto &x){return IU::compare(x,ue);});
		const auto lb = std::lower_bound(begin, end, ue, IU::compare)-begin;
		
		const bool match_u = lb!=n&&!IU::compare(ue,*(begin+lb));
		const auto ub = match_u? lb+1: lb;

		if(n<(1<<10))
		{
			update_impl(get_lch(u), begin, begin+lb, assign);
			update_impl(get_rch(u), begin+ub, end, assign);
		}
		else
		{
			_spawn update_impl(get_lch(u), begin, begin+lb, assign);
			update_impl(get_rch(u), begin+ub, end, assign);
			_sync;
		}

		assign(ue, *(begin+lb));
		if constexpr(has_update<IU,typename IU::type_data&,Iter,Iter>::value)
		{
			const pbbs::range ra(begin,begin+lb);
			if(match_u)
			{
				const auto rb = merge_one(begin+ub,end,ue,IU::type_data::IU::compare);
				const auto res = pbbs::merge(ra, rb, IU::type_data::IU::compare);
				par_for(uint32_t i=0; i<res.size(); ++i)
					*(begin+i) = res[i];
			}
			else
			{
				const pbbs::range rb(begin+ub,end);
				const auto res = pbbs::merge(ra, rb, IU::type_data::IU::compare);
				par_for(uint32_t i=0; i<res.size(); ++i)
					*(begin+i) = res[i];
			}
			IU::update(u->userdata, begin, end);
		}
		else
		{
			node_update(u);
		}
	}

	template<typename Iter>
	void update(Iter begin, Iter end, void (*const assign)(T&, const T&))
	{
		return update_impl(&data_arena[1], begin, end, assign);
	}

	// rank counts from 0
	template<typename Comp>
	bool update_kth(const T &element, void (*const assign)(T&, const T&), uint32_t rank, Comp comp)
	{
		const auto n = size();
		if(rank>=n) return false;

		const auto index = layout.rank2index(rank+1);
		treap_node *u = &data_arena[index];
		assert(!IU::compare(u->key,element)&&!IU::compare(element,u->key));
		//assign(u->key, element);

		for(; u; u=get_parent(u))
		{
			#if !defined(_M_X86)&&!defined(__i386__)
			#warning "Unable to confirm the x86-compatible architecture. \
				Retrieved data may be incorrect."
			#endif
			const auto userdata_new = get_userdata_update_(u);
			if(!pbbs::write_max(&u->userdata, userdata_new, comp))
				return false;
		}
		return true;
	}

	interval_array() :
		layout(0), data_arena(nullptr)
	{
	}

	template<typename Iter>
	interval_array(Iter begin, Iter end) :
		layout(std::distance(begin, end)),
		data_arena(allocator.allocate(layout.size()+1))
	//	data_arena((treap_node*)new char[(layout.size()+1)*sizeof(treap_node)])
	{
		const auto time_begin = high_resolution_clock::now();

		static_assert(std::is_base_of_v<
			std::random_access_iterator_tag,
			typename std::iterator_traits<Iter>::iterator_category
		>);
		const auto n = size();
		if(n==0) return;

		if constexpr(IU::is_nested)
			cerr << "Start Building Treap [" << n << "]\n";

		if constexpr(IU::is_nested)
		{
			auto copyied_input = new typename std::iterator_traits<Iter>::value_type[n];
			par_for(uint32_t i=0; i<n; ++i)
				copyied_input[i] = *(begin+i);
			begin = copyied_input;
		}

		par_for(uint32_t i=1; i<=n; ++i) // i is index
		{
			const auto rank = layout.index2rank(i); // rank starts from 1
			data_arena[i].key = *(begin+rank-1);
		}

		update_array(&data_arena[1]);
		if constexpr(IU::is_nested)
			delete[] begin;

		const auto time_end = high_resolution_clock::now();
		const auto duration = time_end - time_begin;
		const auto msec = duration_cast<milliseconds>(duration).count();

	//	cerr << "Time of make_treap: " << msec << " ms\n";
	}

	interval_array(const T &element);

	interval_array(interval_array &&other) :
		layout(other.layout), data_arena(other.data_arena)
	{
		other.data_arena = nullptr;
	}

	interval_array(const interval_array &other) :
		layout(other.layout),
		data_arena(allocator.allocate(layout.size()+1))

	{
		const auto n = size();
		par_for(uint32_t i=1; i<=n; ++i)
			data_arena[i] = other.data_arena[i];
	}

	interval_array& operator=(interval_array &&other)
	{
		layout = other.layout;
		data_arena = other.data_arena;
		other.data_arena = nullptr;
		return *this;
	}
	// interval_array& operator=(const interval_array&) = default;

	size_t size() const{
		return data_arena? layout.size(): 0;
	}

	void update_array(treap_node *u)
	{
		if(!u) return;

		auto *lch=get_lch(u), *rch=get_rch(u);

		if(layout.get_size(u-data_arena)>(1<<10))
		{
			_spawn update_array(lch);
			update_array(rch);
			_sync;
		}
		else
		{
			update_array(lch);
			update_array(rch);
		}
		// TODO: special handle for nested tree
		static_assert(!IU::is_nested);
		node_update(u);
		/*
		auto &userdata = u->userdata;
		
		userdata = IU::g(u->key);

		u->size = 1;
		if(lch) userdata = IU::f(lch->userdata, userdata), u->size += lch->size;
		if(rch) userdata = IU::f(userdata, rch->userdata), u->size += rch->size;
		*/
	}

	template<typename G, typename F>
	auto query_impl(const treap_node *u, const T &begin, const T &end, const G &retrieve_info, const F &merge_info) const
	{
		if(u==nullptr) return retrieve_info(typename IU::type_data());

		const auto &ue = u->key;
		if(!IU::compare(begin,ue)) return query_impl(get_rch(u), begin, end, retrieve_info, merge_info);
		if(!IU::compare(ue,end)) return query_impl(get_lch(u), begin, end, retrieve_info, merge_info);
		const auto res_left = _spawn query_right(get_lch(u), begin, retrieve_info, merge_info);
		const auto res_right = query_left(get_rch(u), end, retrieve_info, merge_info);
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
		return query_impl(&data_arena[1], begin, end, retrieve_info, merge_info);
	}
	auto query(const T &begin, const T &end) const
	{
		return query_impl(&data_arena[1], begin, end, [](const typename IU::type_data &x){return x;}, IU::f);
	}

	template<class G>
	auto query(const G &retrieve_info) const
	{
		return retrieve_info(data_arena[1].userdata);
	}

	auto query() const
	{
		return query([](const typename IU::type_data &x){return x;});
	}

	~interval_array()
	{
		const auto n = size();
		if(n>0)
			allocator.deallocate(data_arena, n);
	}
};

template<typename T, template<typename> class U, template<typename> class Alloc>
inline interval_array<T,U,Alloc>::interval_array(const T &element)
	: layout(1), data_arena(allocator.allocate(2))
{
	// data_arena->priority = (unsigned int)(hash64((uint64_t)data_arena))>>1;
	auto *root = &data_arena[1];
	root->key = element;
	// root->lchx = root->rchx = nullptr;
	// root->size = 1;
	node_update(root);
}

#endif // __INTERVAL_ARRAY_HPP_