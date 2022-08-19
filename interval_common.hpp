#ifndef __INTERVAL_COMMON_HPP_
#define __INTERVAL_COMMON_HPP_

namespace detail
{
template<class T, class ...Args>
auto has_insert_impl(int) ->
	decltype(void(std::declval<T>().insert(std::declval<Args>()...)), std::true_type{});
template<class, class...>
auto has_insert_impl(...) -> std::false_type;

template<class T, class ...Args>
auto has_remove_impl(int) ->
	decltype(void(std::declval<T>().remove(std::declval<Args>()...)), std::true_type{});
template<class, class...>
auto has_remove_impl(...) -> std::false_type;

template<class T, class ...Args>
auto has_update_impl(int) ->
	decltype(void(std::declval<T>().update(std::declval<Args>()...)), std::true_type{});
template<class, class...>
auto has_update_impl(...) -> std::false_type;

template<class T, class ...Args>
auto has_replace_one_impl(int) ->
	decltype(void(std::declval<T>().replace_one(std::declval<Args>()...)), std::true_type{});
template<class, class...>
auto has_replace_one_impl(...) -> std::false_type;
} // namespace detail

template<class T, class ...Args>
struct has_insert :
	decltype(detail::has_insert_impl<T, Args...>(0)){};

template<class T, class ...Args>
struct has_remove :
	decltype(detail::has_remove_impl<T, Args...>(0)){};

template<class T, class ...Args>
struct has_update :
	decltype(detail::has_update_impl<T, Args...>(0)){};

template<class T, class ...Args>
struct has_replace_one :
	decltype(detail::has_replace_one_impl<T, Args...>(0)){};


#endif // __INTERVAL_COMMON_HPP_