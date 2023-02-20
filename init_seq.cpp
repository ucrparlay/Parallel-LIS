using data_type = unsigned long long;

int initializeRandomWeight(sequence<size_t> &arr, size_t size, size_t limit, size_t seed = 0.01){
  parlay::parallel_for(0, size, [&](size_t i){
    arr[i]=hash64(i+seed*size)%(limit+1);
  });
  return 0;
}

int initializeRandomArray(sequence<data_type> &arr, size_t size, float limit, size_t seed = 0.01){
	const uint32_t cnt_worker = 400;
	auto *subgen = new std::mt19937[cnt_worker];
	{
		std::mt19937 g(seed);
		for(uint32_t ii=0; ii<cnt_worker; ++ii)
			subgen[ii] = std::mt19937(g());
	}
  float gain = limit/size*0.75;
  auto &g = subgen[0];
  parlay::parallel_for(0, size, [&](data_type i){
    auto offset = g()%uint32_t(limit-size*gain);
    auto val = (size-i)*gain + offset;
    arr[i]=(data_type)val;
  });
  return 0;
}