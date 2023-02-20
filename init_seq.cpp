int initializeRandomWeight(sequence<size_t> &arr, size_t size, size_t limit, size_t seed = 0.01){
  parlay::parallel_for(0, size, [&](size_t i){
    arr[i]=hash64(i+seed*size)%limit;
  });
  return 0;
}