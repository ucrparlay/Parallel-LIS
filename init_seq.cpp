#include <cmath>
using data_type = unsigned long long;

int initializeRandomWeight(sequence<size_t> &arr, size_t size, size_t limit, size_t seed = 0.01){
  parlay::parallel_for(0, size, [&](size_t i){
    arr[i]=hash64(i+seed*size)%(limit+1);
  });
  return 0;
}

int initializeRandomArray(sequence<data_type> &arr, size_t size, size_t limit, size_t lisLength, size_t seed = 0.01){
  size_t offset = limit-lisLength;
  if(offset<0)offset=0;
  parlay::parallel_for(0, size, [&](size_t i){
    arr[i]=hash64(i+seed*size)%(lisLength)+offset;
  });
  return 0;
}

int initializeLineArray(sequence<data_type> &arr, size_t size, size_t limit, size_t lisLength, size_t seed = 0.01){
  if(lisLength*lisLength>size){
    assert(size>=lisLength);
    size_t offset = size/lisLength;
    size_t base;
    if(limit<lisLength+offset)base=0;
    else base = limit-lisLength-offset;
    float gain = (float)(lisLength)/size;
    //cout<<"offset = "<<offset<<" base = "<<base<<" gain = "<<gain<<" lisLength-offset = "<<lisLength-offset<<endl;
    if(size>lisLength){
      parlay::parallel_for(0, size, [&](size_t i){
        arr[i]=hash64(i+seed*size)%(offset)+base+gain*i;
      });
    }else{
      parlay::parallel_for(0, size, [&](size_t i){
        arr[i]=base+i;
      });
    }
    return 0;
  }
  return initializeRandomArray(arr,size,limit,lisLength);
}
int initializeSegmentArray(sequence<data_type> &arr, size_t size, size_t limit, size_t lisLength, size_t seed = 0.01){
  parlay::parallel_for(0, size, [&](size_t i){
    arr[i]=hash64(i+seed*size)%(limit+1);
  });
  return 0;
}