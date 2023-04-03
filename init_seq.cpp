#include <cmath>
using data_type = unsigned long long;

int initializeRandomWeight(sequence<size_t> &arr, size_t size, size_t limit, size_t seed = 0.01){
  parlay::parallel_for(0, size, [&](size_t i){
    arr[i]=hash64(i+seed*size)%(limit+1);
  });
  return 0;
}

int initializePureRandomArray(sequence<data_type> &arr, size_t size, size_t limit, size_t lisLength, size_t seed = 0.01){
  parlay::parallel_for(0, size, [&](size_t i){
    arr[i]=hash64(i+seed*size)%(limit);
  });
  return 0;
}

int initializeRandomArray(sequence<data_type> &arr, size_t size, size_t limit, size_t lisLength, size_t seed = 0.01){
  size_t base = limit-lisLength;
  if(base<0)base=0;
  parlay::parallel_for(0, size, [&](size_t i){
    arr[i]=hash64(i+seed*size)%(lisLength)+base;
  });
  return 0;
}

int initializeLineArray(sequence<data_type> &arr, size_t size, size_t limit, size_t lisLength, size_t seed = 0.01, size_t offset = 10){
  assert(size>=lisLength);
  size_t base=0;
  if(limit<lisLength+offset)base=0;
  else base = limit-lisLength-offset;
  float gain = (float)(lisLength)/size;
  cout<<"offset = "<<offset<<" base = "<<base<<" gain = "<<gain<<" lisLength-offset = "<<lisLength-offset<<endl;
  if(offset!=0){
    parlay::parallel_for(0, size, [&](size_t i){
      arr[i]=hash64(i+seed*size)%(offset)+base+gain*i;
    });
  }else{
    parlay::parallel_for(0, size, [&](size_t i){
      arr[i]=base+gain*i;
    });
  }
  return 0;
}
int initializeSegmentArray(sequence<data_type> &arr, size_t size, size_t limit, size_t lisLength, size_t seed = 0.01, size_t offset = 10){
  assert(size>=lisLength);
  if(size==lisLength){
    parlay::parallel_for(0, size, [&](size_t i){
      arr[i]=i;
    });
    return 0;
  }
  size_t seg = size/lisLength; //x-length of each segment
  float gap = (float)limit/lisLength; //y-length of each segment
  float lad = (float)limit/size; //y-length of each ladder
  parlay::parallel_for(0, lisLength, [&](size_t i){
    parlay::parallel_for(0, seg, [&](size_t j){
      arr[i*seg+j]=(i+1)*gap-lad*j+hash64(i*seg+j+seed*size)%((size_t)offset);
    });
  });
  return 0;
}