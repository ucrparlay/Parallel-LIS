#include <iostream>
#include <fstream>
#include <cstdio>
#include <stdlib.h>
#include <cmath>
#include <vector>
#include <algorithm>
#include <ctime>
#include <queue>
#include <string>
#include "parlay/sequence.h"
#include "parlay/parallel.h"
#include "parlay/primitives.h"
#include "parlay/utilities.h"
#include <cstdlib>
#include <cstdint>
#include <cassert>
#include <cstring>
#include <iterator>
#include <functional>
#include <memory>
#include <chrono>
#include <atomic>
#include <typeinfo>
#include "parlay/alloc.h"
#include "interval.hpp"
#include "std.hpp"
#include "interval_array.hpp"
#include "init_seq.cpp"

using namespace std;
using data_type1 = uint64_t;
using data_type = unsigned long long;
#include "basic_tools.h"

size_t ARRAY_SIZE = 1000000000;
float WEIGHT_LIMIT = 10;
size_t ARRAY_LIMIT = 1e6;
size_t LIS_LENGTH = 10;
constexpr size_t SIZE_LIMIT = 4294967295;
constexpr data_type LOWER_LIMIT = 0;
constexpr data_type UPPER_LIMIT = 1e9;
constexpr size_t GRANULARITY = 1024;
constexpr size_t DEPTH_GRANULARITY = 16;//1e9-19
constexpr size_t BATCH_SIZE = 1e4;
enum Pattern { rando = 0, line, segment, pure };

timer t_seq0;
timer t_seq;
timer t_seq1;
timer t_para;
timer t_buildTree;
timer t_findPivot;
timer t_handleWeight;
timer t_queryLeft;
timer t_update;
timer t_prepare;
timer t_unweighted;
timer t_tmp;

typedef uint64_t nid_t;

inline constexpr uint32_t get_y(nid_t u){
	return uint32_t(u>>32);
}

inline constexpr uint32_t get_x(nid_t u){
	return uint32_t(u);
}


inline constexpr nid_t gen_nid(uint32_t x, uint32_t y){
	return (nid_t(y)<<32)|x;
}

inline nid_t gen_nid2(uint32_t x, uint32_t y){
	return (nid_t(y)<<32)|x;
}

long *f = new long[ARRAY_SIZE];
auto * a = new nid_t[ARRAY_SIZE];


template<typename E>
struct U_inner{
	const static constexpr bool is_nested_alloc = false;
	const static constexpr bool is_nested = false;
	typedef long type_data; // the maximum DP value
	static bool compare(const E lhs, const E rhs){
		return lhs<rhs;
	}
	static type_data f(const type_data &l, const type_data &r){
		if(l<=0 && r<=0) return l+r;
		if(l>=0 && r>=0) return std::max(l,r);
		return l<0?l:r;
	}

	static type_data g(const E &element, ...){
		return ::f[get_x(element)]; // TODO
	}
};


template<typename E>
struct U_outer{
	typedef interval_array<nid_t, U_inner> type_data;
	const static constexpr bool is_nested_alloc = true;
	const static constexpr bool is_nested = true;
	static bool compare(const E lhs, const E rhs){
		return get_x(lhs)<get_x(rhs);
	}

	static type_data f(const type_data &ld, const type_data &rd){
		return type_data();
	}

	static type_data g(const E &element){
		return {element};
	}

	static auto g(const E &element, bool){
		return interval_array<nid_t, U_inner>{element};
	}
	template<typename Iter>
	static void update(type_data &data, Iter begin, Iter end)
	{
		data.update(begin, end, [](nid_t&, const nid_t&){});
	}
};

template <typename T>
struct allocator_wrapper/* : pbbs::type_allocator<T>*/
{
	static T* allocate(std::size_t n)
	{
		if(n==1)
			return pbbs::type_allocator<T>::alloc();
		assert(n==1);
		return new T[n];
	}

	void deallocate(T* p, std::size_t n)
	{
		assert(n==1);
		pbbs::type_allocator<T>::free(p);
	}
};


size_t buildTree(sequence<data_type> &arr, sequence<data_type> &arr2, size_t size){
  size_t depth = ceil(log2(size));
  size_t newSize = pow(2, depth+1)-1;
  size_t s = pow(2, depth) - 1;
  
  arr.resize(newSize);
  par_for(size_t i=s; i<s+size; ++i){
    arr[i]=arr2[i-s];
  }
  par_for(size_t i=s+size; i<newSize; ++i){
    arr[i]=ULLONG_MAX;
  }
  
  par_for(size_t i=0; i<s; ++i){
    arr[i]=0;
  }
  return depth;
}

data_type buildMinTreeSeq(sequence<data_type> &nodeArray, size_t root, size_t secondLastEnd){
  size_t leftChild = 2*root+1;
  size_t rightChild = 2*root+2;
  if(root > secondLastEnd){
    return min(nodeArray[leftChild], nodeArray[rightChild]);
  }
  
  nodeArray[leftChild] = buildMinTreeSeq(nodeArray, leftChild, secondLastEnd);
  nodeArray[rightChild] = buildMinTreeSeq(nodeArray, rightChild, secondLastEnd);
  return min(nodeArray[leftChild], nodeArray[rightChild]);
}

data_type buildMinTree(sequence<data_type> &nodeArray, size_t root, size_t secondLastEnd, size_t res){
  size_t leftChild = 2*root+1;
  size_t rightChild = 2*root+2;
  if(root > secondLastEnd){
    return min(nodeArray[leftChild], nodeArray[rightChild]);
  }
  if(root<res){
    parlay::par_do(
      [&]() { nodeArray[leftChild] = buildMinTree(nodeArray, leftChild, secondLastEnd, res); },
      [&]() { nodeArray[rightChild] = buildMinTree(nodeArray, rightChild, secondLastEnd, res); });
  }else{
    nodeArray[leftChild] = buildMinTreeSeq(nodeArray, leftChild, secondLastEnd);
    nodeArray[rightChild] = buildMinTreeSeq(nodeArray, rightChild, secondLastEnd);
  }
  return min(nodeArray[leftChild], nodeArray[rightChild]);
}

data_type findPivotSeq0(sequence<data_type> &nodeArray, sequence<size_t> &rankArray, size_t root, size_t internalEnd, data_type pre, size_t round){
  if(nodeArray[root] == ULLONG_MAX){
    return ULLONG_MAX;
  }
  if(root > internalEnd){
    rankArray[root-internalEnd-1] = round;
    return ULLONG_MAX;
  }
  size_t leftChild = 2*root+1;
  size_t rightChild = 2*root+2;
  
  if(nodeArray[root]==nodeArray[rightChild]){
    if(nodeArray[leftChild]<=pre&&nodeArray[leftChild]!=ULLONG_MAX){
      nodeArray[rightChild] = findPivotSeq0(nodeArray, rankArray, rightChild, internalEnd, nodeArray[leftChild], round);
      nodeArray[leftChild] = findPivotSeq0(nodeArray, rankArray, leftChild, internalEnd, pre, round);
    }
    else{
      nodeArray[rightChild] = findPivotSeq0(nodeArray, rankArray, rightChild, internalEnd, pre, round);
    }
  }else{
    nodeArray[leftChild] = findPivotSeq0(nodeArray, rankArray, leftChild, internalEnd, pre, round);
  }
  return min(nodeArray[rightChild], nodeArray[leftChild]);
}

data_type findPivot(sequence<data_type> &nodeArray, sequence<size_t> &rankArray, size_t root, size_t internalEnd, data_type pre, size_t res, size_t round){
  if(nodeArray[root] == ULLONG_MAX){
    return ULLONG_MAX;
  }
  if(root>=res){
    return findPivotSeq0(nodeArray, rankArray, root, internalEnd, pre, round);
  }
  size_t leftChild = 2*root+1;
  size_t rightChild = 2*root+2;
  if(nodeArray[root]==nodeArray[rightChild]){
    if(nodeArray[leftChild]<=pre&&nodeArray[leftChild]!=ULLONG_MAX){
      data_type lc = nodeArray[leftChild];
      parlay::par_do(
        [&]() { nodeArray[leftChild] = findPivot(nodeArray, rankArray, leftChild, internalEnd, pre, res, round); },
        [&]() { nodeArray[rightChild] = findPivot(nodeArray, rankArray, rightChild, internalEnd, lc, res, round); });
    }
    else{
      nodeArray[rightChild] = findPivot(nodeArray, rankArray, rightChild, internalEnd, pre, res, round);
    }
  }else{
    nodeArray[leftChild] = findPivot(nodeArray, rankArray, leftChild, internalEnd, pre, res, round);
  }
  return min(nodeArray[leftChild], nodeArray[rightChild]);
}



size_t calWeight(sequence<data_type> &initialArray, sequence<size_t> &weightArray, sequence<size_t> &rankArray, size_t size, size_t rounds){
  t_prepare.start();
  auto * workset = new size_t[rounds+1];
  
  parlay::sequence<pair<size_t,size_t>> indexReal(size);
  par_for(size_t i=0; i<size; ++i){
    f[i] = 0;
    a[i] = gen_nid2(i, (uint32_t)initialArray[i]);
    indexReal[i].first = rankArray[i];
    indexReal[i].second = i;
  }
  interval<nid_t, U_outer, allocator_wrapper> q(a, a+size);
  auto less = [&] (pair<size_t,size_t> A, pair<size_t,size_t> B) {
    if(A.first == B.first) return A.second < B.second;
    return A.first < B.first;};
  parlay::sort_inplace(indexReal, less);

  workset[0] = 0;
  workset[rounds] = size;
  uint32_t k = indexReal[0].second;
  a[0] = gen_nid(k, (uint32_t)initialArray[k]);
  
  par_for(size_t i=1; i<size; ++i){
    k = indexReal[i].second;
    a[i] = gen_nid(k, (uint32_t)initialArray[k]);
    if(indexReal[i].first != indexReal[i-1].first){
      workset[indexReal[i].first]=i;
    }
  }
  t_prepare.stop();
  
  for(size_t round = 0; round < rounds; ++round){
    size_t begin = workset[round], end = workset[round+1];
    t_queryLeft.start();
    par_for(size_t i=begin; i<end; ++i){
        const auto u = a[i];
        const auto ind = get_x(u);
        const auto r = q.query_left(u,[&,u](const auto &inner){
            const auto r = inner.query_left(gen_nid(0,get_y(u)));
            return r;
          },
          [](const auto &l, const auto &r){
            if(l<0 && r<0)
              return r;
            if(l>=0 && r>=0) return std::max(l,r);
            return l<0?l:r;
          }
        );
        f[ind] = r + weightArray[ind];
    }
    t_queryLeft.stop();
    t_update.start();
    q.update(a+begin, a+end, [](nid_t&, const nid_t&){});
    t_update.stop();
  }
  t_queryLeft.start();
  const auto ans = q.query_left(size,[&,size](const auto &inner){
      const auto r = inner.query_left(gen_nid(SIZE_LIMIT,SIZE_LIMIT));
      return r;
    },
    [](const auto &l, const auto &r){
      if(l<0 && r<0)
        return r;
      if(l>=0 && r>=0) return std::max(l,r);
      return l<0?l:r;
    }
  );
  t_queryLeft.stop();
  return ans;
}

size_t runWeightedParallel(sequence<data_type> &initialArray, sequence<size_t> weightArray, size_t size, size_t gra=DEPTH_GRANULARITY){
  sequence<data_type> nodeArray;
  sequence<size_t> rankArray(size);
  
  t_buildTree.start();
  size_t depth = buildTree(nodeArray, initialArray, size);
  size_t internalEnd = pow(2, depth) - 2;
  size_t secondLastEnd = pow(2, depth-1) - 2;
  size_t res = 0;
  if(depth>gra) res = pow(2,depth-gra) -1;
  nodeArray[0] = buildMinTree(nodeArray, 0, secondLastEnd, res);
  t_buildTree.stop();
  
  size_t rounds = 0;
  while(nodeArray[0]!=ULLONG_MAX){
    t_findPivot.start();
    nodeArray[0] = findPivot(nodeArray, rankArray, 0, internalEnd, ULLONG_MAX, res, rounds);
    t_findPivot.stop();
    rounds++;
  }
  t_handleWeight.start();
  size_t ans = calWeight(initialArray, weightArray, rankArray, size, rounds);
  t_handleWeight.stop();
  return ans;
}

int runWeightedSequential(sequence<data_type> &initialArray, sequence<size_t> &weightArray, size_t size){
	auto ans = std::numeric_limits<int>::min();
	AVL<zz> avl;
	avl.insert(zz{ans,0});
	for(uint32_t i=0; i<size; ++i)
	{
		int x = initialArray[i];
		int k = avl.get_max(x)+weightArray[i];
		if(k>ans) ans = k;
		avl.insert(zz{x,k});
  }
  return ans;
}

data_type findPivotUnweightedSeq(data_type pre, sequence<data_type> &nodeArray, size_t root, size_t internalEnd){
  if(nodeArray[root] == ULLONG_MAX || root > internalEnd){
    return ULLONG_MAX;
  }
  size_t leftChild = 2*root+1;
  size_t rightChild = 2*root+2;
  
  if(nodeArray[root]==nodeArray[rightChild]){
    if(nodeArray[leftChild]<=pre&&nodeArray[leftChild]!=ULLONG_MAX){
      nodeArray[rightChild] = findPivotUnweightedSeq(nodeArray[leftChild], nodeArray, rightChild, internalEnd);
      nodeArray[leftChild] = findPivotUnweightedSeq(pre, nodeArray, leftChild, internalEnd);
    }
    else{
      nodeArray[rightChild] = findPivotUnweightedSeq(pre, nodeArray, rightChild, internalEnd);
    }
  }else{
    nodeArray[leftChild] = findPivotUnweightedSeq(pre, nodeArray, leftChild, internalEnd);
  }
  return min(nodeArray[rightChild], nodeArray[leftChild]);
}

data_type findPivotUnweighted(data_type pre, sequence<data_type> &nodeArray, size_t root, size_t internalEnd, size_t res){
  if(nodeArray[root] == ULLONG_MAX){
    return ULLONG_MAX;
  }
  if(root>=res){
    return findPivotUnweightedSeq(pre, nodeArray, root, internalEnd);
  }
  size_t leftChild = 2*root+1;
  size_t rightChild = 2*root+2;
  if(nodeArray[root]==nodeArray[rightChild]){
    if(nodeArray[leftChild]<=pre&&nodeArray[leftChild]!=ULLONG_MAX){
      data_type lc = nodeArray[leftChild];
      parlay::par_do(
        [&]() { nodeArray[leftChild] = findPivotUnweighted(pre, nodeArray, leftChild, internalEnd, res); },
        [&]() { nodeArray[rightChild] = findPivotUnweighted(lc, nodeArray, rightChild, internalEnd, res); });
    }
    else{
      nodeArray[rightChild] = findPivotUnweighted(pre, nodeArray, rightChild, internalEnd, res);
    }
  }else{
    nodeArray[leftChild] = findPivotUnweighted(pre, nodeArray, leftChild, internalEnd, res);
  }
  return min(nodeArray[leftChild], nodeArray[rightChild]);
}

size_t runParallelUnweighted(sequence<data_type> nodeArray, size_t depth,  size_t size, size_t gra){
  size_t rounds = 0;
  t_buildTree.start();
  size_t internalEnd = pow(2, depth) - 2;
  size_t secondLastEnd = pow(2, depth-1) - 2;
  size_t res = 0;
  if(depth>gra) res = pow(2,depth-gra) -1;
  nodeArray[0] = buildMinTree(nodeArray, 0, secondLastEnd, res);
  t_buildTree.stop();
  
  while(nodeArray[0]!=ULLONG_MAX){
    t_findPivot.start();
    nodeArray[0] = findPivotUnweighted(ULLONG_MAX, nodeArray, 0, internalEnd, res);
    t_findPivot.stop();
    rounds++;
  }
  return rounds;
}

size_t runUnweightedSequential0(sequence<data_type> nodeArray, size_t depth, size_t size){
  sequence<data_type> lis(size+1, 0);
  size_t leafStart = pow(2, depth)-1;
  lis[1] = nodeArray[leafStart];
  size_t s = 1;
  
  for(size_t i=1; i<size; ++i){
    data_type x = nodeArray[leafStart+i];
    if(x>lis[s])lis[++s]=x;
    else{
      size_t l=1, h=s, m;
      while(l<=h){
        m=(l+h)/2;
        if(x>lis[m])l=m+1;
        else h=m-1;
      }
      lis[l] = x;
    }
  }
  return s;
}
int runUnweightedSequential(sequence<data_type> nodeArray, size_t depth){
  int rounds = 0;
  size_t secondLastEnd = pow(2, depth-1) - 2;
  size_t internalEnd = pow(2, depth) - 2;
  nodeArray[0] = buildMinTreeSeq(nodeArray, 0, secondLastEnd);
  while(nodeArray[0]!=ULLONG_MAX){
    nodeArray[0] = findPivotUnweightedSeq(ULLONG_MAX, nodeArray, 0, internalEnd);
    rounds++;
  }
  return rounds;
}

int run_final_weighted(int ROUND, sequence<data_type> &initialArray, sequence<size_t> &weightArray, size_t size, size_t gra = DEPTH_GRANULARITY, bool seq = false){
  if(seq){
    t_seq.reset();
    t_seq.start();
    int ans1 = runWeightedSequential(initialArray, weightArray, size);
    t_seq.stop();
    //cout<<"ans\ttime"<<endl;
    cout<<ans1<<"\t"<<t_seq.get_total()<<endl;
  }else{
    size_t ans = runWeightedParallel(initialArray, weightArray, size, gra);
    t_buildTree.reset();
    t_findPivot.reset();
    t_handleWeight.reset();
    t_prepare.reset();
    t_queryLeft.reset();
    t_update.reset();
    t_para.reset();
    t_para.start();
    for(int i=0;i<ROUND;++i){
      runWeightedParallel(initialArray, weightArray, size, gra);
    }
    t_para.stop();
    //cout<<"ans\ttotal_time\tfind_pivot\thandle_weight\tprepare\t\tquery_left\tupdate"<<endl;
    cout<< ans <<"\t"<<t_para.get_total()/ROUND;
    cout<<"\t"<<t_findPivot.get_total()/ROUND<<"\t"<<t_handleWeight.get_total()/ROUND;
    cout<<"\t"<<t_prepare.get_total()/ROUND<<"\t"<<t_queryLeft.get_total()/ROUND<<"\t"<<t_update.get_total()/ROUND;
    cout<<endl;
  }
  return 0;
}

int run_final_unweighted(int ROUND, sequence<data_type> &initialArray, size_t size, size_t gra = DEPTH_GRANULARITY, bool seq = false){
  sequence<data_type> nodeArray;
  size_t depth = buildTree(nodeArray, initialArray, size);
  if(seq){
    //array
    t_seq0.reset();
    t_seq0.start();
    int ans0 = runUnweightedSequential0(nodeArray, depth, size);
    t_seq0.stop();
    //work
    t_seq.reset();
    t_seq.start();
    int ans2 = runUnweightedSequential(nodeArray, depth);
    t_seq.stop();
    //cout<<"greedy\t"<< ans0 <<"\t"<<t_seq0.get_total()<<endl;
    //cout<<"seq\t"<< ans2 <<"\t"<<t_seq.get_total()<<endl;
    cout<< ans0 <<"\t"<<t_seq0.get_total()<<"\t"<<t_seq.get_total()<<endl;
  }else{
    size_t ans = runParallelUnweighted(nodeArray, depth, size, gra);
    t_buildTree.reset();
    t_findPivot.reset();
    t_para.reset();
    t_para.start();
    for(int i=0;i<ROUND;++i){
      runParallelUnweighted(nodeArray, depth, size, gra);
    }
    t_para.stop();
    //cout<<"ans\ttotal_time\tfind_pivot\tbuild_tree"<<endl;
    cout<< ans <<"\t"<<t_para.get_total()/ROUND<<"\t"<<t_findPivot.get_total()/ROUND<<"\t"<<t_buildTree.get_total()/ROUND;
    cout<<endl;
  }
  return 0;
}
int main(int argc, char* argv[]){
  int ROUND = 3;
  if (argc == 1) {
    fprintf(
        stderr,
        "Usage: %s [-i input_file] [-a array_size] [-u array_limit] [-l lis_length] [-e weight_limit] [-r test_round] [-g generated_input] [-p data_pattern] [-w] [-s]\n"
        "Options:\n"
        "\t-i,\tif you have generated input file, please type in the file path\n"
        "\t-a,\tsize of input array\n"
        "\t-u,\tthe upper limit of array values, if autogenerated\n"
        "\t-l,\tthe length of LIS, if autogenerated\n"
        "\t-e,\tneed to generate weight, please set weight limit\n"
        "\t-r,\tnumebr of test round\n"
        "\t-g,\tthe filename of the autogenerated input data\n"
        "\t-w,\tweighted LIS\n"
        "\t-s,\trun in sequential\n",
        "\t-p,\tdata pattern: line, segment or random\n",
        argv[0]);
    return 0;
  }
  
  ifstream infile;
  string ofs;
  
  bool gval1=false,gval2=false;
  bool gval3=false;
  bool weighted=false, seq=false;
  Pattern pat = rando;
  size_t offset=10;
  char c;
  while ((c = getopt(argc, argv, "i:a:u:e:r:a:l:g:p:o:ws")) != -1) {
    switch (c) {
      case 'i':
        infile.open(optarg, ifstream::in);
        if(!infile.is_open()){
          cout<<"No Input File"<<endl;
          return 0;
        }
        break;
      case 'a':
        ARRAY_SIZE = atoi(optarg);
        break;
      case 'u':
        ARRAY_LIMIT = atoi(optarg);
        break;
      case 'l':
        gval1 = true;
        LIS_LENGTH = atoi(optarg);
        break;
      case 'e':
        gval2 = true;
        WEIGHT_LIMIT = atof(optarg);
        break;
      case 'r':
        ROUND = atoi(optarg);
        break;
      case 'g':
        gval3 = true;
        ofs = optarg;
        break;
      case 'w':
        weighted = true;
        break;
      case 's':
        seq = true;
        break;
      case 'o':
        offset = atoi(optarg);
        break;
      case 'p':
        gval1 = true;
        if (!strcmp(optarg, "line")) {
          pat = line;
        } else if (!strcmp(optarg, "segment")) {
          pat = segment;
        } else if (!strcmp(optarg, "random")) {
          pat = rando;
        } else if (!strcmp(optarg, "pure")) {
          pat = pure;
        } else {
          fprintf(stderr, "Error: Unknown pattern %s\n", optarg);
          exit(EXIT_FAILURE);
        }
        break;
      default:
        fprintf(stderr, "Error: Unknown option %c\n", optopt);
          exit(EXIT_FAILURE);
    }
  }
  ARRAY_LIMIT = max(ARRAY_LIMIT, LIS_LENGTH);
  pbbs::type_allocator<interval<nid_t, U_outer, allocator_wrapper>::treap_node>::reserve(ARRAY_SIZE);
  sequence<data_type> initialArray(ARRAY_SIZE);
  if(gval1){
    if(pat == rando)initializeRandomArray(initialArray, ARRAY_SIZE, ARRAY_LIMIT, LIS_LENGTH, 0.001);
    if(pat == line)initializeLineArray(initialArray, ARRAY_SIZE, ARRAY_LIMIT, LIS_LENGTH, 0.001, offset);
    if(pat == segment)initializeSegmentArray(initialArray, ARRAY_SIZE, ARRAY_LIMIT, LIS_LENGTH, 0.001, offset);
    if(pat == pure)initializePureRandomArray(initialArray, ARRAY_SIZE, ARRAY_LIMIT, LIS_LENGTH, 0.001);
  }else{
    for(size_t i=0; i<ARRAY_SIZE; ++i)
      infile>>initialArray[i];
  }
  if(weighted){
    sequence<size_t> weightArray(ARRAY_SIZE);
    if(gval2){
      initializeRandomWeight(weightArray, ARRAY_SIZE, WEIGHT_LIMIT, 0.001);
    }else{
      for(size_t i=0; i<ARRAY_SIZE; ++i)
        infile>>weightArray[i];
    }
    run_final_weighted(ROUND, initialArray, weightArray, ARRAY_SIZE, DEPTH_GRANULARITY, seq);
    if(gval3){
      ofstream myfile(ofs, std::ofstream::out);
      if(gval1)for(size_t i=0;i<ARRAY_SIZE;++i)myfile<<initialArray[i]<<"\n";
      if(gval2)for(size_t i=0;i<ARRAY_SIZE;++i)myfile<<weightArray[i]<<"\n";
    }
    return 0;
  }
  if(gval3){
    ofstream myfile(ofs, std::ofstream::out);
    if(gval1)for(size_t i=0;i<ARRAY_SIZE;++i)myfile<<initialArray[i]<<"\n";
  }
  run_final_unweighted(ROUND, initialArray, ARRAY_SIZE, DEPTH_GRANULARITY, seq);
  return 0;
}
