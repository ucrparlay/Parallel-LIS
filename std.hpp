#define ENABLE_PARALLELISM 0
#include <cstdio>
#include <algorithm>
#include <cstring>
#include <limits>
#include <set>
#include <chrono>
#include <cassert>
#include <random>
#include <memory>
// #include "parallelize.h"
#define N 200000005
using namespace std::chrono;

struct zz{
	int a, f;
	bool operator<(const zz &other) const{
		return a<other.a;
	}
	bool operator==(const zz &other) const{
		return a==other.a;
	}
	bool operator>(const zz &other) const{
		return a>other.a;
	}
};

template <typename T>
class AVL{
public:
	struct Tnode{
		T data;
		int h,size;
		Tnode *lch,*rch;
		int max_f;
		inline void recount()
		{
			if(this==lch) return;
			int lh=lch->h;
			int rh=rch->h;
			h=(lh>rh?lh:rh)+1;
			size=lch->size+rch->size+1;
			max_f = std::max(std::max(lch->max_f,rch->max_f), data.f);
		}
	};
private:
	Tnode *lib,*mm,*null,*root;
	int size_;
	void L(Tnode *&t);
	void R(Tnode *&t);
	void Balance_L(Tnode *&t);
	void Balance_R(Tnode *&t);
	void Balance(Tnode *&t);
	bool insert_(const T &x,Tnode *&t);
	int get_max_(const int end, Tnode *t);
public:
	AVL();
	void insert(const T &x) { if(insert_(x,root)) size_++; }
	int get_max(const int end) { return get_max_(end, root); }
	int size(){ return size_; }
};

template <typename T>
AVL<T>::AVL():size_(0)
{
	root=null=mm=lib=new Tnode[N];
	*null=(Tnode){T{0,0},0,0,null,null,0};
}

template <typename T>
void AVL<T>::Balance_L(Tnode *&t)
{
	int lh=t->lch->lch->h;
	int rh=t->lch->rch->h;
	if(lh>=rh) R(t);else L(t->lch),R(t);
}

template <typename T>
void AVL<T>::Balance_R(Tnode *&t)
{
	int lh=t->rch->lch->h;
	int rh=t->rch->rch->h;
	if(rh>=lh) L(t);else R(t->rch),L(t);
}

template <typename T>
void AVL<T>::Balance(Tnode *&t)
{
	int dh=t->lch->h-t->rch->h;
	if(dh>1) Balance_L(t);
	if(dh<-1)Balance_R(t);
}

template <typename T>
void AVL<T>::L(Tnode *&t)
{

	Tnode *x=t;
	t=x->rch;
	x->rch=t->lch;
	t->lch=x;
	x->recount();
	t->recount();
}

template <typename T>
void AVL<T>::R(Tnode *&t)
{

	Tnode *x=t;
	t=x->lch;
	x->lch=t->rch;
	t->rch=x;
	x->recount();
	t->recount();
}

template <typename T>
bool AVL<T>::insert_(const T &x, Tnode *&t)
{
	if(t==null)
	{
		t = ++mm;
		*t = (Tnode){x,1,1,null,null,x.f};
		return true;
	}
	bool did_insertion;
	if(x==t->data){
    if(t->data.f<x.f){
      t->data=x;
    }
    did_insertion=false;
  }
	else if(x<t->data)
		did_insertion=insert_(x,t->lch), Balance(t);
	else did_insertion=insert_(x,t->rch), Balance(t);
	t->recount();
	return did_insertion;
}

template <typename T>
int AVL<T>::get_max_(const int end, Tnode *t)
{
	if(t==null) return 0;

	const int a = t->data.a;
	if(!(a<end))
		return get_max_(end, t->lch);

	const int f_mid = t->data.f;
	int res = std::max(f_mid, get_max_(end, t->rch));
	if(t->lch!=null)
		res = std::max(t->lch->max_f, res);
	return res;
}

uint32_t load_data(int argc, char **argv, uint32_t *output)
{
	assert(argc>=2);
	printf("read from %s\n", argv[0]);
	FILE *fin = fopen(argv[0], "r");

	uint32_t n, n_or=atoi(argv[1]);
	printf("override indicator: %u\n", n_or);

//	scanf("%u", &n);
	printf("read original n: %u\n", n);
	if(n_or!=0)
	{
//		assert(n>=n_or);
		n = n_or;
	}
	printf("n eventually is: %.2e\n", double(n));

	for(uint32_t i=1; i<=n; ++i)
	{
		fscanf(fin, "%u", &output[i]);
	}
	fclose(fin);
	return n;
}

uint32_t generator_line(int argc, char **argv, uint32_t *output)
{
	assert(argc>=5);
	const uint32_t seed = atoi(argv[0]);
	uint32_t n = atoi(argv[1]);
	const uint32_t pattern = atoi(argv[2]);
	const uint32_t range_y = atoi(argv[3]);
	const double gain = double(atoi(argv[4]))*1e-6;
	const uint32_t base_offset = argc>=6? atoi(argv[5]): 0;
	const uint32_t range_offset = argc>=7? atoi(argv[6]): 0;
	const uint32_t cnt_worker = 400;

	printf("seed = %u\n", seed);
	printf("n = %.2e\n", double(n));
	printf("pattern = %c\n", "lo"[pattern]);
	printf("range_y = %.2e\n", double(range_y));
	printf("gain = %.2f\n", gain);
	printf("base_offset = %u\n", base_offset);
	printf("range_offset = %u\n", range_offset);
//	printf("stddev = %.2e\n", double(stddev));
	auto *subgen = new std::mt19937[cnt_worker];
	{
		std::mt19937 g(seed);
		for(uint32_t ii=0; ii<cnt_worker; ++ii)
			subgen[ii] = std::mt19937(g());
	}

	std::unique_ptr<std::normal_distribution<>> distr_normal;
	switch(pattern)
	{
		case 0:
			assert(range_y>n*gain && double(range_y)/n>gain);
			break;
		case 1:
			assert(base_offset+range_offset<range_y);
			assert((n-1)*gain+base_offset>0);
			assert((n-1)*gain+base_offset+range_offset<range_y);
			break;
		//	distr_normal.reset(new std::normal_distribution<>(0,stddev));
	}

	for(uint32_t i=0; i<n; ++i)
	{
		//auto &g = subgen[__cilkrts_get_worker_number()];
		auto &g = subgen[0];
		uint32_t value;
		switch(pattern)
		{
			case 0: // uniform distribution
			{
				auto offset = g()%uint32_t(range_y-n*gain);
				value = (n-i)*gain + offset;
				break;
			}
			case 1: // normal distribution
			{
				auto offset = g()%range_offset;
				value = i*gain + offset + base_offset;
				/*
				auto drift = fabs((*distr_normal)(g));
				if(drift>range_timestamp/2) drift = range_timestamp/2; // clamp
				const uint32_t len = std::round(drift)+range_seg/2;

				begin = g()%(range_timestamp-len)+1;
				end = begin+len;
				*/
				break;
			}
			default:
			{}
		}
		output[i+1] = value;
	}
	delete []subgen;

	return n;
}

uint32_t generator_seg(int argc, char **argv, uint32_t *output)
{
	assert(argc>=4);
	const uint32_t seed = atoi(argv[0]);
	const uint32_t n = atoi(argv[1]);
	const uint32_t k = n/atoi(argv[2]);
	const uint32_t range_y = atoi(argv[3]);
	const uint32_t cnt_worker = 400;

	printf("seed = %u\n", seed);
	printf("n = %.2e\n", double(n));
	printf("k = %u\n", n/k);
	printf("range_y = %.2e\n", double(range_y));
	auto *subgen = new std::mt19937[cnt_worker];
	{
		std::mt19937 g(seed);
		for(uint32_t ii=0; ii<cnt_worker; ++ii)
			subgen[ii] = std::mt19937(g());
	}

	const uint32_t size_seg = n/k;
	const uint32_t range_seg = range_y/k;
	std::set<uint32_t> s;
	while(s.size()<size_seg)
	{
		auto &g = subgen[0];
		const uint32_t num = uint32_t(g()%range_seg+1);
		s.insert(num);
	}

	uint32_t *buffer = new uint32_t[size_seg];
	uint32_t offset = 0;
	for(const auto e : s)
	{
		buffer[offset++] = e;
	}
	for(uint32_t ii=0; ii<k; ++ii)
	{
		for(uint32_t j=0; j<offset; ++j)
		{
			output[offset*ii+j+1] =  buffer[j]*k+(k-ii);
		}
	}
	delete []buffer;

	return s.size()*k;
}

uint32_t *u_backup = new uint32_t[N];

uint32_t prepare_data(int argc, char **argv)
{
	assert(argc>=1);
	decltype(load_data) *data_source[] = {load_data, generator_line, generator_seg};
	const uint32_t source_id = atoi(argv[0]);
	assert(source_id<sizeof(data_source)/sizeof(*data_source));

	auto *output = u_backup;
	uint32_t n = data_source[source_id](argc-1, argv+1, output);

	return n;
}

void output_data(uint32_t n, const char *filename, const uint32_t *data)
{
	FILE *fout = fopen(filename, "w");
	for(uint32_t i=1; i<=n; ++i)
	{
		fprintf(fout, "%u\n", data[i]);
	}
	fclose(fout);
}
/*
int main(int argc, char **argv)
{
	for(int i=0; i<argc; ++i)
		printf("%s ", argv[i]);
	putchar('\n');
	const uint32_t n = prepare_data(argc-2, argv+2);
	uint32_t *u_shrink = new uint32_t[N];
	std::mt19937 g(1206);
	for(uint32_t i=0; i<n/100; ++i)
	{
		u_shrink[i+1] = u_backup[i*100+g()%100+1];
	}
	output_data(n/100, argv[1], u_shrink);
	return 0;

	auto *inputdata = u_backup;

	const auto time_begin = high_resolution_clock::now();
	auto ans = std::numeric_limits<int>::min();
	AVL<zz> a;
	a.insert(zz{ans,0});
	for(uint32_t i=1; i<=n; ++i)
	{
		int x = inputdata[i];
		int f = a.get_max(x)+1;
		if(f>ans) ans = f;
		a.insert(zz{x,f});
	}


	const auto time_end = high_resolution_clock::now();
	const auto duration = time_end - time_begin;
	const auto msec = duration_cast<milliseconds>(duration).count();

	printf("ans: %d\n", ans);
	printf("size: %d\n", a.size());
	printf("%.2f s\n", double(msec)/1000);
	return 0;
}
*/