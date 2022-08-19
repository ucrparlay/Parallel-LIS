#ifndef BASIC_TOOLS_H
#define BASIC_TOOLS_H

#include <iostream>
#include <parlay/sequence.h>
using namespace std;
using namespace parlay;
using data_type = unsigned long long;

void printArray(sequence<data_type> &arr, size_t size, size_t start);
void printArray2(sequence<size_t> &arr, size_t size, size_t start);
void testArray(sequence<data_type> &arr, size_t internalEnd, int rounds);

#endif
