#include "basic_tools.h"

void printArray(sequence<data_type> &arr, size_t size, size_t start = 0){
  for(size_t i = start; i < start+size; ++i){
    if(arr[i]==ULLONG_MAX){
      cout<<"@@\t";
    }
    else{
      cout << arr[i] <<"\t";
    }
  }
  cout<<endl;
}

void printArray2(sequence<size_t> &arr, size_t size, size_t start = 0){
  for(size_t i = start; i < start+size; ++i){
    if(arr[i]==ULLONG_MAX){
      cout<<"@@\t";
    }
    else{
      cout << arr[i] <<"\t";
    }
  }
  cout<<endl;
}


void testArray(sequence<data_type> &arr, size_t internalEnd, int rounds = -1){
    cout<<"ROUND: "<<rounds<<endl;
    printArray(arr, internalEnd+1, 0);
    printArray(arr, internalEnd+2, internalEnd+1);
    cout<<endl;
}