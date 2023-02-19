CC = g++
CFLAGS = -O3 -g3 -mcx16 -march=native -std=c++17 -ggdb -DPARLAY_CILKPLUS -fcilkplus -DCILK -I/home/csmajs/zwan018/.parlaylib/usr/local/include -I/home/csmajs/zwan018/pbbslib -I/home/zwan018/pbbslib -Wall -Wextra

all:	weighted_lis

weighted_lis:	weighted_lis.cpp
	$(CC) $(CFLAGS) weighted_lis.cpp basic_tools.cpp -o weighted_lis
  
clean:
	rm -f weighted_lis
