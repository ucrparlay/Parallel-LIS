CC = g++
CFLAGS = -O3 -g3 -mcx16 -march=native -std=c++17 -ggdb -DPARLAY_CILKPLUS -fcilkplus -DCILK -Wall -Wextra -I.parlaylib/include -Ipbbslib

all:	weighted_lis 

weighted_lis:	weighted_lis.cpp init_seq.cpp
	$(CC) $(CFLAGS) weighted_lis.cpp basic_tools.cpp -o weighted_lis
  
clean:
	rm -f weighted_lis
