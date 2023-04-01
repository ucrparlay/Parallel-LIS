make clean
make

folder="result1"

./weighted_lis -a 100000000 -l 1 -r 1 -p line > $folder/LIS_100m_1_line.log
./weighted_lis -a 100000000 -l 3 -r 1 -p line > $folder/LIS_100m_3_line.log
./weighted_lis -a 100000000 -l 10 -r 1 -p line > $folder/LIS_100m_10_line.log
./weighted_lis -a 100000000 -l 30 -r 1 -p line > $folder/LIS_100m_30_line.log
./weighted_lis -a 100000000 -l 100 -r 1 -p line > $folder/LIS_100m_100_line.log
./weighted_lis -a 100000000 -l 300 -r 1 -p line > $folder/LIS_100m_300_line.log
./weighted_lis -a 100000000 -l 1000 -r 1 -p line > $folder/LIS_100m_1000_line.log
./weighted_lis -a 100000000 -l 3000 -r 1 -p line > $folder/LIS_100m_3000_line.log
./weighted_lis -a 100000000 -l 10000 -r 1 -p line > $folder/LIS_100m_10000_line.log
./weighted_lis -a 100000000 -l 30000 -r 1 -p line > $folder/LIS_100m_30000_line.log
./weighted_lis -a 100000000 -l 100000 -r 1 -p line > $folder/LIS_100m_100000_line.log
./weighted_lis -a 100000000 -l 300000 -r 1 -p line > $folder/LIS_100m_300000_line.log
./weighted_lis -a 100000000 -l 1000000 -r 1 -p line > $folder/LIS_100m_1000000_line.log
./weighted_lis -a 100000000 -l 3000000 -r 1 -p line > $folder/LIS_100m_3000000_line.log
./weighted_lis -a 100000000 -l 10000000 -r 1 -p line > $folder/LIS_100m_10000000_line.log
./weighted_lis -a 100000000 -l 30000000 -r 1 -p line > $folder/LIS_100m_30000000_line.log
./weighted_lis -a 100000000 -l 100000000 -r 1 -p line > $folder/LIS_100m_100000000_line.log


cp $folder/LIS_100m_1_line.log $folder/LIS_100m_line.log
cat $folder/LIS_100m_3_line.log >> $folder/LIS_100m_line.log
cat $folder/LIS_100m_10_line.log >> $folder/LIS_100m_line.log
cat $folder/LIS_100m_30_line.log >> $folder/LIS_100m_line.log
cat $folder/LIS_100m_100_line.log >> $folder/LIS_100m_line.log
cat $folder/LIS_100m_300_line.log >> $folder/LIS_100m_line.log
cat $folder/LIS_100m_1000_line.log >> $folder/LIS_100m_line.log
cat $folder/LIS_100m_3000_line.log >> $folder/LIS_100m_line.log
cat $folder/LIS_100m_10000_line.log >> $folder/LIS_100m_line.log
cat $folder/LIS_100m_30000_line.log >> $folder/LIS_100m_line.log
cat $folder/LIS_100m_100000_line.log >> $folder/LIS_100m_line.log
cat $folder/LIS_100m_300000_line.log >> $folder/LIS_100m_line.log
cat $folder/LIS_100m_1000000_line.log >> $folder/LIS_100m_line.log
cat $folder/LIS_100m_3000000_line.log >> $folder/LIS_100m_line.log
cat $folder/LIS_100m_10000000_line.log >> $folder/LIS_100m_line.log
cat $folder/LIS_100m_30000000_line.log >> $folder/LIS_100m_line.log
cat $folder/LIS_100m_100000000_line.log >> $folder/LIS_100m_line.log
rm $folder/LIS_100m_1_line.log
rm $folder/LIS_100m_3_line.log
rm $folder/LIS_100m_10_line.log
rm $folder/LIS_100m_30_line.log
rm $folder/LIS_100m_100_line.log
rm $folder/LIS_100m_300_line.log
rm $folder/LIS_100m_1000_line.log
rm $folder/LIS_100m_3000_line.log
rm $folder/LIS_100m_10000_line.log
rm $folder/LIS_100m_30000_line.log
rm $folder/LIS_100m_100000_line.log
rm $folder/LIS_100m_300000_line.log
rm $folder/LIS_100m_1000000_line.log
rm $folder/LIS_100m_3000000_line.log
rm $folder/LIS_100m_10000000_line.log
rm $folder/LIS_100m_30000000_line.log
rm $folder/LIS_100m_100000000_line.log