make clean
make
folder="/data5/zshen055/interval"
#folder="../LIS/in"

./weighted_lis -i $folder/LIS_100m_1_line.in -a 100000000 > result1/Para_LIS_100m_1_line.log
<<COMMENT
./weighted_lis -i $folder/LIS_100m_3_line.in -a 100000000 > result1/Para_LIS_100m_3_line.log
./weighted_lis -i $folder/LIS_100m_10_line.in -a 100000000 > result1/Para_LIS_100m_10_line.log
./weighted_lis -i $folder/LIS_100m_30_line.in -a 100000000 > result1/Para_LIS_100m_30_line.log
./weighted_lis -i $folder/LIS_100m_100_line.in -a 100000000 > result1/Para_LIS_100m_100_line.log
./weighted_lis -i $folder/LIS_100m_300_line.in -a 100000000 > result1/Para_LIS_100m_300_line.log
./weighted_lis -i $folder/LIS_100m_1000_line.in -a 100000000 > result1/Para_LIS_100m_1000_line.log
./weighted_lis -i $folder/LIS_100m_3000_line.in -a 100000000 > result1/Para_LIS_100m_3000_line.log
./weighted_lis -i $folder/LIS_100m_10000_line.in -a 100000000 > result1/Para_LIS_100m_10000_line.log
./weighted_lis -i $folder/LIS_100m_30000_line.in -a 100000000 > result1/Para_LIS_100m_30000_line.log
./weighted_lis -i $folder/LIS_100m_100000_line.in -a 100000000 > result1/Para_LIS_100m_100000_line.log
./weighted_lis -i $folder/LIS_100m_300000_line.in -a 100000000 > result1/Para_LIS_100m_300000_line.log
./weighted_lis -i $folder/LIS_100m_1000000_line.in -a 100000000 > result1/Para_LIS_100m_1000000_line.log
./weighted_lis -i $folder/LIS_100m_3000000_line.in -a 100000000 > result1/Para_LIS_100m_3000000_line.log
./weighted_lis -i $folder/LIS_100m_10000000_line.in -a 100000000 > result1/Para_LIS_100m_10000000_line.log
./weighted_lis -i $folder/LIS_100m_30000000_line.in -a 100000000 > result1/Para_LIS_100m_30000000_line.log
COMMENT
./weighted_lis -i $folder/LIS_100m_100000000_line.in -a 100000000 > result1/Para_LIS_100m_100000000_line.log


./weighted_lis -i $folder/LIS_100m_1_line.in 0 -a 100000000 -s > result1/Seq_LIS_100m_1_line.log
<<COMMENT
./weighted_lis -i $folder/LIS_100m_3_line.in 0 -a 100000000 -s > result1/Seq_LIS_100m_3_line.log
./weighted_lis -i $folder/LIS_100m_10_line.in 0 -a 100000000 -s > result1/Seq_LIS_100m_10_line.log
./weighted_lis -i $folder/LIS_100m_30_line.in 0 -a 100000000 -s > result1/Seq_LIS_100m_30_line.log
./weighted_lis -i $folder/LIS_100m_100_line.in 0 -a 100000000 -s > result1/Seq_LIS_100m_100_line.log
./weighted_lis -i $folder/LIS_100m_300_line.in 0 -a 100000000 -s > result1/Seq_LIS_100m_300_line.log
./weighted_lis -i $folder/LIS_100m_1000_line.in 0 -a 100000000 -s > result1/Seq_LIS_100m_1000_line.log
./weighted_lis -i $folder/LIS_100m_3000_line.in 0 -a 100000000 -s > result1/Seq_LIS_100m_3000_line.log
./weighted_lis -i $folder/LIS_100m_10000_line.in 0 -a 100000000 -s > result1/Seq_LIS_100m_10000_line.log
./weighted_lis -i $folder/LIS_100m_30000_line.in 0 -a 100000000 -s > result1/Seq_LIS_100m_30000_line.log
./weighted_lis -i $folder/LIS_100m_100000_line.in 0 -a 100000000 -s > result1/Seq_LIS_100m_100000_line.log
./weighted_lis -i $folder/LIS_100m_300000_line.in 0 -a 100000000 -s > result1/Seq_LIS_100m_300000_line.log
./weighted_lis -i $folder/LIS_100m_1000000_line.in 0 -a 100000000 -s > result1/Seq_LIS_100m_1000000_line.log
./weighted_lis -i $folder/LIS_100m_3000000_line.in 0 -a 100000000 -s > result1/Seq_LIS_100m_3000000_line.log
./weighted_lis -i $folder/LIS_100m_10000000_line.in 0 -a 100000000 -s > result1/Seq_LIS_100m_10000000_line.log
./weighted_lis -i $folder/LIS_100m_30000000_line.in 0 -a 100000000 -s > result1/Seq_LIS_100m_30000000_line.log
COMMENT
./weighted_lis -i $folder/LIS_100m_100000000_line.in 0 -a 100000000 -s > result1/Seq_LIS_100m_100000000_line.log
