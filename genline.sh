make clean
make

folder="result1"
#folder="/data0/zwan018/Parallel-LIS/line/"
round="10"
size_st="100m"
size=100000000
steps=("1" "3" "10" "30" "100" "300" "1000" "3000" "10000" "30000" "100000" "300000" "1000000" "3000000" "10000000" "30000000" "100000000" )

echo "" > $folder/LIS_${size_st}_line.log
for step in ${steps[@]}; do
    ./weighted_lis -a $size -l ${step} -u $size -r $round -p line -g $folder/LIS_${size_st}_${offset}_${step}_segment.log > $folder/LIS_${size_st}_${step}_line.log
    cat $folder/LIS_${size_st}_${step}_line.log >> $folder/LIS_${size_st}_line.log
    rm $folder/LIS_${size_st}_${step}_line.log
done