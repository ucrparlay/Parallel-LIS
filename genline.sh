make clean
make

folder="result1"
#folder="/data0/zwan018/Parallel-LIS/line/"
round="1"
size_st="100m"
size=100000000
offsets=("0" "1" "10" "100" "1000" "10000")
steps=("30000" "100000" "300000" "1000000" "3000000" "10000000" "30000000" "100000000" )

echo "" > $folder/LIS_${size_st}_line.log
for offset in ${offsets[@]}; do
    for step in ${steps[@]}; do
        ./weighted_lis -a $size -l ${step} -u $size -r $round -p line -o $offset > $folder/LIS_${size_st}_${step}_${offset}_line.log
        cat $folder/LIS_${size_st}_${step}_${offset}_line.log >> $folder/LIS_${size_st}_${offset}_line.log
        rm $folder/LIS_${size_st}_${step}_${offset}_line.log
    done
done