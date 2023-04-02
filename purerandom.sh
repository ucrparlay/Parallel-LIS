make clean
make

folder="result1"
#folder="/data0/zwan018/Parallel-LIS/pureRandom/"
round="10"
limit=1000000000
sizes=("1000000" "10000000" "100000000" "1000000000")
cores=("2" "4" "8" "16" "32" )
#"64" "128" "192"
# 
echo "" > $folder/LIS_PureRandom.log
for size in ${sizes[@]}; do
    echo ${size} > $folder/c.txt
    ./weighted_lis -a  ${size} -u $limit -r $round -p pure  -g $folder/LIS_${size}_PureRandom.in > $folder/LIS_tmp_PureRandom.log
    cat $folder/c.txt >> $folder/LIS_PureRandom.log
    cat $folder/LIS_tmp_PureRandom.log >> $folder/LIS_PureRandom.log
    rm $folder/LIS_tmp_PureRandom.log
done
rm $folder/c.txt

for core in ${cores[@]}; do
    echo ${core} > $folder/c.txt
    numactl -i ${core} ./weighted_lis -i $folder/LIS_1000000000_PureRandom.in > $folder/LIS_tmp_PureRandom.log
    cat $folder/c.txt >> $folder/LIS_PureRandom.log
    cat $folder/LIS_tmp_PureRandom.log >> $folder/LIS_PureRandom.log
    rm $folder/LIS_tmp_PureRandom.log
done
rm $folder/c.txt