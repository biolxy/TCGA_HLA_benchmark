for item in $(ls *class2 result.t)
do
    echo $item
    for i in DRB1 DPB1 DQA1 DQB1
    do
        sed -i "s/$i/${i}:/g" $item
    done
done