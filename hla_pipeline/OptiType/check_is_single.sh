##################################################
# File Name   : check_is_single.sh
# Author      : biolxy
# E-mail      : biolxy@aliyun.com
# Created Time: 2018年12月06日 星期四 10时40分31秒
##################################################
#!/bin/bash
# https://www.biostars.org/p/178730/
# https://bedtools.readthedocs.io/en/latest/content/tools/bamtofastq.html
check_is_single(){
    num=$(samtools view -c -f 1 $1)
    if [[ ${num} = 0 ]];then
        echo "single"
    else
        echo "paired"
    fi
}
