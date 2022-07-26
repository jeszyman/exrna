input=$1
threads=$2
output=$(echo $3 | sed 's/.gz$//g')
log=$4

flexbar \
    --adapter-preset SmallRNA \
    --output-log $log \
    --output-reads $output \
    --pre-trim-right 1 \
    --reads $input \
    --threads $threads \
    --htrim-right AT \
    --htrim-min-length 10 \
    --htrim-error-rate 0.1 \
    --zip-output GZ
