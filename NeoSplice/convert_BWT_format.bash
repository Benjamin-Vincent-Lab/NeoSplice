input_dir=$1
output_dir=$2

cat ${input_dir}/bwt_test.0.dat | tr "\000\001\002\003\004\005" "\$ACGNT" | msbwt convert ${output_dir}

