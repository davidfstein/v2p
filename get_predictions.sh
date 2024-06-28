#!/usr/bin/env bash

usage() {
    echo "Parameter usage:
          -i input            path to input vcf (required)
          -o output           path to annotation output parquet file (required)
          -a annotations      absolute path to annotations (required)
          -p precomputed      path to precomputed predictions (optional)
          -c cpu              number of CPUs (optional)
          -g [gene]           annotate specific gene (optional)
          -h help             print this message and exit"
}

input_path=""
output_path=""
c=1
g=""
precomputed_path=""
annotation_path=""

while getopts ":i:o:c::g::p::a:h" opt; do
    case ${opt} in
        i )
            input_path=$OPTARG
            ;;
        o )
            output_path=$OPTARG
            ;;
        c )
            c=$OPTARG
            ;;
        g )
            g=$OPTARG
            ;;
        p )
            precomputed_path=$OPTARG
            ;;
        a )
            annotation_path=$OPTARG
            ;;
        h )
            usage
            exit 0
            ;;
        \? )
            echo "Invalid option: -$OPTARG" 1>&2
            usage
            exit 1
            ;;
        : )
            echo "Option -$OPTARG requires an argument" 1>&2
            usage
            exit 1
            ;;
    esac
done
shift $((OPTIND -1))

# Check required parameters
if [ -z "$input_path" ] || [ -z "$output_path" ] || [ -z "$annotation_path" ]; then
    echo "Error: Input, output, and annotations paths are required."
    usage
    exit 1
fi

PRECOMPUTED_TMPFILE=$(mktemp)
if [[ ${precomputed_path} ]]; then
    python ${V2P_DIR}/scripts/getPrecomputed.py $input_path ${PRECOMPUTED_TMPFILE} ${precomputed_path}
fi

docker run --rm -v $(pwd):/home -v ${annotation_path}:/cadddb dstein96/hpo -i $input_path -o $(basename $input_path .vcf)_annotations.pq -c $c -g $g

python ${V2P_DIR}/scripts/predict.py $(basename $input_path .vcf)_annotations.pq

python ${V2P_DIR}/scripts/merge.py ${PRECOMPUTED_TMPFILE} $(basename $1 .vcf)_annotations_preds.csv 

rm $PRECOMPUTED_TMPFILE
rm $(basename $1 .vcf)_annotations.pq 

mv $(basename $1 .vcf)_annotations_preds.csv ${output_path}