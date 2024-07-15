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
    echo "Retrieving precomputed predictions..."
    python ${V2P_DIR}/scripts/getPrecomputed.py $input_path ${PRECOMPUTED_TMPFILE} ${precomputed_path}
    input_path=$(basename $input_path .vcf)_novel.vcf
fi

if [ -f  ${input_path} ] && [ $(wc -l < ${input_path}) -ge 2 ]; then
    SORT_TMPFILE=$(mktemp)
    echo -e "#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO" > ${SORT_TMPFILE}
    grep -v "^#" ${input_path} | sort -k1,1V -k2,2n >> ${SORT_TMPFILE}
    mv ${SORT_TMPFILE} ${input_path}

    echo "Annotating novel variants..."
    docker run --rm -v $(pwd):/home -v ${annotation_path}:/cadddb dstein96/hpo $input_path $(basename $input_path .vcf)_annotations.pq $c $g

    echo "Predicting novel variant impact..."
    python ${V2P_DIR}/scripts/predict.py $(basename $input_path .vcf)_annotations.pq
fi

if [ ${precomputed_path} ] && [ -f "$(basename $input_path .vcf)_annotations_preds.csv" ]; then
	echo "Merging output..."
	python ${V2P_DIR}/scripts/merge_output.py ${PRECOMPUTED_TMPFILE} $(basename $input_path .vcf)_annotations_preds.csv 
elif [[ ! -f "$(basename $input_path .vcf)_annotations_preds.csv" ]]; then
    mv ${PRECOMPUTED_TMPFILE} $(basename $input_path .vcf)_annotations_preds.csv
fi

rm -f $(basename $input_path .vcf)_annotations.pq $input_path
mv $(basename $input_path .vcf)_annotations_preds.csv ${output_path}
bash /home/david/Projects/pathophenic_back/utils/self_destroy.sh ${output_path} 86400 &
echo "Done"
