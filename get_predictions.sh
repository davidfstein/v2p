#!/usr/bin/env bash

LONGOPTS=input:,output:,cpu:,gene:,help
OPTIONS=i:o:c:g:h

! PARSED=$(getopt --options=$OPTIONS --longoptions=$LONGOPTS --name "$0" -- "$@")
if [[ ${PIPESTATUS[0]} -ne 0 ]]; then
    # e.g. return value is 1
    #  then getopt has complained about wrong arguments to stdout
    echo "Please provide input and output paths."
    exit 2
fi
# read getopt’s output this way to handle the quoting right:
eval set -- "$PARSED"

input_path= output_path= c=1 g= precomputed_path= annotation_path=
# now enjoy the options in order and nicely split until we see --
while true; do
    case "$1" in
        -i|--input)
            input_path="$2"
            shift 2
            ;;
        -o|--output)
            output_path="$2"
            shift 2
            ;;
        -c|--cpu)
            c="$2"
            shift 2
            ;;
        -g|--gene)
            g="$2"
            shift 2
            ;;
        -p|--precomputed)
            precomputed_path="$2"
            shift 2
            ;;
        -a|--annotations)
            annotation_path="$2"
            shift 2
            ;;
        -h|--help)
            echo "Parameter usage: 
                 -i|--input            path to input vcf (required)
                 -o|--output           path to annotation output parquet file (required)
                 -a|--annotations      absolute path to annotations (required)
                 -p|--precomputed      path to precomputed predictions
                 -c|--cpu              # of cpus 
                 -g|--gene             annotate specfic gene
                 -h|--help             print this message and exit"
            shift
            exit 0;
            ;;    
        --)
            shift
            break
            ;;
        *)
            echo "Programming error"
            exit 3
            ;;
    esac
done

if [[ -z $input_path ]]; then
  echo "An input file is required"
  exit 1
fi

if [[ -z $output_path ]]; then
  echo "An output path is required"
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