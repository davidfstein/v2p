#!/usr/bin/env bash

LONGOPTS=input:,output:,cpu:,gene:,help
OPTIONS=i:o:c:g:h

! PARSED=$(getopt --options=$OPTIONS --longoptions=$LONGOPTS --name "$0" -- "$@")
if [[ ${PIPESTATUS[0]} -ne 0 ]]; then
    # e.g. return value is 1
    #  then getopt has complained about wrong arguments to stdout
    echo "Please provide an input file and an output path."
    exit 2
fi
# read getopt’s output this way to handle the quoting right:
eval set -- "$PARSED"

input_path= output_path= c=1 g=
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
        -h|--help)
            echo "Parameter usage: -i|--input path/to/input.vcf (required)
                 -o|--output path/to/output.pq (required)
                 -c|--cpu # of cpus 
                 -g|--gene annotate specfic gene
                 -h|--help print this message and exit"
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

echo Input path: $input_path
echo Output path: $output_path
echo "# cpus: $c"
echo Gene: $g

mkdir /workingdir
cp $input_path /workingdir
cd /workingdir

echo "Running VEP..."
vep_out=$(basename $input_path .vcf)_vep.vcf
/home/scripts/run_vep.sh $input_path $vep_out $c $g
echo "Finished running VEP..."

echo "Processing VEP output..."
data_path=$(basename ${vep_out} .vcf).pq
bash /home/scripts/split_vep_output.sh ${vep_out} $(basename ${vep_out} .vcf).tsv header.txt $ncores
python /home/scripts/add_header.py $(basename ${vep_out} .vcf).tsv header.txt ${data_path}
python /home/scripts/get_positions_info.py ${data_path} ${data_path} ${gene}
rm ${vep_out} ${vep_out}.gz* header.txt 
echo "Finished processing VEP output..."

echo "Getting gene IDs..."
python /home/scripts/get_gene_ids.py ${data_path} lookup_ids.csv
echo "Finished getting gene IDs..."

echo "Getting protein IDs..."
python /home/scripts/get_protein_ids.py ${data_path}
echo "Finished getting protein IDs..."

echo "Mapping genomic to uniprot protein positions..."
python /home/scripts/get_uniprot_pos.py ${data_path}
echo "Finished mapping..."

echo "Retrieving additional variant annotations..."
python /home/scripts/get_all.py ${data_path} $(basename ${data_path} .pq)_additional.pq $ncores
rm lookup_ids.csv 
echo "Done retrieving additional annotations..."

echo "Merging annotations..."
python /home/scripts/merge.py ${data_path} $(basename ${data_path} .pq)_additional.pq $output_path
cp $output_path /home
echo "Annotation finished."
