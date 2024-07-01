#!/usr/bin/env bash

input_path=$1
output_path=$2
c=${3:-1}
g=$4

echo Input path: $input_path
echo Output path: $output_path
echo "# cpus: $c"
if [[ -z $g ]]; then
    echo Gene: $g
fi

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
python /home/scripts/get_all.py ${data_path} $(basename ${data_path} .pq)_additional.pq $c
rm lookup_ids.csv 
echo "Done retrieving additional annotations..."

echo "Merging annotations..."
python /home/scripts/merge.py ${data_path} $(basename ${data_path} .pq)_additional.pq $output_path
cp $output_path /home
echo "Annotation finished."
