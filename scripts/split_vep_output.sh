infile=$1
outfile=$2
outheader=$3
ncores="${4:-1}"
echo $ncores

sed -i "/^#/s/[()]//g" $infile
sed -i "/^#/s/-/_/g" $infile
sed -i "/^#/s/++/plus_plus/g" $infile

bgzip -@ $ncores -c $infile > ${infile}.gz
tabix -p vcf ${infile}.gz

bcftools +split-vep ${infile}.gz -f '%CHROM\t%POS\t%ID\t%REF\t%ALT\t%CSQ\n' -d -A tab -o $outfile
bcftools +split-vep ${infile}.gz -l > $outheader
