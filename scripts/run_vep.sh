WRKDIR=/workingdir/
cachedir=/home/.vep/
ref_fasta=/home/.vep/Homo_sapiens.GRCh38.dna.primary_assembly.fa
plugdir=/home/.vep/Plugins/

myinput=$1
myoutput=$2
ncores="${3:-1}"
gene=$4
pick="--per_gene"
if [ -z "${gene}" ]; then
	pick="--pick_allele"
fi
echo ${pick}

# HG38
vep -i $WRKDIR/$myinput -o $WRKDIR/$myoutput --dir_cache $cachedir --vcf --assembly GRCh38 ${pick} --no_stats --biotype --uniprot --symbol --numbers --canonical --protein --cache --total_length --offline --fasta $ref_fasta  --force_overwrite --fork $ncores --buffer_size 100000 --minimal --allele_number --plugin Blosum62 --plugin Conservation,$plugdir/Conservation/gerp_conservation_scores.homo_sapiens.GRCh38.bw --plugin NearestExonJB --plugin MaxEntScan,$plugdir/maxentscan/ --plugin TSSDistance --plugin dbscSNV,$plugdir/dbscSNV/dbscSNV1.1_GRCh38.txt.gz --plugin dbNSFP,$plugdir/dbNSFP/dbNSFP4.1a_grch38.gz,transcript_match=1,Eigen-PC-raw_coding,Eigen-raw_coding,GERP++_NR,GERP++_RS,GM12878_confidence_value,GM12878_fitCons_score,GenoCanyon_score,H1-hESC_confidence_value,H1-hESC_fitCons_score,HUVEC_confidence_value,HUVEC_fitCons_score,LINSIGHT,LIST-S2_score,LRT_Omega,LRT_score,MPC_score,MutationAssessor_score,SiPhy_29way_logOdds,SiPhy_29way_pi,bStatistic,integrated_confidence_value,integrated_fitCons_score,phastCons100way_vertebrate,phastCons17way_primate,phastCons30way_mammalian,phyloP100way_vertebrate,phyloP17way_primate,phyloP30way_mammalian,VEST4_score,REVEL_score,PrimateAI_score,CADD_phred,gnomAD_exomes_POPMAX_AF,gnomAD_genomes_POPMAX_AF --plugin NMD --dir_plugins $plugdir
