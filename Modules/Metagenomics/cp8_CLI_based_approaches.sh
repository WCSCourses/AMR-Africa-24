
## STEP: 01

# First let's set our working directory
cd ~/course/cp8


conda activate readQC

# create output directory for fastp output
clean_reads=~/course/cp8/clean_reads
mkdir -p $clean_reads

# Provide path to the raw reads directory
raw_reads=~/course/cp8/raw_reads
mkdir -p $raw_reads

rsync ~/course/cp8/downsampled/*.gz $raw_reads/

# Execute the for loop to perform QC on all samples in the raw_reads directory
for fq in $(find $raw_reads -name "*R1.fq.gz"); do
	sampleid=$(basename -s "_R1.fq.gz" $fq)
	read1=$(find $raw_reads -name "${sampleid}*R1*f*q.gz")
	read2=$(find $raw_reads -name "${sampleid}*R2*f*q.gz")
	
	fastp -i "$read1" -I "$read2" \
	-q 20 -l 36 --cut_front -M 10 -W 4 \
	-R "$sampleid" -j $clean_reads/${sampleid}.fastp.json \
	-h $clean_reads/${sampleid}.fastp.html \
	--correction --dedup --overrepresentation_analysis --thread 4 \
	-o $clean_reads/${sampleid}.R1.fq.gz -O $clean_reads/${sampleid}.R2.fq.gz
done >> $clean_reads/qc_step1.log 2>&1



qc_reports=~/course/cp8/multiqc
multiqc -f --no-data-dir $clean_reads --outdir $qc_reports


## STEP:02

# Activate the hocort environment
conda activate hocort


# Provide path to the bowtie2 index files
bwt=~/course/cp8/databases/hocort/human


# Create directory to save decontaminated reads
hocort=~/course/cp8/hocort
mkdir -p $hocort
threads=4

for fq in $(find $clean_reads -name "*R1.fq.gz"); do
	sampleid=$(basename -s ".R1.fq.gz" $fq)
	read1=$(find $clean_reads -name "${sampleid}*R1*f*q.gz")
	read2=$(find $clean_reads -name "${sampleid}*R2*f*q.gz")
	hocort map bowtie2 --threads $threads --filter true \
	-x ${bwt}/grch38 -i $read1 $read2 \
	-o $hocort/${sampleid}.R1.fq $hocort/${sampleid}.R2.fq 2> $hocort/${sampleid}.err

	# compress reads
	gzip $hocort/${sampleid}.R1.fq
	gzip $hocort/${sampleid}.R2.fq
done



## STEP:03
# First let’s create a directory to store our databases
mkdir -p ~/course/cp8/databases/kraken2_8gb
# Next we will decompress the kraken2 database into the path we created above

tar -xvzf ~/course/cp8/databases/kraken2/k2_standard_08gb_20240112.tar.gz -C ~/course/cp8/databases/kraken2_8gb/



## STEP:04

conda activate classify
# create directory for kraken2 output
mkdir -p ~/course/cp8/kraken2

krak=~/course/cp8/kraken2


# set threads
threads=4

# set path to kraken2 database and clean_reads directory

db=~/course/cp8/databases/kraken2_8gb/
clean_reads=~/course/cp8/clean_reads


# execute kraken2
for fq in $(find $hocort -name "*R1.fq.gz"); do
	sampleid=$(basename -s ".R1.fq.gz" $fq)
	read1=$(find $hocort -name "${sampleid}*R1*f*q.gz")
	read2=$(find $hocort -name "${sampleid}*R2*f*q.gz")


	kraken2 --db "$db" --threads $threads --quick --paired \
		--output $krak/${sampleid}.kraken \
		--report $krak/${sampleid}.kraken.report \
		--memory-mapping $read1 $read2 \
		--gzip-compressed \
		--unclassified-out $krak/${sampleid}#_unclassified.fq >> $krak/krak.log
done


# STEP:05
brak=~/course/cp8/bracken
mkdir -p $brak


for file in $(find $krak -name "*kraken.report"); do 
	sampleid=$(basename -s ".kraken.report" $file)
	bracken -d "$db" -i $file -o $brak/${sampleid}.bracken.report -w ${sampleid}.bracken_species.report
done


## STEP:06
krona=~/course/cp8/krona
mkdir -p $krona


## ****** Participants working directly from the server can skip this step ******** ##
# Download and unpack the taxonomy data
copath=$(conda info --base)

cd $copath/envs/classify/opt/krona/taxonomy
wget -c ftp://ftp.ncbi.nih.gov/pub/taxonomy/taxdump.tar.gz

# prepare/build the tax data
ktUpdateTaxonomy.sh --only-build

## ******************************************************************************** ##

cd ~/course/cp8

ktImportTaxonomy -t 5 -m 3 -o $krona/grouped.krona.html $brak
conda deactivate

## STEP:07

# Now let’s get a summary of the taxonomic classification output using multiqc
conda deactivate; conda activate readQC

multiqc -f --no-data-dir $krak --outdir $qc_reports -n krackenres

conda deactivate



## STEP:08
# kma should already be available in your path, and can verify this by using any of the following
# command(s):

conda activate

kma -v
which kma

# You should either get the version of kma or the path to kma executable binary

# set threads
threads=4



# path to resfinder db
res_db=~/course/cp6/resfinder_db
kma_out=~/course/cp8/resistance


mkdir -p $kma_out


# Now let’s run kma

for fq in $(find $hocort -name "*R1.fq.gz"); do
	sampleid=$(basename -s ".R1.fq.gz" $fq)
	read1=$(find $hocort -name "${sampleid}*R1*f*q.gz")
	read2=$(find $hocort -name "${sampleid}*R2*f*q.gz")


	kma -mem_mode -ef -cge -nf -vcf -t $threads -ipe $read1 $read2 -t_db $res_db/all -o $kma_out/amr -1t1
done


