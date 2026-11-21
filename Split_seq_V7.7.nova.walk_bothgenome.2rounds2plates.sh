#!/bin/bash
# Author: Sai Ma <sai.ma2@mssm.edu>
# Last modified date: 2024/09/27
# Designed for processing share-atac/rna/taps/dipC/cite/cellhashing/crop

# input file
# 1) yaml
# 2) BCL or fastq or demultiplexed fastqs
# when there are 4 fastqs per lane, fastqs need to be named as "*S1_L001/2/3/4_R1/R2/I1/I2_001.fastq.gz" or "_S1_R1/R2/I1/I2_001.fastq.gz", R1: bio read1, R2: bio read2, I1: index1, I2: index2
# when there are 2 fastqs per lane, fastqs need to be named as "*S1_L001/2/3/4_R1/2_001.fastq.gz" or "_S1_R1/2_001.fastq.gz". important to rename files as R1 and "R2"
# if using demultiplexed fastq, makes sure use $project.$Type.R1/R2.fastq.gz as file name. The project and Type should match exactly the ones list in this scipt

## Example projects
Project=(Demo_DNA Demo_RNA)
Type=(DipC2 RNA2)

# ATAC RNA DipC DipC2 RNA2 crop cite cellhash notrim (keep all information, in case of any spatial design)
# use ATAC for ChIP and Cut&Tag data
# use crop for tristan guide library
# fill the same information to yaml file
# DipC2 and RNA2 are designed for 2 rounds of barcoding strategy. 
# for ATAC without taps, use ATAC, align R1 read (or both reads), disable taps
# for ATAC + taps, use ATAC, align R1 read (or both reads), enable taps

Genomes=(hg38 hg38)
# both (mm10 and hg38) mm10 guide hg38 noalign (skip alignment) 
ReadsPerBarcode=(100 100)
# reads cutoff for the filtered bam/bed: 100 for full run; 10 for QC run; 1 for crop

Start=Fastq # Bcl, Fastq, or Demultiplexed
Runtype=full # QC or full,  QC only analyze 12M reads
chem=rev # rev or fwd: nova 1.5 & nextseq use rev; nova1.0 uses fwd

###########################
## RNA-seq options, usually don't change these options
removeSingelReadUMI=F # default F; T if seq deep enough or slant is a big concern; def use F for crop-seq
keepIntron=T #default T
cores=32
genename=gene_name # defaultgene_name; gene_name (official gene symbol) or gene_id (ensemble gene name)
refgene=gencode # default gencode; gencode or genes; genes is UCSC genes; gencode also annotate ncRNA
mode=fast # fast or regular; default regular; fast: dedup with my own script; regular: dedup with umitools
# fast mode gives more UMIs because taking genome position into account when dedup
# fast mode doesn't collapse UMIs map to different position, the dup rate for RNA is not accurate

keepMultiMapping=(F F)
# default F; F for species mixing or cell lines, T for low yield tissues (only keep the primarily aligned reads)
# does not matter for DNA assays
# allowing multi-mapping redas will increase the percent of mito reads

N6=(T T)
# default F, T is random primer is used; N6 is not comptiable with regular mode and removeSingelReadUMI mode

#### DipC option:
## Resolutions=1000000 or Resolutions=1000000,2000000,3000000
Resolutions=10000,20000,50000,100000,500000,1000000
BaseResolution=$(echo $Resolutions | awk -F, '{print $1}')

## DNA options
oneread=(F F F F F F F F F F F F F F F F F F F F F F F F F F F F F F F F F F F F F F F F F F \
	   F F F F F F F F F F F F F F F F F F F F F F F \
	   F F F F F F F F F F F F F F F F F F F F F F F) # this option needs to be updated
# option F, Sh, R1; default R1; for DipC2, has to be F
# R1: when R2 is too short, only read1 is used for alignment
# Sh: select shortest fragment using both R1 and R2. only use this if R2 is long enough for ATAC

taps=(F F)
# option F, T. If T, use bwa-meth to align and extract DNA methylation bed file

cleanup=F
# clean up fastq or not

# define path
## User configuration
basedir="${BASEDIR:-/path/to/Split-seq}"

rawdir="${RAWDIR:-${basedir}/demo_fastq}"
dir="${OUTDIR:-${basedir}/output}"
yaml="${CONFIG_YAML:-${basedir}/config_demo.yaml}"

## Pipeline scripts
myPATH="${PIPELINE_SCRIPTS:-${basedir}/scripts/SHARE-seq-nova}"

## References
refdir="${REFERENCE_DIR:-${basedir}/references}"
genomedir="${GENOME_DIR:-${basedir}/Genome}"

tssFilesPATH="${refdir}/TSSfiles"
genomeBed="${refdir}/genomeBed"
gtf="${refdir}/gtf"

bowtieGenome="${genomedir}/bowtie2"
starGenome="${genomedir}/star"
bwaGenome="${genomedir}/bwa"
bwaMem2Genome="${genomedir}/bwa-mem2"
biscuitGenome="${genomedir}/biscuit"

picardPATH="${PICARD_JAR:-${basedir}/picard.jar}"
juicerToolsPATH="${JUICER_TOOLS_JAR:-${basedir}/tools/juicer_tools_1.22.01.jar}"

## For crop option only.
cropGenome="${CROP_GENOME:-${bowtieGenome}/crop/crop}"

export biscuitGenome
################### DO NOT CHANGE BELOW ####################
## version log
# v36 add option for N6 RT primer; remove "matchPolyT" option
# v35 add supoort for saving RNA matrix as h5, the same as 10x v3, this can be used for Running Cellbender
# v4.0 try to speed up pipeline
## remove R1 level barcode
## rewrite ATAC pipeline and demultiplexing, 3x faster
## automatically detect seq chemistry
## removed Fastq_SplitLane
## autimatically detect reads are indexed or not
## try to process RNA in bed format as well
## always remove GGGGGG reads
## removed SkipPlotUMIgene
# v4.1
## remove RawReadsPerBarcode
# v4.3
## add support for multiple lanes
## avoid sorting bam with picard
## removed N6 option since it does not really improve no. of genes
# v4.7
## add fast mode for faster RNA dedup
## sorting with samtools instead of picard
## use pigz to spped up gzip
# v4.8
## add back cutoff for fragment file
# v4.9
## update support for crop-seq
# modify UMI_gene_perCell_plot_v3.R to use tidyR and pivot_wider instead of acast
# v4.10
## add flag to stop alignment if peeps perfer to align by themselves
# v4.11
## changed STAR setting to improve unique alignment rate: removed outFilterMismatchNoverLmax 0.06; added outFilterScoreMinOverLread 0.3 outFilterMatchNminOverLread 0.3
## fixed bug in fast mode, now RNA tsv to bed is correct
# v4.12
## updated R path to lab R
## add notrim option to fastqs
# v4.13
## integrste align tag to Genome (noalign)
## fixed TAPs DNA methylation calling bug
## change no_overlap to --include_overlap, to keep more data in taps
# v5.1
## expand to 192 hybrid oligos
## changed RNA/crop/cite/cellhash demultiplexed-fastq format from readname_UMI_barcode to readname_barcode_UMI
# v5.2
## fixed bug in fast pipeline where i messed up gene names with "_" by using: awk 'gsub(/[_]/,"\t", $1)'
## fixed bug in rm_dup_barcode_UMI_v3.py, the cell barcode was shifted...
## disabled 1bp barcode shift in fastq demultiplexing, which results in mis-assignment of 83 to 143
## fixed bug in Tn5 bias shifting +4/-5 (it was +4/+4)
## fixed bug in lib_size_sc_V5_single_species.R, when calculting ATAC lib size, should not divided by 2
# v 5.3
## add support for calculate lib size of atac with only R2
# v5.4
## updated fastq naming fashing from R1-R4 to R1/R2/I1/I2
## allowing splited lane from bcl2fastq
# v5.5
## add option to keep mito reads in ATAC data, which is useful for footprinting
# v6.0
## add support for random priming
## add support to start with demultiplexed fastq
## remove the discarded fastq option
# v6.2
## change N6 alignment to R1/R2 or both
## add NNNNNN to N6 reads
## update species mixing human genome from hg19 to hg38
## fixed bug in the UMI count generation (acast funtion convert all >1 UMIs to 1)
# v6.3
## add support for dipC
# 6.4
## add fragment file for taps
# 6.5/6.6
## adjust pairtools for dedup
## adjust taps trimming based on the overlap 
# 6.8
## add stats for hic
## use R1 only for taps, since R2 is too short after removing adapter
## switched R1 and R3 barcode annotation
#6.9
## trim taps reads as default
## trim ATAC reads using trim gloare, becasue 150p is too long
#6.10
## fix bug in calcualte lib size for taps with one read only
## fix bug of filtering out singleRead UMIs
#6.11
## fix bug for taps barcode
## change dipC alignment to bowtie local mode
## fix bug using R1 only for ATAC. both reads are aligned, but only R1 is used for making fragments
# 6.12
## align hic by first aligning stiteched reads and unstitched reads separately, then convert both to pairs, then merge them
# 6.14 remove bed.gz for hiC
# 6.15 adjust to 2 rounds of ligation and remove read stitching
# 6.19 use bwameth for dipC2+taps and biscuit to extract methylation
# 7.1 try bwameth for ATAC + taps
## to do: remove rigid.reheader bams
## try use combined genome for ATAC mapping
# 7.2 fixed missing chrososomes in RNA species mixing bug
# 7.4 change bwameth to biscuit, becuase non-directional alignment is needed
# 7.5 remove keepmito option. always keep mito for methylation calling, but remove mito for ATAC fragments
# 7.6 change pairtools setting using walk-policy all
# 7.7 align hiC to combined genome; also generate non-walk pairs
# 2rounds2plates: add barcode 97-192 support in fastq.process.py3

## todo: add mCH option for bsmap bam
## todo: de-multiplexing per barcode

# real code
export SHELL=$(type -p bash)
#source ~/.bashrc
export LC_COLLATE=C
export LANG=C
export TMPDIR=$dir/fastqs
## this may speed up sorting by limit output to bytes 

echo "the number of projects is" ${#Project[@]}
echo "Running $Runtype pipeline"

if [ ! -d $dir ]; then mkdir $dir; fi
if [ ! -d $dir/fastqs ]; then mkdir $dir/fastqs ; fi
if [ ! -d $dir/temp ]; then mkdir $dir/temp ; fi

cp "$0" $dir/
cp $yaml $dir/
cd $dir 
if [ -f $dir/Run.log ]; then rm $dir/Run.log; fi

# Use executables from the active Conda environment or system PATH

# if start with bcl file
if [ "$Start" = Bcl ]; then
    echo "Bcl2fastq"
    if [ -f $rawdir/fastqs/Undetermined_S0_R1_001.fastq.gz ] || [ -f $rawdir/fastqs/Undetermined_S1_R1_001.fastq.gz ] || [ -f $rawdir/fastqs/Undetermined_S0_L001_R1_001.fastq.gz ]; then
	echo "Found Undetermined_S0_L001_I1_001.fastq.gz, skip Bcl2Fastq"
    else
	echo "Converting bcl to fastq"
	mkdir $rawdir/fastqs/
	bcl2fastq -p $cores -R $rawdir --mask-short-adapter-reads 0 -o $rawdir/fastqs/ --use-bases-mask Y*,I8n*,I8n*,Y* --create-fastq-for-index-reads  2>>$dir/Run.log
	cd $rawdir/fastqs/    
    fi
    if [ -f $dir/fastqs/filesi1.xls ]; then rm $dir/fastqs/filler*; fi
    
    rawdir=$rawdir/fastqs/
    Start=Fastq
fi

if [ "$Start" = Demultiplexed ]; then
    mkdir $dir/fastqs/
    noProj=${#Project[@]}
    echo "link fastqs"
    for (( i=0; i< $noProj; i++ )); do
	if [ -f $rawdir/${Project[$i]}.R1.fastq.gz  ]; then 
	    echo "link" ${Project[$i]} fastqs
	    ln -s $rawdir/${Project[$i]}.R1.fastq.gz $dir/fastqs/${Project[$i]}.R1.fastq.gz
	    ln -s $rawdir/${Project[$i]}.R2.fastq.gz $dir/fastqs/${Project[$i]}.R2.fastq.gz
	else
	    echo "Error...couldn't find" ${Project[$i]}.R1.fastq.gz
	fi
    done
fi

if [ "$Start" = Fastq ]; then
    echo "Skip bcltofastq"
    if [ -f $dir/fastqs/filesi1.xls ]; then rm $dir/fastqs/filler*; fi
    cd $rawdir
    if ls *L001_R1_001.fastq.gz 1> /dev/null 2>&1; then
	temp=$(ls *_L001_R1_001.fastq.gz)
	Run=$(echo $temp | sed -e 's/\_S0\_L001\_R1\_001.fastq.gz//')
	singlelane=F
	temp=$(ls *L00*_R1_001.fastq.gz)
	VAR=( $temp )
	nolane=${#VAR[@]}
	echo "Detected $nolane lanes"
    elif ls *S1_R1_001.fastq.gz 1> /dev/null 2>&1; then
	echo "Detected single lane"
	temp=$(ls *S1_R1_001.fastq.gz)
	Run=$(echo $temp | sed -e 's/\_\S1\_\R1\_\001.fastq.gz//')
	singlelane=T
	nolane=1
    else
	echo "No fastq with matched naming format detected; exit..."
	exit
    fi
    
    echo "Run number is:" $Run

    # split fastqs
    mkdir $dir/smallfastqs/
    if [ -f $dir/smallfastqs/0001.1.$Run.R2.fastq.gz ]; then
        echo "Found 0001.$Run.R2.fastq, skip split fastqs"
    else
	if [ ! -f $dir/R1.1.fastq.gz ]; then
            echo "Link fastqs"
	    if [ $singlelane == T ]; then
		ln -s $rawdir/"$Run"_S1_R1_001.fastq.gz $dir/R1.1.fastq.gz
		ln -s $rawdir/"$Run"_S1_R2_001.fastq.gz $dir/R2.1.fastq.gz
		ln -s $rawdir/"$Run"_S1_I1_001.fastq.gz $dir/I1.1.fastq.gz
		ln -s $rawdir/"$Run"_S1_I2_001.fastq.gz $dir/I2.1.fastq.gz
	    else
		parallel 'ln -s '$rawdir'/'$Run'_S0_L00{}_R1_001.fastq.gz '$dir'/R1.{}.fastq.gz' ::: $(seq $nolane)
		parallel 'ln -s '$rawdir'/'$Run'_S0_L00{}_R2_001.fastq.gz '$dir'/R2.{}.fastq.gz' ::: $(seq $nolane)
		parallel 'ln -s '$rawdir'/'$Run'_S0_L00{}_I1_001.fastq.gz '$dir'/I1.{}.fastq.gz' ::: $(seq $nolane)
		parallel 'ln -s '$rawdir'/'$Run'_S0_L00{}_I2_001.fastq.gz '$dir'/I2.{}.fastq.gz' ::: $(seq $nolane)
	    fi
	fi
	
	if [ "$Runtype" = full ]; then
	    # Runing full pipeline
	    echo "Split fastqs to small files"
	    dosplitfull(){
                fastp -i $2/R1.$1.fastq.gz -o $2/smallfastqs/$1.$3.R1.fastq.gz -S 80000000 --thread 1 -d 4 -A -G -L -Q 2>>$2/split.log &
                fastp -i $2/R2.$1.fastq.gz -o $2/smallfastqs/$1.$3.R2.fastq.gz -S 80000000 --thread 1 -d 4 -A -G -L -Q 2>>$2/split.log &
                fastp -i $2/I1.$1.fastq.gz -o $2/smallfastqs/$1.$3.I1.fastq.gz -S 80000000 --thread 1 -d 4 -A -G -L -Q 2>>$2/split.log &
                fastp -i $2/I2.$1.fastq.gz -o $2/smallfastqs/$1.$3.I2.fastq.gz -S 80000000 --thread 1 -d 4 -A -G -L -Q 2>>$2/split.log &
		wait
            }
            export -f dosplitfull
            parallel --delay 1 dosplitfull {} $dir $Run ::: $(seq $nolane)
	elif [ "$Runtype" = QC ]; then
	    # Runing QC pipeline
	    echo "Split fastqs to small files"
	    dosplitQC(){
		fastp -i $2/R1.$1.fastq.gz -o $2/smallfastqs/$1.$3.R1.fastq.gz --reads_to_process $4 -S 4000000 --thread 1 -d 4 -A -G -L -Q 2>>$2/split.log &
		fastp -i $2/R2.$1.fastq.gz -o $2/smallfastqs/$1.$3.R2.fastq.gz --reads_to_process $4 -S 4000000 --thread 1 -d 4 -A -G -L -Q 2>>$2/split.log &
		fastp -i $2/I1.$1.fastq.gz -o $2/smallfastqs/$1.$3.I1.fastq.gz --reads_to_process $4 -S 4000000 --thread 1 -d 4 -A -G -L -Q 2>>$2/split.log &
		fastp -i $2/I2.$1.fastq.gz -o $2/smallfastqs/$1.$3.I2.fastq.gz --reads_to_process $4 -S 4000000 --thread 1 -d 4 -A -G -L -Q 2>>$2/split.log &
		wait
            }
            export -f dosplitQC
	    let reads=12100000/$nolane
	    parallel --delay 1 dosplitQC {} $dir $Run $reads ::: $(seq $nolane)
	else
            echo "Unknown sequencer type, exiting" && exit
        fi
    fi
    
    # trim and index fastq
    if [ -f $dir/fastp.json ]; then rm $dir/fastp.json $dir/fastp.html; fi
    ls $dir/smallfastqs | grep R1 > $dir/filesr1.xls
    ls $dir/smallfastqs | grep R2 > $dir/filesr2.xls
    ls $dir/smallfastqs | grep I1 > $dir/filesi1.xls
    ls $dir/smallfastqs | grep I2 > $dir/filesi2.xls
    cd $dir/
    
    if [ -f $dir/fastqs/Sub.0001.1.$Project.R1.fq.gz ] || [ -f $dir/fastqs/$Project.R1.fastq.gz ]; then
        echo "Found fq.gz, skip updating index"
    else
        echo "Update index and trim fastqs"
	noreadfile=`ls $dir/smallfastqs | grep R1 2>/dev/null | wc -l`
	noindexfile=`ls $dir/smallfastqs | grep I1 2>/dev/null | wc -l`
	if [ $noreadfile == $noindexfile ]; then
	    paste filesr1.xls filesr2.xls filesi1.xls filesi2.xls | awk -v OFS='\t' '{print $1, $2, $3, $4, substr($1,1,7)}'> Filelist2.xls
	    parallel --jobs $cores --colsep '\t' 'if [ -f '$dir'/fastqs/Sub.{5}.'$Project'.R1.fq.gz ]; then echo "found Sub.{5}.'$Project'.R1.fq.gz"; \
	    	     	    	   	  else python3 '$myPATH'/fastq.process.py3.v3.2plates.py \
                                          -a '$dir'/smallfastqs/{1} -b '$dir'/smallfastqs/{2} \
                                          --c '$dir'/smallfastqs/{3} --d '$dir'/smallfastqs/{4} \
					  --out '$dir'/fastqs/Sub.{5} \
					  -t '$chem' -y '$yaml' && pigz --fast -p 4 '$dir'/fastqs/Sub.{5}*fq; fi' :::: Filelist2.xls
	else
	    paste filesr1.xls filesr2.xls | awk -v OFS='\t' '{print $1, $2, substr($1,1,7)}'> Filelist2.xls
	    parallel --jobs $cores --colsep '\t' 'if [ -f '$dir'/fastqs/Sub.{3}.'$Project'.R1.fq.gz ]; then echo "found Sub.{3}.'$Project'.R1.fq.gz"; \ 
	    	     	    	   	  else python3 '$myPATH'/fastq.process.py3.v3.2plates.py -a '$dir'/smallfastqs/{1} -b '$dir'/smallfastqs/{2} \
                                          --out '$dir'/fastqs/Sub.{3} \
					  -t '$chem' -y '$yaml' && pigz --fast -p 4 '$dir'/fastqs/Sub.{3}*fq; fi' :::: Filelist2.xls
	    # --qc # option is available to process even smaller number of reads
	fi
    fi
    rm filesr1.xls filesr2.xls filesi1.xls filesi2.xls
    if [ -f Filelist2.xls ]; then
	rm Filelist2.xls
    fi
fi 

# merge fastq
echo "Merge fastqs"
parallel --jobs 4 'if [ -f '$dir'/fastqs/{}.R1.fastq.gz ] || [ -d '$dir'/{} ]; then echo "found {}.R1.fastq.gz or {} folder"; \
	 	      	   else ls '$dir'/fastqs/Sub*{}*R1.fq.gz | xargs cat > '$dir'/fastqs/{}.R1.fastq.gz && echo "Generated {}.R1.fastq.gz"; fi' ::: ${Project[@]}  
parallel --jobs 4 'if [ -f '$dir'/fastqs/{}.R2.fastq.gz ] || [ -d '$dir'/{} ]; then echo "found {}.R2.fastq.gz or {} folder"; \
                           else ls '$dir'/fastqs/Sub*{}*R2.fq.gz | xargs cat > '$dir'/fastqs/{}.R2.fastq.gz && echo "Generated {}.R2.fastq.gz"; fi' ::: ${Project[@]} 

## clean up small fastqs
if [ $cleanup == "T" ]; then
    rm -r $dir/smallfastqs/* $dir/*/Sub.*.fq.gz
    touch $dir/smallfastqs/0001.1.$Run.R2.fastq.gz
fi


# align
index=0
for Name in ${Project[@]}; do
    echo "project $index : $Name"
    
    if [ ${Type[$index]} == "notrim" ]; then
	if [ ! -d $dir/$Name ]; then
	    mkdir $dir/$Name
	fi
    fi
    if [ -d $dir/$Name ]; then
	echo "Found $Name dir, skip this project"
    else
	
	if [ ${Genomes[$index]} == "noalign" ]; then
	    mkdir $dir/$Name
	    mv $dir/fastqs/$Name*.fastq.gz $dir/$Name
	    ## this will skip aligning this project
	fi
	if [ ${Type[$index]} == "DipC2" ]; then
	    Type[$index]="DipC"
	fi
	if [ ${Type[$index]} == "RNA2" ]; then
	    Type[$index]="RNA"
        fi   
	    #if [ ${Type[$index]} == "DipC" ]; then
	#	Genome1=(hg38 mm10) # Genome1 for alignment, Genome2 for other steps
	#    else
	#	Genome1=(both) # for RNA and cellhash
	    #    fi
	if [ ${Genomes[$index]} == "both" ]; then
	    Genome1=(both)
            Genome2=(hg38 mm10)
	else
	    Genome1=${Genomes[$index]}
	    Genome2=${Genomes[$index]}
	fi
	
	for Species in ${Genome1[@]}; do
	    if [ -d $dir/$Name ]; then
		echo "Found $Name folder, skip alignment"
	    else
		cd $dir/fastqs/
		
		echo "Align $Name to $Species ${Type[$index]} library"
		if [ -f $dir/fastqs/$Name.$Species.bam ] || [ -f $dir/fastqs/$Name.$Species.st.bam ]; then
		    echo "Found $Name.$Species.bam, skip alignment"
		elif [ ${Type[$index]} == "ATAC" ]; then
		    if [ -f $dir/fastqs/$Name.R1.trimmed.fq.gz ]; then
			echo "found $Name.R1.trimmed.fq.gz"
		    else
			echo "Trim fastqs"
			(trim_galore -j $cores --length 15 --paired $dir/fastqs/$Name.R1.fastq.gz $dir/fastqs/$Name.R2.fastq.gz) 2>$Name.$Species.trim.log
			mv $dir/fastqs/$Name.R1_val_1.fq.gz $dir/fastqs/$Name.R1.trimmed.fq.gz
			mv $dir/fastqs/$Name.R2_val_2.fq.gz $dir/fastqs/$Name.R2.trimmed.fq.gz
		    fi
		    if [ ${oneread[$index]} == F ] || [ ${oneread[$index]} == "Sh" ]; then
			if [ ${taps[$index]} == "F" ]; then
			    echo "Align reads with bowtie2"
			    (bowtie2 -X2000 -p $cores --rg-id $Name \
				     -k 5 \
				     -x $bowtieGenome/$Species/$Species \
				     -1 $dir/fastqs/$Name.R1.trimmed.fq.gz \
				     -2 $dir/fastqs/$Name.R2.trimmed.fq.gz | \
				 samtools view -bS -@ $cores - -o $Name.$Species.bam) 2>$Name.$Species.align.log
			else
			    echo "Align reads with biscuit undirectional"
			    (biscuit align -t $cores $biscuitGenome/$Species.fa $Name.R1.trimmed.fq.gz $Name.R2.trimmed.fq.gz | \
				 samtools view -@ 8 -bS > $dir/fastqs/$Name.$Species.bam) 2>$Name.$Species.align.log
			    ## filter -q 1 gives idential results as bismark
			fi
		    elif [ ${oneread[$index]} == "R1" ]; then
			if [ ${taps[$index]} == "F" ]; then
			    echo "Align reads with bowtie2"
                            (bowtie2 -X2000 -p $cores --rg-id $Name \
				     -k 5 \
				     -x $bowtieGenome/$Species/$Species \
				     -U $dir/fastqs/$Name.R1.trimmed.fq.gz | \
				 samtools view -bS -@ $cores - -o $Name.$Species.bam) 2>$Name.$Species.align.log
			    # relax the definitation of unique mapper, now keeps up to 5 multi-mappers
			else
			    echo "Align reads with biscuit undirectional"
			    (biscuit align -t $cores $biscuitGenome/$Species.fa $Name.R1.trimmed.fq.gz | \
                        	 samtools view -@ 8 -bS > $dir/fastqs/$Name.$Species.bam) 2>$Name.$Species.align.log
			    ## filter -q 1 gives almost idential results as bismark
			fi
		    fi
		elif [ ${Type[$index]} == "crop" ]; then
		    # -x $bowtieGenome/viralRNA/viralRNA
                    (bowtie2 -X2000 -p $cores --rg-id $Name \
			     -N 1 -L 20 \
                             -x "$cropGenome" \
			     -U $dir/fastqs/$Name.R1.fastq.gz | \
                         samtools view -bS -@ $cores - -o $Name.$Species.bam) 2>$Name.$Species.align.log
		elif [ ${Type[$index]} == "cite" ]; then
                    (bowtie2 -X2000 -p $cores --rg-id $Name \
                             -N 0 -L 30 \
                             -x $bowtieGenome/cite/cite \
                             -U $dir/fastqs/$Name.R1.fastq.gz | \
                         samtools view -bS -@ $cores - -o $Name.$Species.bam) 2>$Name.$Species.align.log
		elif [ ${Type[$index]} == "cellhash" ]; then
                    echo "Match CITE-seq barcodes, convert fastq to bam"
                    (bowtie2 -p $cores --rg-id $Name \
                             -x $bowtieGenome/cellhash/$Species.cellhash/$Species.cellhash \
                             -U $dir/fastqs/$Name.R1.fastq.gz | \
                         samtools view -bS - -o $Name.$Species.bam) 2>$Name.$Species.align.log
		elif [ ${Type[$index]} == "DipC" ]; then
		     if [ -f $dir/fastqs/$Name.R1.trimmed.fq.gz ]; then
			echo "found $Name.R1.trimmed.fq.gz"
		    else
			(trim_galore -j $cores --length 15 --paired $dir/fastqs/$Name.R1.fastq.gz $dir/fastqs/$Name.R2.fastq.gz) 2>$Name.$Species.trim.log
			mv $dir/fastqs/$Name.R1_val_1.fq.gz $dir/fastqs/$Name.R1.trimmed.fq.gz
			mv $dir/fastqs/$Name.R2_val_2.fq.gz $dir/fastqs/$Name.R2.trimmed.fq.gz
		    fi
		    if [ ${taps[$index]} == "F" ]; then
			echo "Align reads with bwa-mem2"
			(bwa-mem2 mem -t $cores -SP5M $bwaMem2Genome/$Species.fa $Name.R1.trimmed.fq.gz $Name.R2.trimmed.fq.gz | samtools view -@ $cores -bS -o $dir/fastqs/$Name.$Species.bam) 2>$Name.$Species.align.log
		    else
			echo "Align reads with biscuit directional"
			(biscuit align -t $cores -SP5M -b 1 $biscuitGenome/$Species.fa $Name.R1.trimmed.fq.gz $Name.R2.trimmed.fq.gz | \
                                samtools view -bS -@ 8 -o $dir/fastqs/$Name.$Species.bam) 2>$Name.$Species.align.log
		    fi
		elif [ ${Type[$index]} == "RNA" ]; then
		    if [ -f $dir/fastqs/$Name.$Species.align.log ]; then
			echo "found $dir/fastqs/$Name.$Species.align.log, skip alignment for UMI reads"
		    else
			echo "Align UMI reads"
			STAR --chimOutType WithinBAM \
			     --runThreadN $cores \
			     --genomeDir $starGenome/$Species/ \
			     --readFilesIn $dir/fastqs/$Name.R1.fastq.gz  \
			     --outFileNamePrefix $dir/fastqs/$Name.$Species. \
			     --outFilterMultimapNmax 20 \
			     --outFilterScoreMinOverLread 0.3 \
			     --outFilterMatchNminOverLread 0.3 \
			     --outSAMattributes NH HI AS nM MD \
			     --limitOutSJcollapsed 5000000 \
			     --outSAMtype BAM Unsorted \
			     --limitIObufferSize 400000000 \
			     --outReadsUnmapped None \
			     --readFilesCommand zcat

			if [ ${N6[$index]} == "F" ]; then
			    mv $dir/fastqs/$Name.$Species.Aligned.out.bam $dir/fastqs/$Name.$Species.bam
			else
			    echo "Add PT tag to label N6"
			    samtools view -H $dir/fastqs/$Name.$Species.Aligned.out.bam > $dir/fastqs/$Name.$Species.header.sam
			    ## add PT =1 tag to dT reads; PT = 0 to N6 reads
			    cat $dir/fastqs/$Name.$Species.header.sam <(samtools view -@ $cores $dir/fastqs/$Name.$Species.Aligned.out.bam | \
									    awk -v OFS='\t' '{if ($1 ~ /NNNNNNNNNN/) print $0, "PT:i:0"; else print $0, "PT:i:1"}') \
				| samtools view -@ $cores -bS > $dir/fastqs/$Name.$Species.bam
			    rm  $dir/fastqs/$Name.$Species.header.sam $dir/fastqs/$Name.$Species.Aligned.out.bam				
			fi
                        rm -r *_STARtmp *Log.progress.out *SJ.out.tab
			mv $dir/fastqs/$Name.$Species.Log.final.out $dir/fastqs/$Name.$Species.align.log
		    fi
		else
		    echo "Unknown parameter"
		    exit
		fi
		if [ -f $dir/fastqs/$Name.$Species.st.bam ]; then
		    echo "Found $Name.$Species.st.bam, skip sorting bam"
		else
		    echo "Sort $Name.$Species.bam"
		    cd $dir/fastqs/
		    samtools sort -@ 8 -m 8G $Name.$Species.bam > $Name.$Species.st.bam
		    samtools index -@ $cores $Name.$Species.st.bam
		    if [ ${Type[$index]} != "DipC" ]; then
			rm $Name.$Species.bam
		    fi
		fi

		# Update RGID's and Headers
		if [ -f $dir/fastqs/$Name.$Species.unique.st.bam.bai ]; then
		    echo "Found $Name.$Species.unique.st.bam, skip update RGID"
		elif [ ${Type[$index]} == "RNA" ] || [ ${Type[$index]} == "crop" ] || [ ${Type[$index]} == "cite" ] || [ ${Type[$index]} == "cellhash" ]; then
		    if [ -f $Name.$Species.unique.st.bam.bai ]; then
                        echo "Found $Name.$Species.unique.st.bam, skip update RG tag"
		    else
			echo "Update RGID for $Name.$Species.st.bam"
			echo "Remove low quality reads and unwanted chrs"
			if [ ${keepMultiMapping[$index]} == "T" ] && [ ${Type[$index]} == "RNA" ]; then
			    # keep primary aligned reads only
			    # modify chrMT to chrM
			    samtools view -@ $cores -bS -F 256 $Name.$Species.st.bam > $Name.$Species.unique.st.bam
			else
                            samtools view -@ $cores -bS -q 30 $Name.$Species.st.bam  > $Name.$Species.unique.st.bam
			fi
			samtools index -@ $cores $Name.$Species.unique.st.bam
		    fi
		fi

		if [ ${Type[$index]} == "ATAC" ] || [ ${Type[$index]} == "DipC" ]; then
		    if [ ${Genomes[$index]} == "both" ] && [ ! -f $Name.mm10.st.bam ]; then
			echo "Split into hg and mm"
			chrs1=`samtools view -H $Name.$Species.st.bam | grep GRCh38 | cut -f2 | sed 's/SN://g' | awk '{if(length($0)<13)print}'`
			chrs2=`samtools view -H $Name.$Species.st.bam | grep mm10 | cut -f2 | sed 's/SN://g' | awk '{if(length($0)<13)print}'`
			samtools view -@ $cores -h $Name.$Species.st.bam `echo ${chrs1[@]}` | sed 's/GRCh38_//g' | samtools view -@ $cores -b -o $Name.hg38.st.bam &
			samtools view -@ $cores -h $Name.$Species.st.bam `echo ${chrs2[@]}` | sed 's/mm10___//g' | sed 's/mm10_//g' | samtools view -@ $cores -b -o $Name.mm10.st.bam &
			wait
			samtools index -@ $cores $Name.hg38.st.bam &
			samtools index -@ $cores $Name.mm10.st.bam &
			wait
		    fi
		    if [ ${Type[$index]} == "DipC" ] && [ ${Genomes[$index]} == "both" ] && [ ! -f $Name.mm10.bam ]; then
			samtools sort -n -@ 8 -m 8G $Name.hg38.st.bam > $Name.hg38.bam
			samtools sort -n -@ 8 -m 8G $Name.mm10.st.bam > $Name.mm10.bam
		    fi
		fi
		# ATAC processing convert bam to bed
		if [ ${Type[$index]} == "ATAC" ]; then
		    for Species2 in ${Genome2[@]}; do
			if [ -f $Name.$Species2.rmdup.bam ]; then
			    echo "Skip processing bam"
			else
			    if [ -f $Name.$Species2.namesort.bam ]; then
				echo "Skip sort bam on name"
			    else
				echo "Remove low quality reads, unwanted chrs & namesort" $Name.$Species2.st.bam
				chrs=`samtools view -H $Name.$Species2.st.bam | grep chr | cut -f2 | sed 's/SN://g' | grep -v chrM | grep -v Y | awk '{if(length($0)<6)print}'`
				if [ ${oneread[$index]} == F ] || [ ${oneread[$index]} == "Sh" ]; then					    
				    samtools view -b -F 256 -f 0x2 $Name.$Species2.st.bam  `echo $chrs` | samtools sort -@ $cores -m 3G -n -o $Name.$Species2.namesort.bam
				else
                             	    samtools view -b -F 256 -@ $cores $Name.$Species2.st.bam  `echo $chrs` -o $Name.$Species2.namesort.bam
				fi
			    fi
			    if [ -f $Name.$Species2.bed.gz ]; then
				echo "Skip converting bam to bed.gz"
			    else
				echo "Convert namesort.bam to bed.gz & mark duplicates"
				if [ ${oneread[$index]} == F ]; then
				    (bedtools bamtobed -i $Name.$Species2.namesort.bam -bedpe | \
					 sed 's/_/\t/g' | \
					 awk -v OFS="\t" '{if($10=="+"){print $1,$2+4,$6-5,$8}else if($10=="-"){print $1,$2-5,$6+4,$8}}' |\
					 sort --parallel=$cores -S 40G  -k4,4 -k1,1 -k2,2n -k3,3n | \
					 uniq -c | \
					 awk -v OFS="\t" '{print $2, $3, $4, $5, $1}' | \
					 pigz --fast -p $cores > $Name.$Species2.bed.gz)2> temp.log
				    rm temp.log $Name.$Species2.namesort.bam
				elif [ ${oneread[$index]} == "Sh" ]; then
	                            bedtools bamtobed -i $Name.$Species2.namesort.bam -bedpe > $Name.$Species2.bedpe
				    # negative sign indicates the first coordinate is from R2; + indicates the second is from R2
				    # only use 5' on R2 as unique location
				    bedtools bamtobed -i $Name.$Species2.namesort.bam |  awk 'NR%2 == 1 {print $6}' > $Name.$Species2.sign.bed
				    paste $Name.$Species2.bedpe $Name.$Species2.sign.bed | \
					sed 's/_/\t/g' | \
					sed 's/\/1//g' | \
					sed 's/\/2//g' | \
					awk -v OFS="\t" '{if($12=="-") print $1, $2+4, $6-5, $8, $6-$2; else print $1, $6-5, $2+4, $8, $6-$2}' | \
					sort --parallel=$cores -S 40G  -k4,4 -k1,1 -k2,2n | \
					awk -v OFS="\t" 'NR==1 { t1=$1;t2=$2;end=$3; bc=$4;insert=$5; sum=0} {if(t1==$1 && t2==$2 && bc==$4) {sum+=1; if(insert<$5) end=$3; insert=$6} else {print t1, t2, end, bc, sum; t1=$1;t2=$2;end=$3;sum=1;bc=$4; insert=$5}} END {print t1, t2, end, bc, sum}' | \
					awk -v OFS="\t" '{if($2>$3) print $1, $3, $2, $4, $5; else print}' | \
					pigz --fast -p $cores > $Name.$Species2.bed.gz
				    rm $Name.$Species2.bedpe $Name.$Species2.sign.bed $Name.$Species2.namesort.bam
				elif [ ${oneread[$index]} == "R1" ]; then
                                    bedtools bamtobed -i $Name.$Species2.namesort.bam | \
					sed 's/_/\t/g' | \
					sed 's/\/1//g' | \
					awk -v OFS="\t" '{if($7=="+"){print $1,$2+4,$3+4,$5}else if($7=="-"){print $1,$2-5,$3-5,$5}}' | \
					sort --parallel=$cores -S 40G  -k4,4 -k1,1 -k2,2n -k3,3n | \
					uniq -c | \
					awk -v OFS="\t" '{print $2, $3, $4, $5, $1}' | \
					pigz --fast -p $cores > $Name.$Species2.bed.gz
                                    rm $Name.$Species2.namesort.bam
				else
				    echo "Unkown option for deduplicating, exiting..."
				    exit
				fi
			    fi
			    # convert to a bam file for QC
			    bedToBam -i <(zcat $Name.$Species2.bed.gz | \
					      awk -v OFS="\t" '{if($2>$3) print $1, $3, $2, $4, $5; else print}') -g $genomeBed/$Species2.chrom.sizes | \
				samtools sort -@ $cores -m 2G - > $Name.$Species2.rmdup.bam
			    samtools index -@ $cores $Name.$Species2.rmdup.bam
			 fi
			 # count reads for ATAC
			 if [ -f $Name.$Species2.filtered.counts.csv ] || [ -f $Name.$Species2.filtered.counts.csv.gz ]; then
                             echo "found $Name.$Species2.filtered.counts.csv"
			 else
                             echo "count unfiltered reads"
                             zcat $Name.$Species2.bed.gz | \
				 awk -v OFS='\t' '{a[$4] += $5} END{for (i in a) print a[i], i}'| \
				 awk -v OFS='\t' '{if($1 >= '${ReadsPerBarcode[$index]}') print }'> $Name.$Species2.wdup.RG.freq.bed
                             Rscript $myPATH/sum_reads_v2.R $dir/fastqs/ $Name.$Species2.wdup.RG.freq.bed --save
                             mv $Name.$Species2.wdup.RG.freq.bed.csv $Name.$Species2.unfiltered.counts.csv

                             echo "count filtered reads"
                             zcat $Name.$Species2.bed.gz | cut -f4 | uniq -c | \
				 awk -v OFS='\t' '{if($1 >= '${ReadsPerBarcode[$index]}') print }' > $Name.$Species2.rmdup.RG.freq.bed
                             Rscript $myPATH/sum_reads_v2.R $dir/fastqs/ $Name.$Species2.rmdup.RG.freq.bed --save
                             mv $Name.$Species2.rmdup.RG.freq.bed.csv $Name.$Species2.filtered.counts.csv
			     rm $Name.$Species2.wdup.RG.freq.bed $Name.$Species2.rmdup.RG.freq.bed
			 fi

			 # remove barcode with low counts from the fragment file for ATAC
			 if [ -f $Name.$Species2.fragments.tsv.gz ]; then
			     echo "Skip removing low counts barcode combination"
			 else
			     echo "Remove low counts barcode combination"
			     sed -e 's/,/\t/g' $Name.$Species2.filtered.counts.csv | \
				 awk -v OFS=',' 'NR>=2 {if($6 >= '${ReadsPerBarcode[$index]}') print $1,$2,$3,$4,$5} '  > $Name.$Species2.barcodes.txt
			     grep -wFf $Name.$Species2.barcodes.txt  <(zcat $Name.$Species2.bed.gz) | \
				 sort --parallel=$cores -S 40G -k1,1 -k2,2n -k3,3n -k4,4 | bgzip  > $Name.$Species2.fragments.tsv.gz
			     tabix -p bed $Name.$Species2.fragments.tsv.gz
			 fi

			 # get insert-sizes
			 if [ -f $Name.$Species2.rmdup.hist_data.log ]; then
                             echo "Skip checking insert size"
			 else
			     echo "Checking insert size"
			     echo "${oneread[$index]}"
                             echo '' > $Name.$Species2.rmdup.hist_data.log
			     java -jar $picardPATH CollectInsertSizeMetrics VALIDATION_STRINGENCY=SILENT I=$Name.$Species2.st.bam O=$Name.$Species2.rmdup.hist_data.log H=$Name.$Species2.rmdup.hist_data.pdf W=1000  2>>$dir/Run.log
			 fi
			 
			 # get final quality stats
			 if [ ! -f $Name.$Species2.stats.log ]; then
			     echo -e "Chromosome\tLength\tProperPairs\tBadPairs:Raw" > $Name.$Species2.stats.log
			     samtools idxstats $Name.$Species2.st.bam >> $Name.$Species2.stats.log
			     echo -e "Chromosome\tLength\tProperPairs\tBadPairs:Filtered" >> $Name.$Species2.stats.log
			     samtools idxstats $Name.$Species2.rmdup.bam >> $Name.$Species2.stats.log
			 fi
			 
			 # make TSS pileup fig
			 if [ ! -f $Name.$Species2.RefSeqTSS ]; then
                             echo "Create TSS pileup"; set +e
                             python $myPATH/pyMakeVplot_py3.py -a $Name.$Species2.rmdup.bam -b $tssFilesPATH/$Species2.TSS.bed -e 2000 -p ends -v -u -o $Name.$Species2.RefSeqTSS
			 fi
		    done
		fi
			
		#  get dipc pairs
		if [ ${Type[$index]} == "DipC" ] ; then
                    for Species2 in ${Genome2[@]}; do
			cd $dir/fastqs/
			if [ -f $Name.$Species2.pairs.gz ]; then
			    echo "Found $Name.$Species2.pairs.gz, skip parsing pairs"
			else
			    echo "Parse reads $Name.$Species2.pairs.gz"
			    samtools view -@ 8 -h $Name.$Species2.bam \
				| pairtools parse -c $genomeBed/$Species2.chrom.sizes --add-columns mapq --drop-seq --drop-sam \
				| pairtools sort --tmpdir ./ --memory 40G --nproc $cores -o $Name.$Species2.pairs.nonwalk.gz
			    pairtools stats --yaml -o $Name.$Species2.pairs.nonwalk.yaml $Name.$Species2.pairs.nonwalk.gz
			    samtools view -@ 8 -h $Name.$Species2.bam \
				| pairtools parse -c $genomeBed/$Species2.chrom.sizes --add-columns mapq --drop-seq --drop-sam --walks-policy all --min-mapq=1 --no-flip -o $Name.$Species2.pairs.gz
			    ## --no-flip
			fi
			if [ -f $Name.$Species2.pairs.rescue.gz ]; then
			    echo "Found $Name.$Species2.rescue.pairs.gz"
			else
			    echo "Rescue multi-ligated reads"
#			    function simplify {
#				awk -v OFS='\t' '{
#                            	NF = 11		
#                                $9= "."; $10="."	
#                                p = match($1,"_"); $11 = substr($1,p+1,34) ; print}'
#                            }
			    function simplify {
                                awk -v OFS='\t' '{
                                NF = 11
                             	p = match($1,"_"); $11 = substr($1,p+1,34) ; print}'
                            }
#                            function flip {
#				awk -v OFS='\t' -F'[\t]' '
#                                BEGIN { chra2i["!"] = 0 } ARGIND==1 { chra2i[$1] = FNR ; next }
#                                    { if ( chra2i[$2] > chra2i[$4] || ( $2 == $4 && $3 > $5 ) ) {
#                                          print $1 OFS $4 OFS $5 OFS $2 OFS $3 OFS $7 OFS $6 OFS substr($8,2,1) substr($8,1,1) OFS $10 OFS $9 OFS $11 OFS "F"
#                                    } else { print $0 OFS "O" } }
#                          	    ' $chrom_size_file -
#                            }
			    function flip {
                                awk -v OFS='\t' '
                                   { if ( $2 > $4 || ( $2 == $4 && $3 > $5 ) ) {
				         print $1 OFS $4 OFS $5 OFS $2 OFS $3 OFS $7 OFS $6 OFS substr($8,2,1) substr($8,1,1) OFS $10 OFS $9
                                    } else { print $1 OFS $2 OFS $3 OFS $4 OFS $5 OFS $6 OFS $7 OFS $8 OFS $9 OFS $10 }}'
			       }
			    cat <(zcat $Name.$Species2.pairs.gz | head -3000 | grep '^#') <(zcat $Name.$Species2.pairs.gz | grep -v '^#' | simplify | python $myPATH/rescue_walk.py | flip ) | bgzip -@ 8 > $Name.$Species2.pairs.rescue.temp.gz
			    pairtools sort --nproc $cores --memory 40G --tmpdir ./ -o $Name.$Species2.pairs.rescue.gz $Name.$Species2.pairs.rescue.temp.gz
			    rm $Name.$Species2.pairs.rescue.temp.gz
			    pairtools stats --yaml -o $Name.$Species2.pairs.yaml $Name.$Species2.pairs.rescue.gz
			fi
		    
			# remove duplicates
			if [ -f $dir/fastqs/$Name.$Species2.pairs.unique.gz ]; then
			    echo "Skip mark duplicate"
			else
			    echo "Remove duplicates by cutting sites"
                            zcat $Name.$Species2.pairs.rescue.gz | grep -v "#" \
				| awk '{if(length($2)<6 && length($4)<6) print }' \
				| awk 'match($2, "chr") && match($4, "chr")' \
				| sed 's/_/\t/g' \
				| awk -v OFS='\t' '{print $0, $3"_"$4"_"$5"_"$6"_"$2}' \
				| sort --parallel=$cores -S 40G -k12 \
				| awk -v OFS="\t" 'NR==1{ t1=$12;prev=$0;umisum=0} {if(t1==$12) {umisum+=1} else {print prev, umisum ; t1=$12; umisum=1; prev=$0}} END {print prev, umisum}' \
				| tr -s ' ' | sed 's/ /\t/g' | awk -v OFS="\t" '{print $1, $3, $4, $5, $6, $7, $8, $9, $10, $11, $13, $2}' | gzip > $Name.$Species2.pairs.unique.gz
			fi

			# convert to a bam file for QC
			if [ ! -f $Name.$Species2.rmdup.bam ]; then
			    bedToBam -i <(zcat $Name.$Species2.pairs.gz | grep -v \# | \
					  grep UU | \
					  awk '$2 == $4' | \
					  grep chr | \
					  awk -v OFS="\t" '{print $2, $3, $3 + 150, $1}') -g $genomeBed/$Species2.chrom.sizes | \
				samtools sort -@ $cores -m 2G - > $Name.$Species2.rmdup.bam
			    samtools index -@ $cores $Name.$Species2.rmdup.bam
			fi
                
			# count reads for DipC
			if [ -f $Name.$Species2.filtered.counts.csv ] || [ -f $Name.$Species2.filtered.counts.csv.gz ]; then
                            echo "found $Name.$Species2.filtered.counts.csv"
			else
                            echo "count unfiltered contacts"
			    zcat $Name.$Species2.pairs.unique.gz | \
				awk -v OFS='\t' '{a[$12] += $11} END{for (i in a) print a[i], i}'| \
				awk -v OFS='\t' '{if($1 >= '${ReadsPerBarcode[$index]}') print }'> $Name.$Species2.wdup.RG.freq.bed
                            Rscript $myPATH/sum_reads_v2.R $dir/fastqs/ $Name.$Species2.wdup.RG.freq.bed --save
                            mv $Name.$Species2.wdup.RG.freq.bed.csv $Name.$Species2.unfiltered.counts.csv
			    
			    echo "count filtered contacts"
                            zcat $Name.$Species2.pairs.unique.gz | cut -f12 | sort | uniq -c | \
				awk -v OFS='\t' '{if($1 >= '${ReadsPerBarcode[$index]}') print }' > $Name.$Species2.rmdup.RG.freq.bed
                            Rscript $myPATH/sum_reads_v2.R $dir/fastqs/ $Name.$Species2.rmdup.RG.freq.bed --save
                            mv $Name.$Species2.rmdup.RG.freq.bed.csv $Name.$Species2.filtered.counts.csv
                            rm $Name.$Species2.wdup.RG.freq.bed $Name.$Species2.rmdup.RG.freq.bed			
			fi

			# remove barcode with low counts from the contacts file 
			if [ -f $Name.$Species2.pairs.unique.final.gz ]; then
                            echo "Skip removing low counts barcode combination"
			else
                            echo "Remove low counts barcode combination"
			    sed -e 's/,/\t/g' $Name.$Species2.filtered.counts.csv | \
				awk -v OFS=',' 'NR>=2 {if($6 >= '${ReadsPerBarcode[$index]}') print $1,$2,$3,$4,$5}'  > $Name.$Species2.barcodes.txt
                            grep -wFf $Name.$Species2.barcodes.txt  <(zcat $Name.$Species2.pairs.unique.gz) | \
				sort --parallel=$cores -S 40G -k1,1 -k2,2n -k3,3n -k4,4 | bgzip  > $Name.$Species2.pairs.unique.final.gz
			fi

			# get stats for pairs
			if [ -f $dir/fastqs/$Name.$Species2.pairs.stats.txt ]; then
                            echo "Skip calculate stats"
			else
                            echo "Calculate stats for hiC"
                            echo $Name > $Name.$Species2.pairs.stats.txt
                            temp=`zless $Name.$Species2.pairs.unique.final.gz | wc -l`
                            printf "total_pairs: "$temp"\n" >> $Name.$Species2.pairs.stats.txt

                            temp2=`zless $Name.$Species2.pairs.unique.final.gz | awk '$2 == $4' | awk '$3-$5 > 1000 || $3-$5 < -1000' | wc -l`
                            printf "cis_pairs: "$temp2"\n" >> $Name.$Species2.pairs.stats.txt

                            temp3=`zless $Name.$Species2.pairs.unique.final.gz | awk '$2 != $4' | wc -l`
                            printf "trans_pairs: "$temp3"\n" >> $Name.$Species2.pairs.stats.txt
			    
                            temp4=$(($temp2+$temp3))
                            printf "valid_pairs: "$temp4"\n" >> $Name.$Species2.pairs.stats.txt
			    
                            temp5=`echo $temp4/$temp | bc -l`
                            temp5_rounded=`printf "%.3f" $temp5`
                            printf "valid_perc: "$temp5_rounded"\n" >> $Name.$Species2.pairs.stats.txt
			    
                            temp6=`echo $temp2/$temp3 | bc -l`
                            temp6_rounded=`printf "%.3f" $temp6`
                            printf "cis/trans: "$temp6_rounded"\n" >> $Name.$Species2.pairs.stats.txt
			fi

			# filter paris to generate contact files
			if [ -f $dir/fastqs/$Name.$Species2.contacts.final.gz ]; then
                            echo "found $Name.$Species2.contacts.final.gz, skip convert to .contacts file"
			else
                            echo "convert Valid pairs to contacts and make plots"
                            zless $Name.$Species2.pairs.unique.final.gz | \
				awk -v OFS='\t' '{if (($2 == $4 && ($3-$5 > 1000 || $3-$5 < -1000))) print $2, $3, $4, $5, $12; else if ($2 != $4) print $2, $3, $4, $5, $12}' | gzip > $Name.$Species2.contacts.final.gz
                            Rscript $myPATH/sum_contacts.R $dir/fastqs/ $Name.$Species2.contacts.final.gz $Name.$Species2.contacts.pdf --save
			    Rscript $myPATH/plot_contactmap.R $dir/fastqs/ $Name.$Species2.contacts.final.gz $Species2 1000000 $Name.$Species2 --save
			fi
		    
			# convert to .cool
			if [ -f $dir/fastqs/$Name.$Species2.cool ]; then
			    echo "found $Name.$Species2.cool, skip convert to .cool"
			else
			    echo "convert Valid pairs to cool"
			    (zless $Name.$Species2.pairs.unique.final.gz \
				 | cooler cload pairs -c1 2 -p1 3 -c2 4 -p2 5 $genomeBed/$Species2.chrom.sizes:$BaseResolution - $Name.$Species2.cool) 2>>$Name.$Species2.hic.log
			fi

			# convert to .mcool
			if [ -f $dir/fastqs/$Name.$Species2.mcool ]; then
			    echo "found $Name.$Species2.mcool, skip convert to .mcool"
			else
			    echo "convert .cool file to .mcool file"
			    (cooler zoomify --balance --nproc $cores --resolutions $Resolutions --out $Name.$Species2.mcool $Name.$Species2.cool) 2>>$Name.$Species2.hic.log
			fi
			# convert to .hiC
			if [ -f $dir/fastqs/$Name.$Species2.hic ]; then
                            echo "found $Name.$Species2.hiC, skip convert to .hic"
                        else
                            echo "convert Valid pairs to hic"
                            java -Xmx20g -jar $juicerToolsPATH pre $Name.$Species2.pairs.rescue.gz $Name.$Species2.hic $Species2 1>>$Name.$Species2.hic.log
                        fi
		    done
	        fi
								           
		# extract methylation 
		if [ ${taps[$index]} == "T" ]; then
		    doawk() {
                        chr=$1
                        Name=$2
                        Species2=$3
                        zcat "$Name"."$Species2"."$chr".temp.bed.gz | \
                            awk '$3-$8 > 9' | awk '$9-$3 > 9' | \
                            sed 's/\_R1/\tR1/g' | \
                            cut -f2,3,4,5,8 | \
                            awk -v OFS='\t' '{if ($5 == "C") {print $2,$3,$3+1,$1,"1",$4} else if ($5 == "R") {print $2,$3,$3+1,$1,"0",$4}}' | \
			    gzip > "$Name"."$Species2"."$chr".temp.CpG.bed.gz
                    }
		    export -f doawk
		    if [ ${Type[$index]} == "ATAC" ] || [ ${Type[$index]} == "DipC" ]; then
			for Species2 in ${Genome2[@]}; do
			    if [ -f $Name.$Species2.nonCpG.bed.gz ]; then
				echo "Found $Name.$Species2.nonCpG.bed.gz"
			    else
				echo "Extract DNA methylation info, skip extracting DNA methylation"
				chrs=`samtools view -H $Name.$Species2.st.bam | grep chr | cut -f2 | sed 's/SN://g' | grep -v Y | awk '{if(length($0)<6)print}'`
				if [ ${oneread[$index]} == "R1" ]; then
				    parallel 'samtools view -bS '$Name'.'$Species2'.st.bam {} > '$Name'.'$Species2'.{}.bam' ::: `echo ${chrs[@]}`
				else
				    ## only use read1 for methylation calling
				    parallel 'samtools view -bS -f 0x40 '$Name'.'$Species2'.st.bam {} > '$Name'.'$Species2'.{}.bam' ::: `echo ${chrs[@]}`
				fi

				parallel 'biscuit cinread -t cg  -p QNAME,CHRM,CRPOS,BSSTRAND,CRBASE,CQBASE,CRETENTION,QBEG,QEND $biscuitGenome/'$Species2'.fa '$Name'.'$Species2'.{}.bam | gzip > '$Name'.'$Species2'.{}.temp.bed.gz' ::: `echo ${chrs[@]}` > /dev/null 2>&1
				parallel doawk {} $Name $Species2 ::: `echo ${chrs[@]}`
				cat $Name.$Species2.*.temp.CpG.bed.gz > $Name.$Species2.CpG.bed.gz
				
				parallel 'biscuit cinread -t ch  -p QNAME,CHRM,CRPOS,BSSTRAND,CRBASE,CQBASE,CRETENTION,QBEG,QEND $biscuitGenome/'$Species2'.fa '$Name'.'$Species2'.{}.bam | gzip > '$Name'.'$Species2'.{}.temp.bed.gz' ::: `echo ${chrs[@]}` > /dev/null 2>&1
				parallel doawk {} $Name $Species2 ::: `echo ${chrs[@]}`
				cat $Name.$Species2.*.temp.CpG.bed.gz >	$Name.$Species2.nonCpG.bed.gz
			
				rm $Name.$Species2.chr*.temp.bed.gz $Name.$Species2.chr*.temp.CpG.bed.gz
				rm $Name.$Species2.chr*.bam
				zcat  $Name.$Species2.CpG.bed.gz | \
				    awk -v OFS='\t' '{me+=$5}END{print "CpG", me/NR}' > $Name.$Species2.methy.pct.log
			    fi
			    if [ -f $Name.$Species2.CpG.bedGraph ]; then
                                echo "Found $Name.$Species2.CpG.bedGraph"
                            else
				echo "Generate bedgraph"
				echo "track type=bedGraph description=methylation level" >  $Name.$Species2.CpG.bedGraph
				zless $Name.$Species2.CpG.bed.gz | cut -f1,2,3,5 | sort --parallel=$cores -S 40G -k1,1 -k2n,2 -k4,3 | uniq -c | \
				    awk -v OFS="\t" 'NR==1 { t1=$2;t2=$3;t3=$4;count=$1; sum=$5*$1} NR > 1 {if(t1==$2 && t2==$3) {count=count+$1; sum=sum+$5*$1} else {print t1,t2,t3,sum/count*100,sum,count-sum; t1=$2;t2=$3;t3=$4;count=$1; sum=$5*$1}} END {print t1,t2,t3,sum/count*100,sum,count-sum}' >> $Name.$Species2.CpG.bedGraph
			    fi
			done
		    fi
		fi
		
		# RNA processing
		if [ ${Type[$index]} == "RNA" ] && [ ! -f $Name.mm10.unique.st.bam ]; then
		    cd $dir/fastqs/
		    if [ ${Genomes[$index]} == "both" ]; then
			echo "Split into hg and mm"
			chrs1=`samtools view -H $Name.$Species.unique.st.bam | grep GRCh38 | cut -f2 | sed 's/SN://g' | awk '{if(length($0)<13)print}'`
			chrs2=`samtools view -H $Name.$Species.unique.st.bam | grep mm10 | cut -f2 | sed 's/SN://g' | awk '{if(length($0)<13)print}'`
			samtools view -@ $cores -h $Name.$Species.unique.st.bam `echo ${chrs1[@]}` | sed 's/GRCh38_//g' | samtools view -@ $cores -b -o $Name.hg38.unique.st.bam &
			samtools view -@ $cores -h $Name.$Species.unique.st.bam `echo ${chrs2[@]}` | sed 's/mm10___//g' | samtools view -@ $cores -b -o $Name.mm10.unique.st.bam &
			wait
			samtools index -@ $cores $Name.hg38.unique.st.bam &
			samtools index -@ $cores $Name.mm10.unique.st.bam &
			wait
		    else
			echo "Single species is aligned"
		    fi		
		fi

		# cellhashing processing
		if [ ${Type[$index]} == "cellhash" ] && [ ! -f $Name.$Species.unique.st.bam ]; then
		    cd $dir/fastqs/
                    echo "Split into hg and mm"
                    chrs1=`samtools view -H $Name.$Species.unique.st.bam | grep hg38 | cut -f2 | sed 's/SN://g'`
                    chrs2=`samtools view -H $Name.$Species.unique.st.bam | grep mm10 | cut -f2 | sed 's/SN://g'`
                    samtools view -@ $cores -b $Name.$Species.unique.st.bam -o temp1.bam `echo ${chrs1[@]}`
                    samtools view -@ $cores -b $Name.$Species.unique.st.bam -o temp2.bam `echo ${chrs2[@]}`
                    mv temp1.bam $Name.hg38.unique.st.bam
                    mv temp2.bam $Name.mm10.unique.st.bam
                    samtools index -@ $cores $Name.hg38.unique.st.bam &
                    samtools index -@ $cores $Name.mm10.unique.st.bam &
		    wait
                    rm temp1.bam temp2.bam
		fi

		# assign feature to reads
		cd  $dir/fastqs/
		if [ ${Type[$index]} == "RNA" ]; then
		    for Species2 in ${Genome2[@]}; do
			if [ -f $dir/fastqs/$Name.$Species2.exon.featureCounts.bam ] || [ -f $dir/fastqs/$Name.$Species2.wdup.bam.bai ]; then
			    echo "Skip exon feasure count"
			else
		    	    # excliude multimapping, uniquely mapped reads only, -Q 30, for real sample, might consider include multi-mapping
			    echo "Feature counting on exons"
			    # count exon
                            if [ ${keepMultiMapping[$index]} == "T" ]; then
				featureCounts -T $cores -Q 0 -M -a $gtf/$Species2.$refgene.gtf -t exon -g $genename \
                                    -o $Name.$Species2.feature.count.txt -R BAM $Name.$Species2.unique.st.bam >>$dir/Run.log
                            else
				featureCounts -T $cores -Q 30 -a $gtf/$Species2.$refgene.gtf -t exon -g $genename \
                                    -o $Name.$Species2.feature.count.txt -R BAM $Name.$Species2.unique.st.bam >>$dir/Run.log
                            fi
			    mv $Name.$Species2.unique.st.bam.featureCounts.bam $Name.$Species2.exon.featureCounts.bam
			fi
			if [ -f $dir/fastqs/$Name.$Species2.wdup.bam.bai ]; then
			    echo "Skip intron and exon feasure count"
			else
			    # count both intron and exon
			    echo "Count feature on both intron and exon"
			    if [ $keepIntron == "T" ]; then
				if [ ${keepMultiMapping[$index]} == "T" ]; then
				    featureCounts -T $cores -Q 0 -M -a $gtf/$Species2.$refgene.gtf -t gene -g $genename \
						  -o $Name.$Species2.feature.count.txt -R BAM $Name.$Species2.exon.featureCounts.bam >>$dir/Run.log
				# then for reads mapped to multiple genes, only keep if all its alignments are within a single gene
				else
				    featureCounts -T $cores -Q 30 -a $gtf/$Species2.$refgene.gtf -t gene -g $genename \
					-o $Name.$Species2.feature.count.txt -R BAM $Name.$Species2.exon.featureCounts.bam >>$dir/Run.log
				fi
				samtools sort -@ $cores -m 2G -o $Name.$Species2.wdup.bam  $Name.$Species2.exon.featureCounts.bam.featureCounts.bam
				rm $Name.$Species2.exon.featureCounts.bam.featureCounts.bam
			    else
				samtools sort -@ $cores -m 2G -o $Name.$Species2.wdup.bam $Name.$Species2.exon.featureCounts.bam
			    fi
			    samtools index -@ $cores $Name.$Species2.wdup.bam
			fi
		    done
		fi
		
		# add crop barcode to reads
                if [ ${Type[$index]} == "crop" ] || [ ${Type[$index]} == "cite" ] || [ ${Type[$index]} == "cellhash" ]; then
                    for Species2 in ${Genome2[@]}; do
                        if [ -f $dir/fastqs/$Name.$Species2.wdup.bam.bai ]; then
                            echo "Skip adding aligned guide barcode to reads"
                        else
			    samtools view -H $Name.$Species2.unique.st.bam > $Name.$Species2.header.sam
			    samtools view $Name.$Species2.unique.st.bam | awk '{print $0"\tXT:Z:"$3}' > $Name.$Species2.temp.sam
			    cat $Name.$Species2.header.sam $Name.$Species2.temp.sam | samtools view -b -@ $cores -o $Name.$Species2.wdup.bam
                            samtools index $Name.$Species2.wdup.bam
                            rm $Name.$Species2.header.sam $Name.$Species2.temp.sam
			fi
		    done
		fi
		
		# group UMIs
		if [ ${Type[$index]} == "RNA" ] || [ ${Type[$index]} == "cellhash" ] || [ ${Type[$index]} == "crop" ] || [ ${Type[$index]} == "cite" ]; then
		    for Species2 in ${Genome2[@]}; do
			if [ -f $dir/fastqs/$Name.$Species2.grouped.bam ] || [ -f $dir/fastqs/$Name.$Species2.groups.tsv.gz ] || [ -f $dir/fastqs/$Name.$Species2.groups.tsv ] || [ -f $dir/fastqs/$Name.$Species2.groups.dT.tsv ] || [ -f $dir/fastqs/$Name.$Species2.groups.N6.tsv ]; then
			    echo "Found $Name.$Species.groups.bam, skip grouping UMIs"
			else
			    echo "Group reads to unique UMIs"
			    if [ $mode == "regular" ]; then
				## this is previously used. Seems to get more slant and fewer UMIs, but get accurate lib size estimation
				umi_tools group --extract-umi-method=read_id \
                                          --per-gene --gene-tag=XT --per-cell \
                                          -I $dir/fastqs/$Name.$Species2.wdup.bam \
                                          --output-bam -S $dir/fastqs/$Name.$Species2.grouped.bam \
                                          --group-out=$dir/fastqs/$Name.$Species2.groups.tsv --skip-tags-regex=Unassigned >>$dir/Run.log
				
			    else
				## own UMI dedup by matching bc-umi-align position
				if [ ${N6[$index]} == "F" ]; then
                                    samtools view -@ $cores $Name.$Species2.wdup.bam | grep XT:Z: | \
					sed 's/Unassigned_Ambiguity/discard/g' | \
					sed 's/Unassigned_MappingQuality/discard/g' | \
					awk 'gsub(/[_]/,"\t", $1)' | \
					awk -v OFS='\t' '{if($NF ~/discard/){$NF=$(NF-1)} print $5, $6, $2, $3, $NF}' | \
					sed 's/XT:Z://g' > $Name.$Species2.wdup.bed
                                    # remove dup reads
                                    python3 $myPATH/rm_dup_barcode_UMI_v3.py \
                                            -i $Name.$Species2.wdup.bed -o $Name.$Species2.groups.tsv --m 1
                                    rm $Name.$Species2.wdup.bed
				else
				    ## group poly T reads by UMI and cutting position
				    samtools view -@ $cores $Name.$Species2.wdup.bam | \
					grep PT:i:1 | \
					grep XT:Z: | \
                                        sed 's/Unassigned_Ambiguity/discard/g' | \
                                        sed 's/Unassigned_MappingQuality/discard/g' | \
                                        awk 'gsub(/[_]/,"\t", $1)' | \
                                        awk -v OFS='\t' '{if($NF ~/discard/){$NF=$(NF-1)} print $5, $6, $2, $3, $NF}' | \
                                        sed 's/XT:Z://g' > $Name.$Species2.wdup.dT.bed

				    # remove dup reads
                                    python3 $myPATH/rm_dup_barcode_UMI_v3.py \
                                            -i $Name.$Species2.wdup.dT.bed -o $Name.$Species2.groups.dT.tsv --m 1
                                    rm $Name.$Species2.wdup.dT.bed

				     ## group N6 reads by gene, barcode and cutting position
                                    samtools view -@ $cores $Name.$Species2.wdup.bam | \
					grep PT:i:0 | \
					grep XT:Z: | \
                                        sed 's/Unassigned_Ambiguity/discard/g' | \
                                        sed 's/Unassigned_MappingQuality/discard/g' | \
                                        awk 'gsub(/[_]/,"\t", $1)' | \
                                        awk -v OFS='\t' '{if($NF ~/discard/){$NF=$(NF-1)} print $5"_"$6, $2, $NF}' | \
                                        sed 's/XT:Z://g' > $Name.$Species2.wdup.N6.bed
				    
                                    # remove dup reads
				    cat $Name.$Species2.wdup.N6.bed | \
					sort --parallel=$cores -S 24G -k1,1 -k2,2 -k3,3 | \
					uniq -c | \
					awk -v OFS='\t' '{print $3, $4, $1, $2}' > $Name.$Species2.groups.N6.tsv
				    rm $Name.$Species2.wdup.N6.bed
				fi
			    fi
			fi
			# convert groupped UMI to bed file
			if [ -f $dir/fastqs/$Name.$Species2.bed.gz ]; then
			    echo "Skip removing single-read UMIs and GGGG"
			else
			    if [ $removeSingelReadUMI == "T" ]; then
				echo "Filter UMI that has only 1 read and convert to bed file"
				if [ $mode == "regular" ]; then
				    # in regular mode, duplicate UMI reads are kept in the tsv, so I remove those reads first (awk 'NR==1{ id="N"} {if(id != $11 ) {id = $11; print $2, $6, $10}}')
				    less $dir/fastqs/$Name.$Species2.groups.tsv | \
					awk 'gsub(/[_]/,"\t", $1)' | \
                                       	awk 'FNR>1 {if($10 >1){if($9 != "GGGGGGGGGG"){print}}}' | \
                                        awk 'NR==1{ id="N"} {if(id != $11 ) {id = $11; print $2, $6, $10}}' | \
                                        awk -v OFS="\t" 'NR==1{ t1=$1;t2=$2;readsum=0; umisum=0} {if(t1==$1 && t2==$2) {readsum+=$3; umisum+=1} \
                                            else {print t1, t2, umisum, readsum; t1=$1;t2=$2;umisum=1;readsum=$3}}' | \
                                        pigz --fast -p $cores > $Name.$Species2.bed.gz
				else
				    if [ ${N6[$index]} == "F" ]; then
					less $dir/fastqs/$Name.$Species2.groups.tsv | \
				            awk -v OFS="\t" '{if($3 >1){print}}' | \
					    sort  --parallel=$cores -S 24G -k1,1 -k2,2 | \
					    awk -v OFS="\t" 'NR==1 { t1=$1;t2=$2;readsum=0; umisum=0} {if(t1==$1 && t2==$2) {readsum+=$3; umisum+=1} else {print t1, t2, umisum, readsum; t1=$1;t2=$2;umisum=1;readsum=$3}} END {print t1, t2, umisum, readsum}' | \
					    pigz --fast -p $cores > $Name.$Species2.bed.gz
				    else
                                        cat $dir/fastqs/$Name.$Species2.groups.dT.tsv $dir/fastqs/$Name.$Species2.groups.N6.tsv | \
					    awk -v OFS="\t" '{if($3 >1){print}}' | \
                                            sort  --parallel=$cores -S 24G -k1,1 -k2,2 | \
                                            awk -v OFS="\t" 'NR==1 { t1=$1;t2=$2;readsum=0; umisum=0} {if(t1==$1 && t2==$2) {readsum+=$3; umisum+=1} else {print t1, t2, umisum, readsum; t1=$1;t2=$2;umisum=1;readsum=$3}} END {print t1, t2, umisum, readsum}' | \
                                            pigz --fast -p $cores > $Name.$Species2.bed.gz
				    fi
				fi
			    else
				echo "Convert to bed file"
				if [ $mode == "regular" ]; then
				    ## count reads per gene
				    less $dir/fastqs/$Name.$Species2.groups.tsv | \
					awk 'gsub(/[_]/,"\t", $1)' | \
					awk 'FNR>1 {if($9 != "GGGGGGGGGG"){print}}' | \
					awk 'NR==1 { id="N"} {if(id != $11 ) {id = $11; print $2, $6, $10}}' | \
					awk -v OFS="\t" 'NR==1 { t1=$1;t2=$2;readsum=0; umisum=0} {if(t1==$1 && t2==$2) {readsum+=$3; umisum+=1} else {print t1, t2, umisum, readsum; t1=$1;t2=$2;umisum=1;readsum=$3}}' | \
					pigz --fast -p $cores > $Name.$Species2.bed.gz
					
					## 3rd column is umi, 4th column is read
					## note: difference between these two versions of groups.tsv
					## 1) umitools output keep all the barcode-UMI and don't collapse them
					## 2) my script already collpsed them at alignment position level
				else
				    if [ ${N6[$index]} == "F" ]; then
					less $dir/fastqs/$Name.$Species2.groups.tsv | \
					    sort --parallel=$cores -S 24G -k1,1 -k2,2 | \
					    awk -v OFS="\t" 'NR==1 { t1=$1;t2=$2;readsum=0; umisum=0} {if(t1==$1 && t2==$2) {readsum+=$3; umisum+=1} else {print t1, t2, umisum, readsum; t1=$1;t2=$2;umisum=1;readsum=$3}} END {print t1, t2, umisum, readsum}' | \
					    pigz --fast -p $cores > $Name.$Species2.bed.gz
				    else
					## count both dT and N6 UMIs
					cat $dir/fastqs/$Name.$Species2.groups.dT.tsv $dir/fastqs/$Name.$Species2.groups.N6.tsv | \
					    sort --parallel=$cores -S 24G -k1,1 -k2,2 | \
                                            awk -v OFS="\t" 'NR==1 { t1=$1;t2=$2;readsum=0; umisum=0} {if(t1==$1 && t2==$2) {readsum+=$3; umisum+=1} else {print t1, t2, umisum, readsum; t1=$1;t2=$2;umisum=1;readsum=$3}} END {print t1, t2, umisum, readsum}' | \
                                            pigz --fast -p $cores > $Name.$Species2.bed.gz
					
				    fi
				fi
			    fi
			fi

			# count reads for RNA
			if [ -f $Name.$Species2.filtered.counts.csv ] || [ -f $Name.$Species.filtered.counts.csv.gz ]; then
                            echo "found $Name.$Species.filtered.counts.csv"
			else
                            echo "count unfiltered reads"
			    zcat $Name.$Species2.bed.gz | \
				awk -v OFS='\t' '{a[$1] += $4} END{for (i in a) print a[i], i}' | \
				awk -v OFS='\t' '{if($1 >= '${ReadsPerBarcode[$index]}') print }'> $Name.$Species2.wdup.RG.freq.bed
                            Rscript $myPATH/sum_reads_v2.R $dir/fastqs/ $Name.$Species2.wdup.RG.freq.bed --save
                            mv $Name.$Species2.wdup.RG.freq.bed.csv $Name.$Species2.unfiltered.counts.csv

                            echo "count filtered reads"
                            zcat $Name.$Species2.bed.gz | \
				awk -v OFS='\t' '{a[$1] += $3} END{for (i in a) print a[i], i}' | \
				awk -v OFS='\t' '{if($1 >= '${ReadsPerBarcode[$index]}') print }' > $Name.$Species2.rmdup.RG.freq.bed
                            Rscript $myPATH/sum_reads_v2.R $dir/fastqs/ $Name.$Species2.rmdup.RG.freq.bed --save
                            mv $Name.$Species2.rmdup.RG.freq.bed.csv $Name.$Species2.filtered.counts.csv
			    
                            rm $Name.$Species2.wdup.RG.freq.bed $Name.$Species2.rmdup.RG.freq.bed
			fi		
			
			# remove barcode combination that has less then N reads
			if [ -f $Name.$Species2.cutoff.bed.gz ]; then
			    echo "Skip removing low counts barcode combination"
			else
			    echo "Remove low counts barcode combination"
			    sed -e 's/,/\t/g' $Name.$Species2.filtered.counts.csv | \
				awk -v OFS=',' 'NR>=2 {if($6 >= '${ReadsPerBarcode[$index]}') print $1,$2,$3,$4,$5} '  > $Name.$Species2.barcodes.txt
				grep -wFf $Name.$Species2.barcodes.txt <(zcat $Name.$Species2.bed.gz) | pigz --fast -p $cores > $Name.$Species2.cutoff.bed.gz
			fi
		    done
		fi
		
		# Gene body coverage and reads distribution
		if [ ${Type[$index]} == "RNA" ]; then
		    cd $dir/fastqs/
		    for Species2 in ${Genome2[@]}; do
			if [ -f $dir/fastqs/$Name.$Species2.geneBodyCoverage.curves.pdf ]; then
			    echo "Skip calculate gene coverage"
			else
			    echo "Calculate genebody converage and reads distribution of" $Name.$Species2.wdup.bam
			    geneBody_coverage.py -i $dir/fastqs/$Name.$Species2.wdup.bam -r <(head -n 5000 $genomeBed/$Species2.UCSC_RefSeq.bed) -o $Name.$Species2 2>>$dir/Run.log
			fi

			if [ -f $dir/fastqs/$Name.$Species2.read_distribution.txt ]; then
			    echo "Skip calculate read disbution"
			else
			    echo "Calculate read distribution"
			    if [ $Runtype = QC ]; then
				read_distribution.py -i $dir/fastqs/$Name.$Species2.wdup.bam -r $genomeBed/$Species2.UCSC_RefSeq.bed > $Name.$Species2.read_distribution.txt 2>>$dir/Run.log
			    else
				# only use 1% of reads
				samtools view -s 0.01 -o $dir/fastqs/temp.bam $dir/fastqs/$Name.$Species2.wdup.bam
				read_distribution.py -i $dir/fastqs/temp.bam -r $genomeBed/$Species2.UCSC_RefSeq.bed > $Name.$Species2.read_distribution.txt 2>>$dir/Run.log
				rm $dir/fastqs/temp.bam
			    fi
			    # plot reads disbution
			    tail -n +5 $Name.$Species2.read_distribution.txt | head -n -1 > temp1.txt
			    head -n 3  $Name.$Species2.read_distribution.txt | grep Assigned | sed 's/Total Assigned Tags//g' | sed 's/ //g' > temp2.txt
			    Rscript $myPATH/Read_distribution.R $dir/fastqs/ $Name.$Species2 --save
			    rm temp1.txt temp2.txt
			fi
		    done
		fi
		
                if [ ${Type[$index]} == "DipC" ]; then
		    for Species2 in ${Genome2[@]}; do
			# get final quality stats
			echo -e "Chromosome\tLength\tProperPairs\tBadPairs:Raw" > $Name.$Species2.stats.log
			samtools idxstats $Name.$Species2.st.bam >> $Name.$Species2.stats.log
			echo -e "Chromosome\tLength\tProperPairs\tBadPairs:Filtered" >> $Name.$Species2.stats.log
			samtools idxstats $Name.$Species2.rmdup.bam >> $Name.$Species2.stats.log
			if [ ! -f $Name.$Species2.RefSeqTSS ]; then		
			    # make TSS pileup fig
			    echo "Create TSS pileup"; set +e
			    python $myPATH/pyMakeVplot_py3.py -a $Name.$Species2.rmdup.bam -b $tssFilesPATH/$Species2.TSS.bed -e 2000 -p ends -v -u -o $Name.$Species2.RefSeqTSS
			fi
		    done
		fi
		
		if [ -d "$dir/fastqs/tmp" ]; then rm -r $dir/fastqs/tmp; fi
	     fi						       
	done
	
	for Species2 in ${Genome2[@]}; do
	    if [ ${Type[$index]} == "RNA" ] || [ ${Type[$index]} == "crop" ] || [ ${Type[$index]} == "cite" ] || [ ${Type[$index]} == "cellhash" ]; then
	  	# plot UMI/cell or gene/cell
		echo "Convert bed to UMI count matrix, plot UMI_gene_perCell, generate h5"
		Rscript $myPATH/UMI_gene_perCell_plot_v3.R $dir/fastqs/ $Name.$Species2 --save
	    fi
	done

	# estimate lib size
	if [ -f $Name.counts.csv ] || [ -f $Name.counts.csv.gz ]; then
	    echo "Found $Name.counts.csv, skip calculate lib size"
	else
	    echo "Estimate lib size"
            if [ ${Genomes[$index]} == "both" ]; then
		if [ -f $Name.hg38.unfiltered.counts.csv ] && [ ! -f $Name.hg38.filtered.counts.csv ]; then
		    cp $Name.hg38.unfiltered.counts.csv $Name.hg38.filtered.counts.csv
		    echo "Error: Could locate $Name.hg38.filtered.counts.csv"
		fi
		if [ -f $Name.mm10.unfiltered.counts.csv ] && [ ! -f $Name.mm10.filtered.counts.csv ]; then
		    cp $Name.mm10.unfiltered.counts.csv $Name.mm10.filtered.counts.csv
		    echo "Error: Could locate $Name.mm10.filtered.counts.csv"
		fi
		if [ -f $Name.hg38.unfiltered.counts.csv ] && [ -f $Name.mm10.unfiltered.counts.csv ]; then
                    echo "Calcuating library size for $Name"
		    Rscript $myPATH/lib_size_sc_V5_species_mixing.R ./ $Name ${ReadsPerBarcode[$index]} ${Type[$index]} ${oneread[$index]} --save
		fi
	    else
		if [ -f $Name.${Genomes[$index]}.unfiltered.counts.csv ] && [ ! -f $Name.${Genomes[$index]}.filtered.counts.csv ]; then
		    cp $Name.${Genomes[$index]}.unfiltered.counts.csv $Name.${Genomes[$index]}.filtered.counts.csv
		else
		    Rscript $myPATH/lib_size_sc_V5_single_species.R ./ $Name ${ReadsPerBarcode[$index]} ${Genomes[$index]} ${Type[$index]} ${oneread[$index]} --save
		fi
	    fi
	fi
	if [ ! -d $dir/$Name/ ]; then 
	    mkdir $dir/$Name/ && mv $dir/fastqs/$Name.* $dir/$Name
	fi
     fi
     index=$((index + 1))
done

cd $dir
# get stats for ATAC
echo "get stats for ATAC"
echo "Name" > Names.atac.xls
index=0; count=0
if [ -d $dir/temp/ ]; then rm -r $dir/temp/; fi
mkdir $dir/temp/
for Name in ${Project[@]}; do
    if [ ${Type[$index]} == "ATAC" ] && [ ${taps[$index]} == "F" ]; then
	if [ ${Genomes[$index]} == "both" ]; then
	    Genome2=(hg38 mm10)
	    for Species in ${Genome2[@]}; do
		echo $Name.$Species >> Names.atac.xls
		mkdir $dir/temp/$Name.$Species/
		cp $dir/$Name/$Name.$Species.dups.log $dir/temp/$Name.$Species/
		cp $dir/$Name/$Name.$Species.align.log $dir/temp/$Name.$Species/
		cp $dir/$Name/$Name.$Species.stats.log $dir/temp/$Name.$Species/
		cp $dir/$Name/$Name.$Species.rmdup.hist_data.log $dir/temp/$Name.$Species/
		cp $dir/$Name/$Name.$Species.RefSeqTSS $dir/temp/$Name.$Species/
		let "count=count+1"
            done
        else
            Species=${Genomes[$index]}
	    echo $Name.$Species >> Names.atac.xls
	    mkdir $dir/temp/$Name.$Species/
	    cp $dir/$Name/$Name.$Species.dups.log $dir/temp/$Name.$Species/
	    cp $dir/$Name/$Name.$Species.align.log $dir/temp/$Name.$Species/
	    cp $dir/$Name/$Name.$Species.stats.log $dir/temp/$Name.$Species/
	    cp $dir/$Name/$Name.$Species.rmdup.hist_data.log $dir/temp/$Name.$Species/
	    cp $dir/$Name/$Name.$Species.RefSeqTSS $dir/temp/$Name.$Species/
	    let "count=count+1"
	fi
    fi
    index=$((index + 1))
done


if [ $count -eq 0 ]; then
    rm Names.atac.xls
else
    cd $dir/temp/
    cp $dir/Names.atac.xls $dir/temp/
    python $myPATH/pySinglesGetAlnStats_py3.py -i ./Names.atac.xls -d ./
    mv $dir/temp/Names.merged.xls $dir/Names.atac.merged.xls
fi
rm -r $dir/temp/

## get hic stats
echo "get stats for HiC"
echo "Name" > $dir/Names.hic.xls
index=0; count=0
if [ -d $dir/temp/ ]; then rm -r $dir/temp/; fi
mkdir $dir/temp/
cd  $dir
for Name in ${Project[@]}; do
#    echo ${Type[$index]}
    if [ ${Type[$index]} == "DipC" ] || [ ${Type[$index]} == "DipC2" ]; then
#	echo $Name
        if [ ${Genomes[$index]} == "both" ]; then
            Genome2=(hg38 mm10)
            for Species in ${Genome2[@]}; do
                echo $Name.$Species >> Names.hic.xls
                mkdir $dir/temp/$Name.$Species/
                cp $dir/$Name/$Name.$Species.dups.log $dir/temp/$Name.$Species/
#                cp $dir/$Name/$Name.$Species.align.log $dir/temp/$Name.$Species/
                cp $dir/$Name/$Name.$Species.stats.log $dir/temp/$Name.$Species/
                #cp $dir/$Name/$Name.$Species.rmdup.hist_data.log $dir/temp/$Name.$Species/
                cp $dir/$Name/$Name.$Species.RefSeqTSS $dir/temp/$Name.$Species/
                cp $dir/$Name/$Name.$Species.pairs.stats.txt $dir/temp/$Name.$Species/
                let "count=count+1"
            done
        else
            Species=${Genomes[$index]}
            echo $Name.$Species >> Names.hic.xls
            mkdir $dir/temp/$Name.$Species/
            cp $dir/$Name/$Name.$Species.dups.log $dir/temp/$Name.$Species/
            cp $dir/$Name/$Name.$Species.align.log $dir/temp/$Name.$Species/
            cp $dir/$Name/$Name.$Species.stats.log $dir/temp/$Name.$Species/
            #cp $dir/$Name/$Name.$Species.rmdup.hist_data.log $dir/temp/$Name.$Species/
            cp $dir/$Name/$Name.$Species.RefSeqTSS $dir/temp/$Name.$Species/
            cp $dir/$Name/$Name.$Species.pairs.stats.txt $dir/temp/$Name.$Species/
            let "count=count+1"
        fi
    fi
    index=$((index + 1))
done

if [ $count -eq 0 ]; then
    rm Names.hic.xls
else
    cd $dir/temp/
    cp $dir/Names.hic.xls $dir/temp/
    python3 $myPATH/py3HiCAlnStats.py -i $dir/temp/Names.hic.xls -d ./
    mv $dir/temp/Names.merged.xls $dir/Names.hic.merged.xls
#    mv Names.merged.xls.iSize.txt.mat Names.hic.iSize.txt.mat
fi
#rm -r $dir/temp/


# gather useful files
cd $dir
mkdir $dir/Useful
cp $dir/*/*.png $dir/Useful
cp $dir/*/*.pdf $dir/Useful

# generate bigwig
let i=0 
unset Subproject
for (( Index=0; Index<${#Type[@]}; Index++ ))
do
    if [ ${Type[$Index]} == "ATAC" ] || [ ${Type[$Index]} == "DipC" ] || [ ${Type[$Index]} == "DipC2" ]; then
	Subproject[$i]=${Project[$Index]}
	let "i+=1"
    fi
done

count=`ls -1 $dir/Useful/*bw 2>/dev/null | wc -l`
Genome=(hg19 mm10 hg38) # define genome

if [ ! -z $Subproject ]; then
    if [ $count != 0 ]; then
	echo "Found bigwig files, skip converting bam to bigwig"
    else
	echo "Genarate bigwigs"
	for Species in ${Genome[@]}; do
	    parallel --jobs 6 --delay 1 'igvtools count -w 25 -f mean -e 250 '$dir'/{}/{}.'$Species'.rmdup.bam '$dir'/Useful/{}.'$Species'.rmdup.wig '$genomeBed'/'$Species'.chrom.sizes &>>'$dir'/Run.log' ::: `echo ${Subproject[@]}`
	    cd $dir/Useful
	    parallel --jobs 6 --delay 1 'wigToBigWig {} '$genomeBed'/'$Species'genome.bed {.}.bw 2>>'$dir'/Run.log' ::: *wig
	    rm *wig
	done
    fi
fi

# clean up
echo "Cleaning up folder"
cp -f $dir/*/*.csv $dir/Useful/
rm $dir/Useful/*filtered.counts.csv 

pigz --fast -p $cores $dir/*/*.csv
pigz --fast -p $cores $dir/*/*.groups.tsv

rm $dir/*/*wdup.all.bam* $dir/*/*namesort.bam $dir/igv.log $dir/*/*exon.featureCounts.bam 
rm $dir/*/*grouped.bam $dir/*/*unique.st.bam*

if [ $cleanup == "T" ]; then
    rm $dir/*/*.groups.tsv.gz
fi

echo "The pipeline is completed!! Author: Sai Ma <sai.ma2@mssm.edu>"

exit
