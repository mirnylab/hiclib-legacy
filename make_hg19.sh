#!/bin/sh

# This script is a modification of the make_hg19.sh script from the bowtie2
# distribution. The script can be modified to download the data for the other 
# vertebrate species by the change in GENOME_NAME.

# Downloads sequence for the hg19 version of H. spiens (human) from
# UCSC.

GENOME_NAME=hg19
UCSC_BASE=http://hgdownload.cse.ucsc.edu/goldenPath/${GENOME_NAME}/chromosomes
CHRS_TO_INDEX="chr1,chr1 \
chr2,chr2 \
chr3,chr3 \
chr4,chr4 \
chr5,chr5 \
chr6,chr6 \
chr7,chr7 \
chr8,chr8 \
chr9,chr9 \
chr10,chr10 \
chr11,chr11 \
chr12,chr12 \
chr13,chr13 \
chr14,chr14 \
chr15,chr15 \
chr16,chr16 \
chr17,chr17 \
chr18,chr18 \
chr19,chr19 \
chr20,chr20 \
chr21,chr21 \
chr22,chr22 \
chrX,chrX \
chrY,chrY \
chrM,chrM"

BOWTIE_FOLDER=`find ./bin -iname "*bowtie*" -type d`

get() {
    echo ${1} ${2}
	file=$1
    outfile=$2
	if ! wget --version >/dev/null 2>/dev/null ; then
		if ! curl --version >/dev/null 2>/dev/null ; then
			echo "Please install wget or curl somewhere in your PATH"
			exit 1
		fi
		curl -o `basename $1` $2
		return $?
	else
		wget $1 -O $2
		return $?
	fi
}

BOWTIE_BUILD_EXE=${BOWTIE_FOLDER}/bowtie2-build
if [ ! -x "$BOWTIE_BUILD_EXE" ] ; then
	if ! which bowtie2-build ; then
		echo "Could not find bowtie2-build in current directory or in PATH"
		exit 1
	else
		BOWTIE_BUILD_EXE=`which bowtie2-build`
	fi
fi

OLDIFS=$IFS
INPUTS=
for c in $CHRS_TO_INDEX ; do
    IFS=","
    set $c
    F=${1}.fa.gz
    OUT_F=${2}.fa
	if [ ! -f ${OUT_F} ] ; then
		get ${UCSC_BASE}/$F $F|| (echo "Error getting $F" && exit 1)
		zcat ${F} > ${OUT_F}  || (echo "Error unzipping $F" && exit 1)
        rm ${F}
	fi
	[ -n "$INPUTS" ] && INPUTS=$INPUTS,${OUT_F}
	[ -z "$INPUTS" ] && INPUTS=${OUT_F}
done
IFS=$OLDIFS

CMD="${BOWTIE_BUILD_EXE} ${INPUTS} ${GENOME_NAME}"
echo Running $CMD
if $CMD ; then
	echo "${GENOME_NAME} index built"
else
	echo "Index building failed; see error message"
fi

echo "Move the bowtie index to the bowtie folder"
mkdir -p ${BOWTIE_FOLDER}/index
mv *.bt2 ${BOWTIE_FOLDER}/index

echo "Move the ${GENOME_NAME} FASTA genome to the fasta folder"
mkdir -p fasta/${GENOME_NAME}
mv *.fa fasta/${GENOME_NAME}

echo "Download the gap file"
wget http://hgdownload.cse.ucsc.edu/goldenPath/${GENOME_NAME}/database/gap.txt.gz
zcat gap.txt.gz > fasta/${GENOME_NAME}/gap.txt
rm gap.txt.gz
