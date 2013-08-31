#!/bin/sh

# This script is a modification of the make_hg19.sh script from the bowtie2
# distribution.


# Downloads sequence for the sacCer3 version of S. cerevisiae(baker's yeast) from
# UCSC.

GENOME_NAME=sacCer3
UCSC_BASE=http://hgdownload.cse.ucsc.edu/goldenPath/${GENOME_NAME}/chromosomes
CHRS_TO_INDEX="chrI,chr1 \
chrII,chr2 \
chrIII,chr3 \
chrIV,chr4 \
chrV,chr5 \
chrVI,chr6 \
chrVII,chr7 \
chrVIII,chr8 \
chrIX,chr9 \
chrX,chr10 \
chrXI,chr11 \
chrXII,chr12 \
chrXIII,chr13 \
chrXIV,chr14 \
chrXV,chr15 \
chrXVI,chr16 \
chrM,chrM"

BOWTIE_FOLDER=`find ./bin -iname "*bowtie*" -type d | head -1`
echo ${BOWTIE_FOLDER}

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
		zcat ${F} | sed 's/'${1}/${2}'/' > ${OUT_F}  || (echo "Error unzipping $F" && exit 1)
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

echo "Make the gap file"
python make_sacCer3_gap.py fasta/${GENOME_NAME}

