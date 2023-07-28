#!/bin/bash

# m4_ignore(
echo "This is just a script template, not the script (yet) - pass it to 'argbash' to fix this." >&2
exit 11  #)Created by argbash-init v2.10.0
# ARG_POSITIONAL_SINGLE([primer_file], [File containing the primers used for the amplification, forward and reverse separated by a space])
# ARG_POSITIONAL_SINGLE([run_directory], [Directory in which the fastq files from sequencing can be found])
# ARG_OPTIONAL_SINGLE([length], [l], [Minimum length of reads (default: 50)], [50])
# ARG_OPTIONAL_SINGLE([outdir], [o], [Directory for the output of cutadapt], [3_analysis/3.1_cutadapt])
# ARG_DEFAULTS_POS
# ARG_HELP([Script to trim the primers from the sequencing reads.])
# ARGBASH_GO

# [ <-- needed because of Argbash

# vvv  PLACE YOUR CODE HERE  vvv
# For example:
PRIMERS=$_arg_primer_file
RUNDIR=$_arg_run_directory
LENGTH=$_arg_length
OUTDIR=$_arg_outdir
SCRIPT_DIR=$( cd -- "$( dirname -- "${BASH_SOURCE[0]}" )" &> /dev/null && pwd )

### Get primers
FWD=`cat $PRIMERS | cut -d" " -f1`
REV=`cat $PRIMERS | cut -d" " -f2`
### Get there reverse-complements
FWD_RC=`python $SCRIPT_DIR/rev_complement.py $FWD`
REV_RC=`python $SCRIPT_DIR/rev_complement.py $REV`

if [ ! -d $OUTDIR ]; then
  mkdir -p $OUTDIR;
fi
mkdir -p $OUTDIR/untrimmed/
mkdir -p $OUTDIR/tooshort/

## Check program is installed, and if not propose way to install it
command -v cutadapt >/dev/null 2>&1 || { echo -e >&2 "I require cutadapt but it's not installed in your environment. \nTry installing it with conda: 'conda install -c bioconda cutadapt'.  Aborting."; exit 1; }

## Run cutadapt in a for loop
for line in `cat $MANIFEST`; do 
    $SAMPLE=`echo $line | cut -f1`
    $R1=`echo $line | cut -f2`
    $R2=`echo $line | cut -f3`
    cutadapt -g "Fwd_primer=^$FWD;max_error_rate=0.1...Rev_RC=$REV_RC;max_error_rate=0;rightmost" \
            -G "Rev_primer=^$REV;max_error_rate=0.1...Fwd_RC=$FWD_RC;max_error_rate=0;rightmost" \
            --minimum-length $LENGTH \ 
            --too-short-output $OUTDIR/tooshort/${SAMPLE}_R1_tooshort.fastq.gz \
            --too-short-paired-output $OUTDIR/tooshort/${SAMPLE}_R2_tooshort.fastq.gz \
            --untrimmed-output $OUTDIR/untrimmed/${SAMPLE}_R1_untrimmed.fastq.gz \
            --untrimmed-paired-output $OUTDIR/untrimmed/${SAMPLE}_R2_untrimmed.fastq.gz \
            -o $OUTDIR/${SAMPLE}_L001_R1_001.fastq.gz \
            -p $OUTDIR/${SAMPLE}_L001_R2_001.fastq.gz \
            $R1 \ 
            $R2
done

# ^^^  TERMINATE YOUR CODE BEFORE THE BOTTOM ARGBASH MARKER  ^^^

# ] <-- needed because of Argbash
