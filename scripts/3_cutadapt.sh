#!/bin/bash

# Created by argbash-init v2.10.0
# ARG_POSITIONAL_SINGLE([primer_file],[File containing the primers used for the amplification, forward and reverse separated by a space])
# ARG_POSITIONAL_SINGLE([run_directory],[Directory in which the fastq files from sequencing can be found])
# ARG_POSITIONAL_SINGLE([metadata],[Path to metadata file])
# ARG_OPTIONAL_SINGLE([length],[l],[Minimum length of reads (default: 50)],[50])
# ARG_OPTIONAL_SINGLE([outdir],[o],[Directory for the output of cutadapt],[3_analysis/3.1_cutadapt])
# ARG_DEFAULTS_POS()
# ARG_HELP([Script to trim the primers from the sequencing reads.])
# ARGBASH_GO()
# needed because of Argbash --> m4_ignore([
### START OF CODE GENERATED BY Argbash v2.10.0 one line above ###
# Argbash is a bash code generator used to get arguments parsing right.
# Argbash is FREE SOFTWARE, see https://argbash.io for more info


die()
{
	local _ret="${2:-1}"
	test "${_PRINT_HELP:-no}" = yes && print_help >&2
	echo "$1" >&2
	exit "${_ret}"
}


begins_with_short_option()
{
	local first_option all_short_options='loh'
	first_option="${1:0:1}"
	test "$all_short_options" = "${all_short_options/$first_option/}" && return 1 || return 0
}

# THE DEFAULTS INITIALIZATION - POSITIONALS
_positionals=()
_arg_primer_file=
_arg_run_directory=
_arg_metadata=
# THE DEFAULTS INITIALIZATION - OPTIONALS
_arg_length="50"
_arg_outdir="3_analysis/3.1_cutadapt"


print_help()
{
	printf '%s\n' "Script to trim the primers from the sequencing reads."
	printf 'Usage: %s [-l|--length <arg>] [-o|--outdir <arg>] [-h|--help] <primer_file> <run_directory> <metadata>\n' "$0"
	printf '\t%s\n' "<primer_file>: File containing the primers used for the amplification, forward and reverse separated by a space"
	printf '\t%s\n' "<run_directory>: Directory in which the fastq files from sequencing can be found"
	printf '\t%s\n' "<metadata>: Path to metadata file"
	printf '\t%s\n' "-l, --length: Minimum length of reads (default: 50) (default: '50')"
	printf '\t%s\n' "-o, --outdir: Directory for the output of cutadapt (default: '3_analysis/3.1_cutadapt')"
	printf '\t%s\n' "-h, --help: Prints help"
}


parse_commandline()
{
	_positionals_count=0
	while test $# -gt 0
	do
		_key="$1"
		case "$_key" in
			-l|--length)
				test $# -lt 2 && die "Missing value for the optional argument '$_key'." 1
				_arg_length="$2"
				shift
				;;
			--length=*)
				_arg_length="${_key##--length=}"
				;;
			-l*)
				_arg_length="${_key##-l}"
				;;
			-o|--outdir)
				test $# -lt 2 && die "Missing value for the optional argument '$_key'." 1
				_arg_outdir="$2"
				shift
				;;
			--outdir=*)
				_arg_outdir="${_key##--outdir=}"
				;;
			-o*)
				_arg_outdir="${_key##-o}"
				;;
			-h|--help)
				print_help
				exit 0
				;;
			-h*)
				print_help
				exit 0
				;;
			*)
				_last_positional="$1"
				_positionals+=("$_last_positional")
				_positionals_count=$((_positionals_count + 1))
				;;
		esac
		shift
	done
}


handle_passed_args_count()
{
	local _required_args_string="'primer_file', 'run_directory' and 'metadata'"
	test "${_positionals_count}" -ge 3 || _PRINT_HELP=yes die "FATAL ERROR: Not enough positional arguments - we require exactly 3 (namely: $_required_args_string), but got only ${_positionals_count}." 1
	test "${_positionals_count}" -le 3 || _PRINT_HELP=yes die "FATAL ERROR: There were spurious positional arguments --- we expect exactly 3 (namely: $_required_args_string), but got ${_positionals_count} (the last one was: '${_last_positional}')." 1
}


assign_positional_args()
{
	local _positional_name _shift_for=$1
	_positional_names="_arg_primer_file _arg_run_directory _arg_metadata "

	shift "$_shift_for"
	for _positional_name in ${_positional_names}
	do
		test $# -gt 0 || break
		eval "$_positional_name=\${1}" || die "Error during argument parsing, possibly an Argbash bug." 1
		shift
	done
}

parse_commandline "$@"
handle_passed_args_count
assign_positional_args 1 "${_positionals[@]}"

# OTHER STUFF GENERATED BY Argbash

### END OF CODE GENERATED BY Argbash (sortof) ### ])
# [ <-- needed because of Argbash


# vvv  PLACE YOUR CODE HERE  vvv
# For example:
PRIMERS=$_arg_primer_file
RUNDIR=$_arg_run_directory
LENGTH=$_arg_length
OUTDIR=$_arg_outdir
METADATA=$_arg_metadata
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

echo "Forward primer: $FWD"
echo "RC of Forward: $FWD_RC"
echo "Reverse primer: $REV"
echo "RC of Reverse: $REV_RC"

## Check program is installed, and if not propose way to install it
command -v cutadapt >/dev/null 2>&1 || { echo -e >&2 "I require cutadapt but it's not installed in your environment. \nTry installing it with conda: 'conda install -c bioconda cutadapt'.  Aborting."; exit 1; }

## Run cutadapt in a for loop
echo "Starting cutadapt on each sample"
for line in `cat $METADATA`; do
    $SAMPLE=`echo $line | cut -f1`
    $R1=`echo ${SAMPLE}_R1_0001.fastq.gz`
    $R2=`echo ${SAMPLE}_R2_0001.fastq.gz`
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
