#!/bin/bash

# Created by argbash-init v2.10.0
# ARG_OPTIONAL_SINGLE([outpath],[o],[Complete path to output directory for concatenated runs (default: current ./)],[./])
# ARG_POSITIONAL_DOUBLEDASH([])
# ARG_POSITIONAL_SINGLE([Metadata],[Metadata file for the samples])
# ARG_POSITIONAL_INF([RunDir],[Directories containing sequencing runs to concatenate],[2])
# ARG_DEFAULTS_POS()
# ARG_HELP([<Help message>])
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
	local first_option all_short_options='oh'
	first_option="${1:0:1}"
	test "$all_short_options" = "${all_short_options/$first_option/}" && return 1 || return 0
}

# THE DEFAULTS INITIALIZATION - POSITIONALS
_positionals=()
_arg_metadata=
_arg_rundir=('' '' )
# THE DEFAULTS INITIALIZATION - OPTIONALS
_arg_outpath="./"


print_help()
{
	printf '%s\n' "<Help message>"
	printf 'Usage: %s [-o|--outpath <arg>] [-h|--help] [--] <Metadata> <RunDir-1> <RunDir-2> [<RunDir-3>] ... [<RunDir-n>] ...\n' "$0"
	printf '\t%s\n' "<Metadata>: Metadata file for the samples"
	printf '\t%s\n' "<RunDir>: Directories containing sequencing runs to concatenate"
	printf '\t%s\n' "-o, --outpath: Complete path to output directory for concatenated runs (default: current ./) (default: './')"
	printf '\t%s\n' "-h, --help: Prints help"
}


parse_commandline()
{
	_positionals_count=0
	while test $# -gt 0
	do
		_key="$1"
		if test "$_key" = '--'
		then
			shift
			test $# -gt 0 || break
			_positionals+=("$@")
			_positionals_count=$((_positionals_count + $#))
			shift $(($# - 1))
			_last_positional="$1"
			break
		fi
		case "$_key" in
			-o|--outpath)
				test $# -lt 2 && die "Missing value for the optional argument '$_key'." 1
				_arg_outpath="$2"
				shift
				;;
			--outpath=*)
				_arg_outpath="${_key##--outpath=}"
				;;
			-o*)
				_arg_outpath="${_key##-o}"
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
	local _required_args_string="'Metadata' and 'RunDir' (2 times)"
	test "${_positionals_count}" -ge 3 || _PRINT_HELP=yes die "FATAL ERROR: Not enough positional arguments - we require at least 3 (namely: $_required_args_string), but got only ${_positionals_count}." 1
}


assign_positional_args()
{
	local _positional_name _shift_for=$1
	_positional_names="_arg_metadata _arg_rundir[0] _arg_rundir[1] "
	_our_args=$((${#_positionals[@]} - 3))
	for ((ii = 0; ii < _our_args; ii++))
	do
		_positional_names="$_positional_names _arg_rundir[$((ii + 2))]"
	done

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
# Script that concatenates sequencing runs when there are several.
# Written by Stefania Concetta Quattro
# Modified by Sinaeda Anderssen

### Get arguments
RUNDIRS=${_arg_rundir[@]}
OUTDIR=$_arg_outpath
SAMPLE_METADATA=$_arg_metadata

### Communicate about your work
echo '#############'
echo "Run directories:"
perl -E 'say join "\n", @ARGV' ${RUNDIRS[@]}
echo 'Output directory of merged libraries is ' $OUTDIR
echo 'Working on samples included in the file ' $SAMPLE_METADATA

if [ ! -d ./$OUTDIR ]; then
  mkdir -p ./$OUTDIR;
fi
MAINDIR=$PWD
cd $OUTDIR

### Loop through samples and concatenate corresponding ones from different runs
for sample in $(cat ../$SAMPLE_METADATA | cut -f1)
do
	if [[ $sample == \#* ]]; then continue; fi
	echo $sample
	mkdir $sample
    cd $sample
	echo "Working directory $PWD"
    cat `printf "$MAINDIR/%s/$sample-*/*_R1_001.fastq.gz " ${RUNDIRS[@]}` > ${sample}_R1_001.fastq.gz
    cat `printf "$MAINDIR/%s/$sample-*/*_R2_001.fastq.gz " ${RUNDIRS[@]}` > ${sample}_R2_001.fastq.gz
	mkdir fastqc
	fastqc --outdir fastqc --extract -t 15 *R1_001.fastq.gz *R2_001.fastq.gz
	cd ../
	echo 'Quality analysis (FASTQC) finished for sample' $sample
done

### Make symlinks to fastqc/ directory and run multiqc on fastqc/
mkdir fastqc/
cd fastqc/
for f in `cat ../sample-metadata.tsv | cut -f1`; do for qcfile in `ls ../1_cat_reads/$f/fastqc/`; do ln -s ../1_cat_reads/$f/fastqc/$qcfile ./$qcfile; done; done
cd -
multiqc fastqc/
# ^^^  TERMINATE YOUR CODE BEFORE THE BOTTOM ARGBASH MARKER  ^^^

# ] <-- needed because of Argbash
