#!/bin/bash

# Created by argbash-init v2.10.0
# ARG_OPTIONAL_SINGLE([threads],[t],[Number of threads to use for the analysis (default: 10)],[10])
# ARG_OPTIONAL_BOOLEAN([paired],[p],[Boolean indicating if the sequences are paired (default: on)],[on])
# ARG_OPTIONAL_BOOLEAN([unclassified],[],[Report unclassified reads to unclassified_reads#.fastq (default: on)],[on])
# ARG_OPTIONAL_BOOLEAN([classified],[],[Report classified reads to classified_reads#.fastq (default: on)],[on])
# ARG_OPTIONAL_SINGLE([outpath],[o],[Path to output directory (default: ./)],[./])
# ARG_POSITIONAL_SINGLE([indir],[Directory in which fastq files are stored.])
# ARG_POSITIONAL_SINGLE([manifest],[Path to manifest file.])
# ARG_DEFAULTS_POS()
# ARG_HELP([<Loops over fastq files to assign taxonomy using kraken2 and maxikraken database (2019 version)>])
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
	local first_option all_short_options='tpoh'
	first_option="${1:0:1}"
	test "$all_short_options" = "${all_short_options/$first_option/}" && return 1 || return 0
}

# THE DEFAULTS INITIALIZATION - POSITIONALS
_positionals=()
_arg_indir=
_arg_manifest=
# THE DEFAULTS INITIALIZATION - OPTIONALS
_arg_threads="10"
_arg_paired="on"
_arg_unclassified="on"
_arg_classified="on"
_arg_outpath="./"


print_help()
{
	printf '%s\n' "<Loops over fastq files to assign taxonomy using kraken2 and maxikraken database (2019 version)>"
	printf 'Usage: %s [-t|--threads <arg>] [-p|--(no-)paired] [--(no-)unclassified] [--(no-)classified] [-o|--outpath <arg>] [-h|--help] <indir> <manifest>\n' "$0"
	printf '\t%s\n' "<indir>: Directory in which fastq files are stored."
	printf '\t%s\n' "<manifest>: Path to manifest file."
	printf '\t%s\n' "-t, --threads: Number of threads to use for the analysis (default: 10) (default: '10')"
	printf '\t%s\n' "-p, --paired, --no-paired: Boolean indicating if the sequences are paired (default: on) (on by default)"
	printf '\t%s\n' "--unclassified, --no-unclassified: Report unclassified reads to unclassified_reads#.fastq (default: on) (on by default)"
	printf '\t%s\n' "--classified, --no-classified: Report classified reads to classified_reads#.fastq (default: on) (on by default)"
	printf '\t%s\n' "-o, --outpath: Path to output directory (default: ./) (default: './')"
	printf '\t%s\n' "-h, --help: Prints help"
}


parse_commandline()
{
	_positionals_count=0
	while test $# -gt 0
	do
		_key="$1"
		case "$_key" in
			-t|--threads)
				test $# -lt 2 && die "Missing value for the optional argument '$_key'." 1
				_arg_threads="$2"
				shift
				;;
			--threads=*)
				_arg_threads="${_key##--threads=}"
				;;
			-t*)
				_arg_threads="${_key##-t}"
				;;
			-p|--no-paired|--paired)
				_arg_paired="on"
				test "${1:0:5}" = "--no-" && _arg_paired="off"
				;;
			-p*)
				_arg_paired="on"
				_next="${_key##-p}"
				if test -n "$_next" -a "$_next" != "$_key"
				then
					{ begins_with_short_option "$_next" && shift && set -- "-p" "-${_next}" "$@"; } || die "The short option '$_key' can't be decomposed to ${_key:0:2} and -${_key:2}, because ${_key:0:2} doesn't accept value and '-${_key:2:1}' doesn't correspond to a short option."
				fi
				;;
			--no-unclassified|--unclassified)
				_arg_unclassified="on"
				test "${1:0:5}" = "--no-" && _arg_unclassified="off"
				;;
			--no-classified|--classified)
				_arg_classified="on"
				test "${1:0:5}" = "--no-" && _arg_classified="off"
				;;
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
	local _required_args_string="'indir' and 'manifest'"
	test "${_positionals_count}" -ge 2 || _PRINT_HELP=yes die "FATAL ERROR: Not enough positional arguments - we require exactly 2 (namely: $_required_args_string), but got only ${_positionals_count}." 1
	test "${_positionals_count}" -le 2 || _PRINT_HELP=yes die "FATAL ERROR: There were spurious positional arguments --- we expect exactly 2 (namely: $_required_args_string), but got ${_positionals_count} (the last one was: '${_last_positional}')." 1
}


assign_positional_args()
{
	local _positional_name _shift_for=$1
	_positional_names="_arg_indir _arg_manifest "

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
# Written by Sinaeda Anderssen

### Get arguments
THREADS=$_arg_threads
OUTDIR=$_arg_outpath
UNCL=$_arg_unclassified
CLASS=$_arg_classified
PAIRED=$_arg_paired
INDIR=$_arg_indir
MANIFEST=$_arg_manifest

if [ ! -d ./$OUTDIR ]; then
  mkdir -p ./$OUTDIR;
fi

### Loop through samples and concatenate corresponding ones from different runs
for sample in $(cat $MANIFEST | cut -f1)
do
	if [[ $sample == \#* ]]; then continue; fi
    R1=$INDIR/${sample}_*_R1*.fastq.gz
	R2=$INDIR/${sample}_*_R2*.fastq.gz
	kraken2 --db /space/databases/maxikraken2_1903_140GB/ \
		--paired $R1 $R2 \
		--unclassified-out kraken2_fungi/${sample}_unclassified_reads#.fastq \
		--classified-out kraken2_fungi/${sample}_classified_reads#.fastq \
		--output kraken2_fungi/${sample}_output.txt \
		--report kraken2_fungi/${sample}_report.txt \
		--report-minimizer-data \
		--use-names \
		--gzip-compressed \
		--memory-mapping \
		--threads $THREADS
done

# ^^^  TERMINATE YOUR CODE BEFORE THE BOTTOM ARGBASH MARKER  ^^^

# ] <-- needed because of Argbash
