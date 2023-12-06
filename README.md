# Description of the Amplicon Sequencing analysis workflow.

The experimental procedure can be found in the PDF file (Sample_Prep_Protocol_16S_ITS1-2_AS.pdf) written by E. Yashiro.
For the sequencing, we ended up asking for two runs (in order to get a correct number of total reads) per experiment (2x bacterial samples, 2x fungal samples, 2x oomycete samples). 

Steps in the analysis:
0. Make the directory branching (so the output is uniform, and the scripts follow eachother efficiently)
1. Quality control (`fastqc`, `multifastqc`) + manual inspection
2. Cut adapters from the reads, and remove adapters (`trim_galore`) and primers (`cutadapat`)
3. Denoising the reads (`dada2`)
4. Taxonomic assignment and count table generation
5. Data analysis:
   5.1. Filter low-ASV samples and very low-frequency taxa
   5.2. Alpha-diversity (Shannon and/or Simpson)
   5.3. Beta-diversity
   5.4. 

**NOTE**: The 2 runs are not merged at the beginning of the workflow. _Why?_ Because during denoising, the software (`dada2`) creates a model which represents the sequencing error-rate, and this model (like the error-rate) can change from one run to another. Hence, merging the 2 separate runs can influence the error-model, and as we don't know how problematic this can be, we just run the pipeline on the two runs separately, then add the ASV counts at the end.

Necessary input files:
1. **sample-metadata.tsv** : containing the metadata for the experiment, for each sample
2. **primer_file.txt** : containing the forward and reverse primer sequence, separated by a :space:.
3. **reads** : fastq file containing the sequencing reads

The scripts have been tested, and work on macOS and Linux (Ubuntu). 
Each of them accept `--argument` type arguments (optional) and positional arguments. To know the usage, either follow the workflow, or run `script_name --help`. These help messages were created using argbash **(ref)**.

# Clone the repository to get the code and metadata

Use one of these commands to clone the repository to your computer, or download the repo as a zip file, then unzip it.
```bash
git clone git@github.com:SinaedaA/ASP.git
gh repo clone SinaedaA/ASP
```
The data can be downloaded at [xxx]. 
For each experiment, move that data into the corresponding folder (bacteria, fungi or oomycetes).

```bash
mv path/to/data bacteria/
mv path/to/data fungi/
mv path/to/data oomycetes/
```

# Actual workflow and notes
An example can be found in `1_4_steps.sh`, but here's a detailed step-by-step. This will follow the workflow for bacteria, but specify when things are different for the ITS of fungi and oomycetes.
For every step, a logfile is created, with the TIMESTAMP in the filename. The `2>&1 | tee` specifies that the standard output AND the standard error should be written to the logfile.
## 0. Create the directory branches
Go into the bacteria directory, and make a run1 directory, as well as a logfiles sub-directory, then move into run1. From here, you should never have to change directories again, everything can be run from inside run1. 

```bash
mkdir -p run1/logfiles/
cd run1
```

Now, we can run `0_make_dir_branches.sh`, which is a simple script to create the branching inside run1 (or other base directory of your choice). Of course, replace `path/to/read_dir` with the name of your read directory (for bacteria, the data for run1 is inside AZU03_1-374628257/FASTQ_Generation_2022-12-04_16_22_43Z-636427794).

```bash
mkdir logfiles/
TIMESTAMP=$(date '+%Y%m%d-%H%M%S')
# To get Usage:
../../scripts/0_make_dir_branches.sh --help
../../scripts/0_make_dir_branches.sh path/to/read_dir 2>&1 | tee logfiles/0_make_dir_branches_${TIMESTAMP}.log
```

This will create the following directories: 
- 0_raw_reads: containing symbolic links for all the read files to their location.
- 1_fastqc
- 
- 3_analysis, with subdirectories : 3.1_cutadapt, 3.2_trimming, and 3.3_taxonomy, and 3.4_data_analysis

Check the symbolic links have been created inside 0_raw_reads, by running:
```bash
ls -lrth 0_raw_reads/
```

## 1. Run FastQC and MultiQC on your raw reads
The next step is to check the quality of the sequencing reads, by running `fastqc` and `multiqc`. This is done by the `1_fastqc.sh` script.
This script takes in 2 positional arguments (metadata and read_directory), as well as one optional argument, the output path. If this is not specified, it will use 1_fastqc, and create it if it doesn't exist yet.

```bash
TIMESTAMP=$(date '+%Y%m%d-%H%M%S')
## Usage
../../scripts/1_fastqc.sh --help
../../scripts/1_fastqc.sh ../sample-metadata.tsv 0_raw_reads/ --outpath 1_fastqc/ 2>&1 | tee logfiles/1_fastqc_${TIMESTAMP}.log
```

Then, manually inspect the FastQC output and/or the multiQC output in your browser.

## 2. Optional: Make a manifest
(maybe remove this part, as it was mainly used for QIIME2, which I don't use extensively anymore.)
```bash
TIMESTAMP=$(date '+%Y%m%d-%H%M%S')
../../scripts/2_make_manifest.sh ../sample-metadata.tsv 0_raw_reads/ 2>&1 | tee logfiles/2_make_manifest_${TIMESTAMP}.log
```

## 3. Remove adapters and primers from the sequencing reads
This uses `cutadapt` as well as `trim_galore`, and trims the primers and sequencing adapters (resp.) from each read. Please check both of these programs are installed (althought the script will fail and return possible options for download if it can't find these software). 

- Cutadapt can be installed with conda or micromamba from the bioconda channel (conda install -c bioconda cutadapt). The most recent version as of this writing is 4.4
-  Trim_Galore requires Cutadapt (as it is a Perl wrapper script) and can be downloaded from https://github.com/FelixKrueger/TrimGalore 

```bash
cutadapt --version
~/bin/TrimGalore-0.6.10/trim_galore --version ## this is my path to the executable
```

Because the ITS can be of variable length, it is technically possible that during amplification, the fwd amplicon overlaps with the reverse primers, and so we get the reverse complement of the REV primer on the 3' end. So, as a safety measure, we'll remove the FWD-primer **and** the REV-RC-primer (rev-complement of REV-primer) from the forward read, and the opposite from the reverse read.

It needs 2 positional arguments:
1. *primer file*: can be found in the bacteria, fungi or oomycete directory, and has the format "forward_primer reverse_primer" on a single line, both primers separated by a space.
2. *run directory*: raw_reads directory.
3. *path/to/trim_galore*: path to the Trim_Galore executable file
4. *sample-metadata.tsv*: path to the metadata

And accepts 2 optional arguments:
1. `--length` : minimum length of the reads to be retained in the output (the default is 50 bases)
2. `--outdir` : output directory containing the trimmed reads.

**WHAT DOES IT DO?**
First, it will create subdirectories inside the output directory: `trim_galore` (containing adapter-trimmed sequences), `tooshort` (containing for each sample, reads that were shorter than defined threshold after trimming), `untrimmed` (containing read pairs in which the FWD or the REV reads were not found). *Note*: the reverse-complement of the reads are not a 'mandatory' find, meaning that, if cutadapt doesn't find REV-RC inside the forward reads, but does find FWD at the beginning, the read is considered trimmed anyway. 

Then, it will quickly check if it can find the softwares that are used, and return an error if this is not the case.

Finally, it will remove the sequencing adapters using Trim_Galore (a perl wrapper for cutadapt which already contains commonly used adapters). Right now, the script is made to remove **Nextera** adapters. This can be changed in the script, trim_galore can also accept to detect which adapters are present. Trim_Galore also re-executes FastQC and stores the results in `trim_galore/fastqc`. You can execute `multiqc` on that directory to get a summary report.

On these files, stored in `trim_galore` subdirectory, it will apply cutadapt on these files and remove the primer sequences. 

------ Need to change this in the future -------
**NB:** I added the Trim_Galore step after extensive testing with all cutadapt options to find the optimal ones for our data. In the future, I will check whether all of the options I use are available in Trim_Galore as well, in order to bypass having two *cutting* steps (trim_galore followed by cutadapt). But I'm not sure all options are maintained, which is why it is done in 2 steps. 

```bash
TIMESTAMP=$(date '+%Y%m%d-%H%M%S')
## Usage:
../../scripts/3_cutadapt.sh --help
../../scripts/3_cutadapt.sh ../primer_file.txt 0_raw_reads/ ../sample-metadata.tsv 2>&1 | tee logfiles/3_cutadapt_${TIMESTAMP}.log

## Check how many reads were discarded in each file, and how many passed filter
grep "Pairs discarded as untrimmed:" logfiles/3_cutadapt_${TIMESTAMP}.log
grep "Pairs written (passing filters):" logfiles/3_cutadapt_${TIMESTAMP}.log
```
------Write about trim_galore output (fastqc to check adapter removal)-------

Finally, 3_cutadapt re-executes FastQC on all output files, so we should re-examine the quality of the reads after removing the primers.
We can re-run `multiqc` on this output, for easier interpretation.

```bash
## Check how many reads were discarded in each file, and how many passed filter
grep "Pairs discarded as untrimmed:" logfiles/3_cutadapt_${TIMESTAMP}.log
grep "Pairs written (passing filters):" logfiles/3_cutadapt_${TIMESTAMP}.log
multiqc 3_analysis/3.1_cutadapt/fastqc/ --filename 3_cutadapt_MultiQC.html --outdir 3_analysis/3.1_cutadapt/
```

## 4. Denoise and join
### 4.1. Denoising
Usually, `dada2` is used for denoising reads from amplicon sequencing experiments. *explain what denoising is*

### 4.2. Joining
Once the reads are denoised, they need to be joined/merged, and for this, we can either use `dada2` or another software `flash2`.
*explain difference between the two*

Both of these options are included in the Rscript `4_dada2.R`, and using flash2 can be specified with the `--join` option. 

The script needs 2 positional arguments:
1. *input directory* : absolute path to the output of 3_cutadapt.sh
2. *output directory* : absolute path to where the output should be written

It also accepts a number of optional arguments:
- `--join` : which software to use for joining reads (dada2 or flash2)
- `--bin` : path to local installation of flash2
- `--flashout` : output directory for flash2 output
- `--overlap` : minimum overlap between FOR and REV reads for joining to happen
- `--primerfile` : path to the primer file (same as in the cutadapt step) to check the presence of primers, and the success of cutadapt
As well as a number of trim and trunc options (`Rscript ../../4_dada2.R --help` for more details).

```bash
TIMESTAMP=$(date '+%Y%m%d-%H%M%S')
## Usage:
../../scripts/4_dada2.R --help
## Using dada2 for joining
Rscript ../../scripts/4_dada2.R 3_analysis/3.1_cutadapt/ 3_analysis/3.2_trimming/dada2/ --join dada2 --primerfile ../primer_file.txt --truncQ 2 2>&1 | tee logfiles/4_dada2_${TIMESTAMP}.log
## Using flash2 for joining
Rscript ../../scripts/4_dada2.R 3_analysis/3.1_cutadapt/ 3_analysis/3.2_trimming/dada2_flash2/ --join flash2 --primerfile ../primer_file.txt --truncQ 2 2>&1 | tee logfiles/4_flash2_dada2_${TIMESTAMP}.log
```

Some information on the default parameters in 4_dada2.R (more detailed explanation of each parameter is given when asking for `--help`):
- truncQ = 2 (will truncate when base is encountered with score < 2)
- minLen = 50 (will discard reads shorter than 50 bases)
- maxEE = 0 (allows max 0 error rate)
- trimLeft = 0, trimRight = 0 (no trimming on either side)
- truncLen = 0 (no truncation)
- maxN = 0 (allows 0 Ns, this is *necessary* for dada2 and thus is not an option for the user)

With the default `truncQ=2` we might discard many many reads because of one bad quality base. Therefore, it's usually better to check the quality of the reads after the `cutadapt` step, to be able to specify `trimLeft` and `trimRight`, and/or `truncLen`, based on your actual data.

------Write about output-------

## 5. Taxonomic assignment
For taxonomic assignment, I decided to use IDTAXA (https://microbiomejournal.biomedcentral.com/articles/10.1186/s40168-018-0521-5). This is available in an R package, and the dada2 tutorial described a way to integrate the output from dada2 into IDTAXA.

The R package to be installed in `DECIPHER`, and is described here: http://www2.decipher.codes/Classification.html 

Under Downloads on this website, several pre-modified databases for classification are available, two of which we will need: 	`SILVA SSU r138 (modified)` and 	`UNITE 2021 (unmodified)`. These are RData objects, which I personally stored in my home directory (`~/lib/IDTAXA`). Remember where you store these files, as this will be one of the input parameters to the 5_taxonomy.R script.

Remember to change the `PATHO_TO_SEQTAB` to either dada2, or flash2_dada2 (depending on your choice of joining software in step4).

```bash
TIMESTAMP=$(date '+%Y%m%d-%H%M%S')
## Usage:
Rscript path/to/scripts/5_taxonomy.R --help
## Bacteria:
Rscript ../scripts/5_taxonomy.R 3_analysis/3.2_trimming/PATH_TO_SEQTAB/seqtab_nochim.tsv 3_analysis/3.3_taxonomy/ $HOME/lib/IDTAXA/SILVA_SSU_r138_2019.rdata 2>&1 | tee logfiles/5_taxonomy_${TIMESTAMP}.log
## Fungi
Rscript ../scripts/5_taxonomy.R 3_analysis/3.2_trimming/PATH_TO_SEQTAB/seqtab_nochim.tsv 3_analysis/3.3_taxonomy/ $HOME/lib/IDTAXA/UNITE_v2021_May2021.rdata 2>&1 | tee logfiles/5_taxonomy_${TIMESTAMP}.log
```

------Write about output-------

## 6. Data analysis
Now that we have taxonomy count tables for each run, we can add them up for corresponding samples.
*Reminder*: we kept the runs separate until now, because dada2 estimates run-specific error rates, thus adding everything together in the beginning of the pipeline can bias the results.

Once we have the sample-metadata and the joint count table in the `data_analysis/` directory, we can proceed:

```bash
Rscript ../scripts/6_data_analysis.R new_count_table.tsv sample-metadata.tsv --organism Hv
```