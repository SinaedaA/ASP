### Get all orientations of primers
allOrients <- function(primer) {
    # Create all orientations of the input sequence
    require(Biostrings)
    dna <- DNAString(primer) # The Biostrings works w/ DNAString objects rather than character vectors
    orients <- c(
        Forward = dna, Complement = Biostrings::complement(dna), Reverse = Biostrings::reverse(dna),
        RevComp = Biostrings::reverseComplement(dna)
    )
    return(sapply(orients, toString)) # Convert back to character vector
}

### Count the number of times our primers appear in our fwd and rev reads, considering all primer orientations
primerHits <- function(primer, fn) {
    # Counts number of reads in which the primer is found
    nhits <- vcountPattern(primer, sread(readFastq(fn)), fixed = FALSE)
    return(sum(nhits > 0))
}

#### Wrap dada2 step (without the merge)
dada2_wrap <- function(inpath, filenames = list(), maxEE, truncQ, truncLen, trimLeft, trimRight, minLen) {
    if (length(filenames) == 0) {
        cli_abort("No files given for dada2 analysis")
    }
    else if (length(filenames) == 1) {
        cli_alert_info("Performing dada2 denoising on 'single-end' reads (probably coming from Flash2 merging).")
        ## Define output filenames
        filtered_files <- file.path(inpath, "filtered", basename(filenames[[1]]))
        filtered <- filterAndTrim(filenames[[1]], filtered_files,
            maxN = 0, maxEE = maxEE, truncQ = truncQ, truncLen = truncLen,
            trimLeft = trimLeft, trimRight = trimRight, minLen = minLen,
            rm.phix = TRUE, compress = TRUE, multithread = TRUE
        )
        ## Learning error rates
        errors <- learnErrors(filtered_files, multithread = TRUE)
        pdf(file = paste0(outpath, "/errors_estimated.pdf"))
        plotErrors(errors, nominalQ = TRUE)
        dev.off()
        ## Sample inference
        dada <- dada(filtered_files, err = errors, multithread = TRUE)
        ## Return
        return(list(filtered, filtered_files, dada))
    }
    else if (length(filenames) == 2) {
        cli_alert_info("Performing dada2 denoising on 'paired-end' reads (not yet merged).")
        ## Define output filenames (after real filtering)
        ## this part assumes that all the files will remain after filtering, but some can end up empty, if quality checks were not or badly performed
        ## for example, one of my files from bacteria/run1 has 153 reads, and none pass the filter.
        filtered_F <- file.path(inpath, "filtered", basename(filenames[[1]]))
        filtered_R <- file.path(inpath, "filtered", basename(filenames[[2]]))
        filtered <- filterAndTrim(filenames[[1]], filtered_F, filenames[[2]], filtered_R,
            maxN = 0, maxEE = maxEE, truncQ = truncQ, truncLen = truncLen,
            trimLeft = trimLeft, trimRight = trimRight, minLen = minLen,
            rm.phix = TRUE, compress = TRUE, multithread = TRUE
        )
        print(filtered)
        
        ## Learning error rates
        errors_F <- learnErrors(filtered_F, multithread = TRUE)
        errors_R <- learnErrors(filtered_R, multithread = TRUE)
        pdf(file = paste0(outpath, "/errors_estimated.pdf"))
        plotErrors(errors_F, nominalQ = TRUE)
        plotErrors(errors_R, nominalQ = TRUE)
        dev.off()
        ## Sample inference
        ## so, in this part, dada is looking for some files that don't exist (BactV5V7-BC45-144_S45_L001)
        ## add part to remove filenames from filtered_F and filtered_R when files don't exist.
        ## THIS should be done BEFORE learnErrors !
        dada_F <- dada(filtered_F, err = errors_F, multithread = TRUE)
        dada_R <- dada(filtered_R, err = errors_R, multithread = TRUE)
        ## Return
        return(list(filtered, filtered_F, filtered_R, dada_F, dada_R))
    }
    else { 
        cli_abort("Unsupported number of input files (> 2)")
    }
}

getN <- function(x) sum(getUniques(x))