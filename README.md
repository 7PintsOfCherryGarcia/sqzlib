# sqzlib 
# fastA/Q sequence encoding and compression


## What is it?
At its core sqzlib is small and flexible fastA/Q DNA sequence encoding library that uses either [zlib](https://github.com/madler/zlib) or [zstd](https://github.com/facebook/zstd) for compression. For encoding, sqzlib uses a combination of bitpacking, runlength encoding, and quality binning to store fastX DNA data. Afterwards, batches of sequences are compressed into "blocks". These blocks can be independently accessed from each other when reading an sqz encoded file.
Most importantly, sqzlib is designed to be compatible with [Heng Li's kseq.h](https://github.com/attractivechaos/klib) header only fastA/Q parser. This means that any application that uses kseq for reading sequences in the fasta format ( Eg. [minimap2](https://github.com/lh3/minimap2), [seqstats](https://github.com/clwgg/seqstats), [bwa](https://github.com/lh3/bwa) ) can be easily patched to use sqzlib.


## What is contained here?

This repo contains the library developed in C and targeting x86 GNU/Linux (no other platform has been tested), as well as a CLI tool for encoding and decoding fasta and fastq DNA sequence files.

If you want to know how well sqzlib compares to gzip or pigz please head to [benchmarks]()


## A super quick getting started

### Requirements

To build sqzlib you will need:

    make
    gcc
    zlib
    zstd

### Get it

    git clone https://github.com/7PintsOfCherryGarcia/sqzlib.git

### Compile it

    cd sqzlib
    make

This will produce the following files:

    sqz      - Command line tool for encoding and compressing of fasta or fastq files
    
    libsqz.a - sqz library
    
    sqzlib.h - Header file that needs to be included if developing with sqzlib

### sqz CLI encoder

    sqz -h
    
This will print how to use sqz

    [sqz]: Usage:

	sqz [options] <file>

		Options:
		-d 		Decompress previously sqz compressed file.
		-o <file>	Write to outputfile <file>.
		-t <n>		Use <n> compression/decompression threads.
		-l <lib>        Use <lib> compression library.
		                        Default: -l zlib
					Options: zlib, zstd
		-h 		This help message.



### Run sqz on a fastA/Q file

    //Compress
    sqz <input>
    
This will create a new file with the sufix .sqz

    <input>.sqz


    //Decompress
    sqz -d <input>.sqz
    
This will print the decoded/decompressed data to standard output


### Using kseq.h with sqzlib

sqzlib provides a very similar interface to zlib's gz API. For example:

    //Open a gzip formatted file with zlib
    gzFile fp = gzopen("filename", "r");
    
    //Open a sqz formatted file with sqzlib
    sqzFile fp = sqzopen("filename", "r");
    
This means that in order parse fastA/Q data stored in sqz format, these functions must be passed to the preprocessor macros. Below is a simple C program that loops over an sqz encoded fastX file and outputs the number of sequences.

    //countseqsqz.c
    //Compile with
    // gcc -o countseqsqz countseqsqz.c libsqz.a -lz -lzstd
    #include <stdio.h>
    #include "sqzlib.h"
    #include "kseq.h"
    //Indicate kseq we are providing an sqz file that is read with sqzread
    KSEQ_INIT(sqzFile, sqzread)

    int main(int argc, char *argv[]) {
        if (argc < 2) return -1;
        //Open sqzfile
        sqzFile fp = sqzopen(argv[1], "r");
        kseq_t *seq = kseq_init(fp);
        int l;
        size_t n = 0;
        while ( (l = kseq_read(seq)) >= 0 )
            n++;
        fprintf(stdout, "%lu sequences\n", n);
        kseq_destroy(seq);
        //Close sqz file
        sqzclose(fp);
    }


### Patching kseq.h based applications with sqzlib
Patching applications that use kseq.h for fastX parsing requires very few code changes. Basically, the only thing that needs to change are the file type and read function that are provided to the KSEQ_INIT() macro.

    KSEQ_INIT(sqzFile, sqzread)

Additionally, the functions used to open and close the input file will also need to be changed.

    sqzFile fp = sqzopen(argv[1], "r");
    sqzclose(fp);

For an example of a kseq.h based tool patched to use sqzlib, you can look at a [patch](https://github.com/7PintsOfCherryGarcia/sqzlib/tree/master/patches/seqstats) for the [seqstats](https://github.com/clwgg/seqstats) fasta summary statistics program. Patches for other applications can be found [here](https://github.com/7PintsOfCherryGarcia/sqzlib/tree/master/patches).


### sqzlib low level API
In addition to the kseq.h compatibility functions, a low level API will be available for version 1 in the comming weeks. In summary, this API permits the creation of arbitrary sequence "blocks" as well as the extraction of random sequences. In addition, the low level API permits access to sequence blocks in a multihtreaded fashion as blocks are independent of each other. As mentioned before, low level API is not documented yet but you can see how it is being implemented in [sqz_pthread.c](https://github.com/7PintsOfCherryGarcia/sqzlib/blob/master/src/pthread/sqz_pthread.c).

# Upcomming updates
* Retrieval of random sequences from sqzfile
* Flexible quality score binning
* Conservation of masked sequences
* Addition of other additional compression engines
