# sqzlib 
# fastA/Q sequence encoding and compression


## What is it?
At its core sqzlib is small and flexible fastA/Q DNA sequence encoding library that uses either [zlib](https://github.com/madler/zlib) or [zstd](https://github.com/facebook/zstd) for compression. Most importantly, sqzlib is designed to be compatible with [Heng Li's kseq.h](https://github.com/attractivechaos/klib) header only fastA/Q parser. This means that any aplication that uses kseq for reading sequences in the fasta format ( Eg. [minimap2](https://github.com/lh3/minimap2), [seqstats](https://github.com/clwgg/seqstats), [bwa](https://github.com/lh3/bwa) ) can be easily patched to use sqzlib.


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

### Command line options

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

### Link sqzlib 

    //Don't forget zlib and zstd
    gcc yourcode.c libsqz.a -lz -lzstd

### Using kseq.h with libsqz

libsqz tries to emulate zlibs gz API. For example:

    //Open a gzip formatted file with zlib
    gzFile fp = gzopen("filename", "rb");
    
    //Open a sqz formatted file
    sqzFile fp = sqzopen("filename", "rb");
    
This means that in order parse fastA/Q data stored in sqz format, these emulated functions must be passed to the preprocessor macors.
    
    //includes
    
    KSEQ_INIT(sqzFile, sqzread)
    
    int main(...) {
        ...
        sqzFile fp = sqzopen("file", "mode");
        kseq_t *seq = kseq_init(seq);
    
        //Reading with kseq_read()
        ...
        //Close file
        sqzclose(fp);
    }
    
For a complete example, you can look at a [patch](https://github.com/7PintsOfCherryGarcia/sqzlib/tree/master/patches/seqstats) for the [seqstats](https://github.com/clwgg/seqstats) fasta summary statistics program.
