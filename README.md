# sqzlib 
# fastA/Q sequence encoding with zlib compression

## What is it?
At its core sqzlib is small and flexible DNA sequence encoding library that uses zlib "deflate" compression algorithm for storing fasta and fastq files.
Most importantly, sqzlib is designed to be compatible with [Heng Li's kseq.h](https://github.com/attractivechaos/klib) header only fasta parser. This means that any aplication that uses kseq for reading sequences in the fasta format ( Eg. [minimap2](https://github.com/lh3/minimap2), [seqstats](https://github.com/clwgg/seqstats), [bwa](https://github.com/lh3/bwa) ) can be easily patched to use sqzlib.

For more information on how to patch a kseq using aplication with sqzlib, head to the [examples](https://github.com/7PintsOfCherryGarcia/sqzlib/tree/master/examples) directory.

## What is contained here?
This repo contains the library developed in C and targeting any GNU/Linux OS, as well as a CLI tool for encoding and decoding fasta and fastq DNA sequence files.

If you want to know how well sqzlib compares to zlib or VGzip please head to [benchmarks]()

If you want to know how to write code with sqzlib head to [develop]()


## A super quick getting started

### Requirements

To build sqzlib you will need:

    make
    gcc
    zlib

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

    //Don't forget zlib
    gcc yourcode.c libsqz.a -lz
