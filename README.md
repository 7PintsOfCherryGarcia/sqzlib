# sqzlib DNA fastA/Q encoding of sequences + compression

## What is it?
At its core sqzlib is small and flexible DNA sequence encoding library that can use a variety of compression algorithms for handling fasta and fastq files.
Most importantly, sqzlib is designed to be compatible with Heng [Li's kseq.h](www.example.com) header only fasta parser. This means that any aplication that uses
kseq for reading sequences in the fasta format ( Eg. [minimap2](www.example.com), [seqstats](www.example.com), [yak](www.example.com) ) can be easily patched to
use sqzlib.

For more information on how to patch a kseq using aplication with sqzlib, head to the [examples](www.example.com) directory. If you clone this repo recurrsively,
you will find patched versions of minimap2, seqstats, and bwa that will readily take sqzlib encoded fastA/Q files. 

## What is contained here?
This repo contains the library developed in C and targeting any GNU/Linux OS, as well as a CLI tool for encoding and decoding fasta and fastq DNA sequence files.

If you want to know how well sqzlib compares to zlib or VGzip please head to [benchmarks](www.example.com)

If you want to know how to write code with sqzlib head to [develop](www.example.com)

If you want to know more about the internals of sqzlib head to [How it works?](www.example.com)

## A super quick getting started 

### Get it

    git clone --recursive https://github.com/7PintsOfCherryGarcia/sqzlib.git

### Compile it

    cd squeezma
    make

This will produce the following files:

    sqz      - Command line tool for encoding and compressing of fasta or fastq files
    
    libsqz.a - sqz library
    
    sqzlib.h, sqz_data.h - Header files that need to be included if using sqzlib

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


