# sqzlib DNA fastA/Q encoding of sequences + compression

At its core sqzlib is small and flexible DNA sequence encoding library that can use a variety of compression algorithms for storing fasta and fastq files.

This repo contains the C library as well as a CLI tool for encoding and decoding sequence files.

If you want to know how well sqzlib compares to zlib or VGzim please head to [benchmarks](www.example.com)

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
    
This will create a new file with the sufux .sqz

    <input>.sqz

    //Decompress
    sqz -d <input>.sqz
    
This will print the decoded/decompressed data to standard output

### Link sqzlib 

    //Don't forget zlib
    gcc yourcode.c libsqz.a -lz


