### Patching a kseq.h using application

Here you will find patches for commonly used progrms that utilize kseq.h as their fastA/Q parser.

## Applying a patch

To apply a patch just run:

    git apply <patchfile>
    
# For example:

seqstats patch:

    git clone https://github.com/clwgg/seqstats
    cd seqstats
    git patch /path/to/sqzlib/patches/seqstats/compat
    make
    ./seqstats <sqzfile>
    
# Compatibility with gzip files

sqzlib will fall back to zlib routines if provided by a file not in the sqz format.
