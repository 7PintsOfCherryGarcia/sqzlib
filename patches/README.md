# Patching a kseq.h using application

Here you will find patches for commonly used progrms that utilize kseq.h as their fastA/Q parser.

## Applying a patch

To apply a patch run:

    git apply <patchfile>
    
Withn the repository you want to patch
    
### For example:

seqstats patch:

    git clone https://github.com/clwgg/seqstats
    cd seqstats
    git patch /path/to/sqzlib/patches/seqstats/compat
    make
    ./seqstats <sqzfile>
