# Patching a kseq.h using application

Here you will find patches for commonly used progrms that utilize kseq.h as their fastA/Q parser.

Additionally, I have forked these repositories and patched them to use sqzlib. You can clone the patched repositories from [seqstats](https://github.com/7PintsOfCherryGarcia/seqstats), [minimap2](https://github.com/7PintsOfCherryGarcia/minimap2), and [bwa-mem2](https://github.com/7PintsOfCherryGarcia/bwa-mem2).

## Applying a patch

To apply a patch run:

    git apply <patchfile>
    
From within the repository you want to patch.
    
### For example:

seqstats patch:

    git clone https://github.com/clwgg/seqstats
    cd seqstats
    git patch /path/to/sqzlib/patches/seqstats/compat
    make
    ./seqstats <sqzfile>
