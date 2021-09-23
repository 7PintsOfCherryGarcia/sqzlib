#!/bin/bash

#this is a comment

runsqz() {
    /usr/bin/time -f %e\\t%M sqz -t 1 $1 2> $1.sqz.time
    echo "Done sqz"
    /usr/bin/time -f %e sqz -d -t 1 $1.sqz > dcpsqz 2> $1.sqz.dtime
    echo "Done sqz -d"

    /usr/bin/time -f %e\\t%M sqz -t 2 $1 2> $1.2.sqz.time
    echo "Done sqz 2"
    /usr/bin/time -f %e sqz -d -t 2 $1.sqz > dcpsqz 2> $1.2.sqz.dtime
    echo "Done sqz -d 2"

    /usr/bin/time -f %e\\t%M sqz -t 4 $1 2> $1.4.sqz.time
    echo "Done sqz 4"
    /usr/bin/time -f %e sqz -d -t 4 $1.sqz > dcpsqz 2> $1.4.sqz.dtime
    echo "Done sqz -d 4"

    /usr/bin/time -f %e\\t%M sqz -t 8 $1 2> $1.8.sqz.time
    echo "Done sqz 8"
    /usr/bin/time -f %e sqz -d -t 8 $1.sqz > dcpsqz 2> $1.8.sqz.dtime
    echo "Done sqz -d 8"

    /usr/bin/time -f %e\\t%M sqz -t 16 $1 2> $1.16.sqz.time
    echo "Done sqz 16"
    /usr/bin/time -f %e sqz -d -t 16 $1.sqz > dcpsqz 2> $1.16.sqz.dtime
    echo "Done sqz -d 16"

    ls -l $1.sqz | cut -f5 -d ' ' > $1.sqz.size
    ls -l $1 | cut -f5 -d ' ' > $1.size
    rm $1.sqz
}

rungzip() {
    /usr/bin/time -f %e\\t%M gzip -9 -c dcpsqz > $1.9.gz 2> $1.9.gzip.time
    echo "Done gzip 9"
    ls -l $1.9.gz | cut -f5 -d ' ' > $1.9.gzip.size
    /usr/bin/time -f %e gzip -d -c $1.9.gz > gz.9 2> $1.9.gzip.dtime
    echo "Done gzip -d 9"
    rm $1.9.gz
    rm gz.9

    /usr/bin/time -f %e\\t%M gzip -c dcpsqz > $1.d.gz 2> $1.d.gzip.time
    echo "Done gzip"
    ls -l $1.d.gz | cut -f5 -d ' ' > $1.d.gzip.size
    /usr/bin/time -f %e gzip -d -c $1.d.gz > gz.1 2> $1.d.gzip.dtime
    echo "Done gzip -d"
    rm $1.d.gz

    /usr/bin/time -f %e\\t%M gzip -1 -c dcpsqz > $1.1.gz 2> $1.1.gzip.time
    echo "Done gzip 1"
    ls -l $1.1.gz | cut -f5 -d ' ' > $1.1.gzip.size
    /usr/bin/time -f %e gzip -d -c $1.1.gz > gz.1 2> $1.1.gzip.dtime
    echo "Done gzip -d 9"
    rm $1.1.gz

}

runpigz() {
    /usr/bin/time -f %e\\t%M pigz -p 1 -9 -c dcpsqz > $1.pigz 2> $1.9.pigz.time
    ls -l $1.pigz | cut -f5 -d ' ' > $1.9.pigz.size
    /usr/bin/time -f %e pigz -d -c $1.pigz > gz.9 2> $1.9.pigz.dtime
    rm $1.pigz
    rm gz.9
    echo "Done pigz 9"

    /usr/bin/time -f %e\\t%M pigz -9 -p 2 -c dcpsqz > $1.pigz 2> $1.9.2.pigz.time
    rm $1.pigz
    echo "Done pigz 9 2"

    /usr/bin/time -f %e\\t%M pigz -9 -p 4 -c dcpsqz > $1.pigz 2> $1.9.4.pigz.time
    rm $1.pigz
    echo "Done pigz 9 4"

    /usr/bin/time -f %e\\t%M pigz -9 -p 8 -c dcpsqz > $1.pigz 2> $1.9.8.pigz.time
    rm $1.pigz
    echo "Done pigz 9 8"

    /usr/bin/time -f %e\\t%M pigz -9 -p 16 -c dcpsqz > $1.pigz 2> $1.9.16.pigz.time
    rm $1.pigz
    echo "Done pigz 9 16"


    /usr/bin/time -f %e\\t%M pigz -p 1 -c dcpsqz > $1.pigz 2> $1.d.pigz.time
    ls -l $1.pigz | cut -f5 -d ' ' > $1.d.pigz.size
    /usr/bin/time -f %e pigz -d -c $1.pigz > gz.9 2> $1.d.pigz.dtime
    rm $1.pigz
    rm gz.9
    echo "Done pigz d"

    /usr/bin/time -f %e\\t%M pigz -p 2 -c dcpsqz > $1.pigz 2> $1.d.2.pigz.time
    rm $1.pigz
    echo "Done pigz d 2"

    /usr/bin/time -f %e\\t%M pigz -p 4 -c dcpsqz > $1.pigz 2> $1.d.4.pigz.time
    rm $1.pigz
    echo "Done pigz d 4"

    /usr/bin/time -f %e\\t%M pigz -p 8 -c dcpsqz > $1.pigz 2> $1.d.8.pigz.time
    rm $1.pigz
    echo "Done pigz d 8"

    /usr/bin/time -f %e\\t%M pigz -p 16 -c dcpsqz > $1.pigz 2> $1.d.16.pigz.time
    rm $1.pigz
    echo "Done pigz d 16"


    /usr/bin/time -f %e\\t%M pigz -p 1 -1 -c dcpsqz > $1.pigz 2> $1.1.pigz.time
    ls -l $1.pigz | cut -f5 -d ' ' > $1.1.pigz.size
    /usr/bin/time -f %e pigz -d -c $1.pigz > gz.9 2> $1.1.pigz.dtime
    rm $1.pigz
    rm gz.9
    echo "Done pigz 1"

    /usr/bin/time -f %e\\t%M pigz -1 -p 2 -c dcpsqz > $1.pigz 2> $1.1.2.pigz.time
    rm $1.pigz
    echo "Done pigz 1 2"

    /usr/bin/time -f %e\\t%M pigz -1 -p 4 -c dcpsqz > $1.pigz 2> $1.1.4.pigz.time
    rm $1.pigz
    echo "Done pigz 1 4"

    /usr/bin/time -f %e\\t%M pigz -1 -p 8 -c dcpsqz > $1.pigz 2> $1.1.8.pigz.time
    rm $1.pigz
    echo "Done pigz 1 8"

    /usr/bin/time -f %e\\t%M pigz -1 -p 16 -c dcpsqz > $1.pigz 2> $1.1.16.pigz.time
    rm $1.pigz
    echo "Done pigz 1 16"

}

runcmp() {
    diff <(seqstats dcpsqz) <(seqstats gz.1)
    cat $1.d.gzip.time $1.9.gzip.time $1.1.gzip.time $1.d.pigz.time $1.d.2.pigz.time $1.d.4.pigz.time $1.d.8.pigz.time $1.d.16.pigz.time $1.9.pigz.time $1.9.2.pigz.time $1.9.4.pigz.time $1.9.8.pigz.time $1.9.16.pigz.time $1.1.pigz.time $1.1.2.pigz.time $1.1.4.pigz.time $1.1.8.pigz.time $1.1.16.pigz.time $1.sqz.time $1.2.sqz.time $1.4.sqz.time $1.16.sqz.time $1.16.sqz.time > cmptimes
    paste <(printf "gzip\ngzip.9\ngzip.1\npigz\npigz.2\npigz.4\npigz.8\npigz.16\npigz.9\npigz.9.2\npigz.9.4\npigz.9.8\npigz.9.16\npigz.1\npigz.1.2\npigz.1.4\npigz.1.8\npigz.1.16\nsqz\nsqz.2\nsqz.4\nsqz.8\nsqz.16\n") cmptimes > tmp; mv tmp cmptimes
    cat $1.d.gzip.dtime $1.9.gzip.dtime $1.1.gzip.dtime $1.d.pigz.dtime $1.9.pigz.dtime $1.1.pigz.dtime $1.sqz.dtime $1.2.sqz.dtime $1.4.sqz.dtime $1.8.sqz.dtime $1.16.sqz.dtime > dcptimes
    paste <(printf "gzip\ngzip.9\ngzip.1\npigz\npigz.9\npigz.1\nsqz\nsqz.2\nsqz.4\nsqz.8\nsqz.16\n") dcptimes > tmp; mv tmp dcptimes
    cat $1.d.gzip.size $1.9.gzip.size $1.1.gzip.size $1.d.pigz.size $1.9.pigz.size $1.1.pigz.size $1.sqz.size  > cmpsizes
    paste <(printf "gzip\ngzip.9\ngzip.1\npigz\npigz.9\npigz.1\nsqz\n") cmpsizes > tmp; mv tmp cmpsizes
    ogsize=$(ls -l $1 | cut -f5 -d ' ')
    cat <(printf "#$ogsize\n") cmpsizes > tmp; mv tmp cmpsizes
    #cat <(printf "#FMT\tcmp.time\tmax.mem\tdcp.time\t$ogsize\n") benchmark > tmp
    #mv tmp benchmark
    rm *time *.size
    rm gz.1
    rm dcpsqz
    #rm cmptimes dcptimes cmpsizes
}

if [ -f $1 ]
then
    runsqz $1
    rungzip $1
    runpigz $1
    runcmp $1
else
    echo "file $1 not found"
fi



