#!/bin/bash

#this is a comment

runsqz() {
    /usr/bin/time -f %e\\t%M sqz -t 1 $1 2> $1.sqz.time
    /usr/bin/time -f %e sqz -d $1.sqz > dcpsqz 2> $1.sqz.dtime
    ls -l $1.sqz | cut -f5 -d ' ' > $1.sqz.size
    ls -l $1 | cut -f5 -d ' ' > $1.size
    #rm $1.sqz
}

rungzip() {
    /usr/bin/time -f %e\\t%M gzip -9 -c $1 > $1.9.gz 2> $1.9.gzip.time
    /usr/bin/time -f %e gzip -d -c $1.9.gz > dcpgz.9 2> $1.9.gzip.dtime
    ls -l $1.9.gz | cut -f5 -d ' ' > $1.9.gzip.size
    #rm $1.9.gz
    #rm dcpgz.9
    /usr/bin/time -f %e\\t%M gzip -1 -c $1 > $1.1.gz 2> $1.1.gzip.time
    /usr/bin/time -f %e gzip -d -c $1.1.gz > dcpgz.1 2> $1.1.gzip.dtime
    ls -l $1.1.gz | cut -f5 -d ' ' > $1.1.gzip.size
    #rm $1.1.gz
}

runcmp() {
    diff <(seqstats dcpsqz) <(seqstats dcpgz.1)
    cat $1.9.gzip.time $1.sqz.time $1.1.gzip.time > cmptimes
    cat $1.9.gzip.dtime $1.sqz.dtime $1.1.gzip.dtime > dcptimes
    cat $1.9.gzip.size $1.sqz.size $1.1.gzip.size > cmpsizes
    paste <(printf "gzip.9\nsqz\ngzip.1\n") cmptimes dcptimes cmpsizes > tmp
    mv tmp benchmark
    ogsize=$(ls -l $1 | cut -f5 -d ' ')
    cat <(printf "#FMT\tcmp.time\tmax.mem\tdcp.time\t$ogsize\n") benchmark > tmp
    mv tmp benchmark
    #rm *time *.size
    #rm dcpgz.1
    #rm dcpsqz
    #rm cmptimes dcptimes cmpsizes
}

if [ -f $1 ]
then
    runsqz $1
    rungzip $1
    runcmp $1
else
    echo "file $1 not found"
fi



