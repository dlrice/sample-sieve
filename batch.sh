#!/usr/bin/env bash
directory=$1
for file in $( ls $directory ); do
    fullpath=$directory/$file
    echo $fullpath
    ./sample-sieve.sh $fullpath
done
