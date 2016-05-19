#!/usr/bin/env zsh


for t in *(/) ; do
    echo $t
    ./exec-test.zsh $t
done
