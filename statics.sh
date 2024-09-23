#!/bin/bash

files="src/*.c src/cuda/*.c src/cvode/*.c src/sbml/*.c"
grep -E "^[a-zA-Z0-9_]+ [a-zA-Z0-9_]+\([^)]+\);$" "src/functions.h" \
    | while read sig; do
        name="$(echo "$sig" | sed -E 's/^[a-zA-Z0-9_]+ //; s/\([^)]+\);$//')"
        file="$(grep -l "\<${name}\>" $files)"
        used=$(echo "$file" | wc -l)
        if [ $used -eq 1 ]; then
            echo "${sig}::::${file}"
        elif [ $used -eq 0 ]; then
            echo "${sig}" >> unused_functions.txt
        fi
    done \
| while read work; do 
    sig2="$(echo "$work" | awk -F"::::" '{print $1}')"
    file2="$(echo "$work" | awk -F"::::" '{print $2}')"

    grep -Fv "$sig2" "src/functions.h" > tmp.h
    mv tmp.h "src/functions.h"

    find src -iname "$file2" \
        | while read file3; do
        sed -i "1istatic $sig2" "${file3}"
    done
done
