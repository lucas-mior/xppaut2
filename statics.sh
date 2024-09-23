#!/bin/bash

files="src/*.c src/cuda/*.c src/cvode/*.c src/sbml/*.c"
grep -E "^[a-zA-Z0-9_]+ [a-zA-Z0-9_]+\([^)]+\);$" "src/functions.h" \
    | while read sig; do
        name="$(echo "$sig" | sed -E 's/^[a-zA-Z0-9_]+ //; s/\([^)]+\);$//')"
        file="$(grep -l "\<${name}\>" $files)"
        used=$(echo "$files" | wc -l)
        if [ $used -eq 1 ]; then
            echo "${sig}::::${file}"
        elif [ $used -eq 0 ]; then
            echo "${sig}" >> unused_functions.txt
        fi
    done \
| while read work; do 
    sig="$(echo "$work" | awk -F"::::" '{print $1}')"
    file="$(echo "$work" | awk -F"::::" '{print $2}')"

    grep -Fv "$sig" "src/functions.h" > tmp.h
    mv tmp.h "src/functions.h"

    find src -iname "$file" \
        | while read file2; do
        sed -i "1istatic $sig" "${file2}"
    done
done
