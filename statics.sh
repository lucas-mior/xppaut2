#!/bin/bash

files="src/*.c src/cuda/*.c src/cvode/*.c src/sbml/*.c"
grep -E "^[a-zA-Z0-9_]+ [a-zA-Z0-9_]+\([^)]+\);$" "src/functions.h" \
    | while read sig; do
        name="$(echo "$sig" | sed -E 's/^[a-zA-Z0-9_]+ //; s/\([^)]+\);$//')"

        instances="$(grep "\<${name}\>" $files | wc -l)"

        file="$(grep -l "\<${name}\>" $files)"
        used=$(echo "$file" | wc -l)

        if [ $used -eq 1 ]; then
            echo "${sig}::::${file}"
        fi
        if [ $instances -eq 1 ]; then
            echo "${sig}" >> unused_functions.txt
        fi
    done \
| while read work; do 
    s2="$(echo "$work" | awk -F"::::" '{print $1}')"
    f2="$(echo "$work" | awk -F"::::" '{print $2}')"

    grep -Fv "$s2" "src/functions.h" > tmp.h
    mv tmp.h "src/functions.h"

    sed -i "1istatic $s2" "$f2";
done
