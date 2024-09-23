#!/bin/bash

files="src/*.c src/cuda/*.c src/cvode/*.c src/sbml/*.c"
grep --color=auto -E "^[a-zA-Z0-9_]+ [a-zA-Z0-9_]+\([^)]+\);$" "src/functions.h" \
    | while read -r signature; do
        name="$(echo "$signature" | sed -E 's/^[a-zA-Z0-9_]+ //; s/\([^)]+\);$//')"
        files="$(grep -l "\<${name}\>" $files)"
        used=$(echo "$files" | wc -l)
        [ $used -le 1 ] && echo "${signature}::::${file}"
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
