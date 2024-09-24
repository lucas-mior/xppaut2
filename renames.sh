#!/bin/bash

files="src/*.c src/cuda/*.c src/cvode/*.c src/sbml/*.c"
grep -E "^extern [[:alnum:]_]+ [[:alnum:]_]+;$" $files \
    | while read match; do
        type="$(awk '{print $2}' <<< "$match")"
        file="$(awk '{file = gensub("([[:alnum:]_]+):.+", "\\1", "g", $1); print file}' <<< "$match")"
        name="$(awk '{print $3}' <<< "$match")"
        grep -q "^$type $name" "$file" || sed -i "/^extern $type $name/d" "$file"
    done
