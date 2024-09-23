#!/bin/bash

files="src/*.c src/cuda/*.c src/cvode/*.c src/sbml/*.c"
grep -E -l "^[a-zA-Z0-9_]+ (\*?[a-zA-Z0-9_]+, ?)+\*?[a-zA-Z0-9_]+;$" src/auto_x11.c \
    | while read file; do
awk '
/^[a-zA-Z0-9_]+ (*?[a-zA-Z0-9_]+, ?)+*?[a-zA-Z0-9_]+;$/ {
    type = $1
    $1 = "";
    for (i = 2; i <= NF; i += 1) {
        var = gensub("(*?[a-zA-Z0-9_]+)[,;] ?", "\\1", "g", $i);
        printf("%s %s;\n", type, var);
    }
    for (i = 2; i <= NF; i += 1) {
        getline
    }
}
{
    print
}' "$file" | tee "${file}.2"
mv "${file}.2" "$file"
done
