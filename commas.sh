#!/bin/bash

files="src/*.c src/cuda/*.c src/cvode/*.c src/sbml/*.c"
grep -E -l "^[a-zA-Z0-9_]+ (\*?[a-zA-Z0-9_]+, ?)+\*?[a-zA-Z0-9_]+;$" $files \
    | while read file; do
awk '
/^[a-zA-Z0-9_]+ (*?[a-zA-Z0-9_]+, ?)+*?[a-zA-Z0-9_]+;$/ {
    type = $1
    $1 = "";
    for (i = 2; i <= NF; i += 1) {
        var = gensub("(*?[a-zA-Z0-9_]+)[,;] ?", "\\1", "g", $i);
        printf("%s %s;NEWLINE", type, var);
    }
}
{
    print
}' "$file" > "${file}.2"
mv "${file}.2" "$file"
done
