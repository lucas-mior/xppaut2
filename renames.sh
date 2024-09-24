#!/bin/bash

grep -E "^[[:alnum:]_]+ \**[[:alnum:]_]+\([^)]+\);$" "src/functions.h" \
    | sed -E 's/^[[:alnum:]_]+ //; s/\(.+\);$//;' \
    | sed -E 's|\*|\\*|g' \
    | while read function; do
    file="$(grep -E -l "^${function}\(.+\) {$" src/*.c | sed 's/\.c$//; s|src/||')"
    grep -E -q "_" <<< "$function" || echo "$function $file"
    # echo "$function $function"
done
