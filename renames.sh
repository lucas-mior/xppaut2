#!/bin/bash

grep -E "^[[:alnum:]_]+ [[:alnum:]_]+\([^)]+\);$" "src/functions.h" \
    | sed -E 's/^[[:alnum:]_]+ //; s/\(.+\);$//;' \
    | while read function; do
    file="$(grep -E -l "^${function}\(.+\) {$" src/*.c | sed 's/\.c$//; s|src/||')"
    grep -q "$file" <<< "$function" || echo "$function $file"
done
