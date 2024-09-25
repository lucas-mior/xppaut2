#!/bin/bash

rm -- *.2 *.3

TYPE="[[:alnum:]_]+"
FUNC="\**[[:alnum:]_]+"
files="src/*.c"

grep -E "^$FUNC\(" $files \
| while read match; do
    file="${match%:*}"
    function="${match##*:}"

    file="$(echo "$file" | sed -E 's|^src/(cude/\|cvode/)?||g; s/\.c$//;')"
    function="$(echo "$function" | sed -E 's/\(.+//')"
    
    [[ "$file" =~ "autlib" ]] && continue

    grep -q "$file" <<< "$function" || echo "$function" >> "${file}.2"
    grep -q "_" <<< "$function" || echo "$function" >> "${file}.2"
done

for file2 in *.2; do
    sort "$file2" | uniq | tee "${file2/.2/.3}"
done

rm -- *.2
