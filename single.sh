#!/bin/bash

FUNC_NAME="[[:alnum:]_]+"
files="src/*.c src/*.h"

grep -E "^$FUNC_NAME\([[:alnum:]_]+ [[:alnum:]_]+\)" $files \
| while read match; do
    file="${match%:*}"
    function="${match##*:}"
    function="${function%(*}"

    n="$(grep -F -l "$function" $files | wc -l)"
    if [ $n -eq 1 ]; then
        m="$(grep -F "$function" "$file" | wc -l)"
        [ $m -le 3 ] && echo "$function"
    fi
done
