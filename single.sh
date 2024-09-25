#!/bin/bash

FUNC_NAME="\**[[:alnum:]_]+"
files="src/*.[ch] src/cvode/*.[ch]"

grep -E "^$FUNC_NAME\(" $files \
| while read match; do
    file="${match%:*}"
    function="${match##*:}"
    function="${function%(*}"

    n="$(grep -F -l "$function" $files | wc -l)"
    if [ $n -eq 2 ]; then
        # m="$(grep -F "$function" "$file" | wc -l)"
        # [ $m -le 2 ] && echo "${function%(*}"
        echo "$function"
    fi
done
