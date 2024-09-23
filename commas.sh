#!/bin/bash

find src -iname "*.[ch]" | while read file; do
awk '
/^[[:alnum:]_]+ (\*?[[:alnum:]_]+, )+(\*?[[:alnum:]_]+);$/ {
    type = $1
    for (i = 2; i <= NF; i += 1) {
        var = gensub("(\*?[[:alnum:]_]+)[,;]", "\\1", "g", $i);
        printf("%s %s;NEWLINELINE", type, var);
    }
    getline
}
{
    print
}' "$file" > "${file}.2"
mv "${file}.2" "$file"

sed -i ':a;N;$!ba; s/NEWLINELINE/\n/g' "$file"
done
