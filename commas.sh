#!/bin/bash

find src -iname "*.[ch]" | while read file; do
# int32 avsymfonts[5], avromfonts[5];
#
IDENT="\*?[[:alnum:]_]+"

awk " /^    [[:alnum:]_]+ ($IDENT), ?($IDENT);\$/ {
    print
    type = \$2
    for (i = 2; i <= NF; i += 1) {
        var = gensub(\"(    $IDENT)[,;]/\", \"\\1\", \"g\", \$i);
        printf(\"%s %s;NEWLINELINE\", type, var);
    }
    getline
}{
    print
}" "$file" > "${file}.2"
mv "${file}.2" "$file"

sed -i ':a;N;$!ba; s/NEWLINELINE/\n/g' "$file"
sed -i 's/,;$/;/' "$file"
sed -i 's/;;$/;/' "$file"
done
