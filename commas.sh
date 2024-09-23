#!/bin/bash

find src -iname "*.[ch]" | while read file; do
# extern double X_LO[10], Y_LO[10], X_HI[10], Y_HI[10];
IDENT="\*?[[:alnum:]_]+"
BRACKETS='\[.*\]'

awk " /^[[:alnum:]_]+ (($IDENT)($BRACKETS)?, )+($IDENT)($BRACKETS)?;\$/ {
# print
    type = \$1
    for (i = 2; i <= NF; i += 1) {
        var = gensub(\"($IDENT)($BRACKETS)?[,;]/\", \"\\1\\2\", \"g\", \$i);
        printf(\"%s %s;NEWLINELINE\", type, var);
    }
    getline
}{
    print
}" "$file" | tee "${file}.2"
mv "${file}.2" "$file"

sed -i ':a;N;$!ba; s/NEWLINELINE/\n/g' "$file"
sed -i 's/,;$/;/' "$file"
sed -i 's/;;$/;/' "$file"
done
