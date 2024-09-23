#!/bin/bash

find src -iname "*.[ch]" | while read file; do
# int32 avsymfonts[5], avromfonts[5];
#
IDENT="\*?[[:alnum:]_]+"
BRACKETS='\[..?.?.?\]'

awk " /^extern [[:alnum:]_]+ ($IDENT)($BRACKETS)?, ?($IDENT)($BRACKETS)?;\$/ {
    type = \$2
    for (i = 2; i <= NF; i += 1) {
        var = gensub(\"($IDENT)($BRACKETS)?[,;]/\", \"\\1\\2\", \"g\", \$i);
        printf(\"extern %s %s;NEWLINELINE\", type, var);
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
