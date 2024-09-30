#!/bin/bash

IDENT="[[:alnum:]_]+"
ASSIGN=" = \S+"
BRACKET="\[[^]]*\]"

find src -iname "*.[ch]" | while read file; do

awk \
" /^    +$IDENT ($IDENT($BRACKET)?($ASSIGN)?, )+$IDENT($BRACKET)?($ASSIGN)?;\$/ {
# print; exit
    type = \$1
    \$1 = \"\"

    split(\$0, array, \",\");
    for (i in array) {
        printf(\"%s %s;NEWLINELINE\", type, array[i]);
    }
    getline
}{
    print
}" \
"$file" | tee "${file}.2"
mv "${file}.2" "$file"

sed -i ':a;N;$!ba; s/NEWLINELINE/\n/g' "$file"
sed -i 's/,;$/;/' "$file"
sed -i 's/;;$/;/' "$file"
done
