#!/bin/bash

IDENT="[[:alnum:]_]+"

find src -iname "*.[ch]" | while read file; do

awk \
" /^extern $IDENT ($IDENT(\[[^]]+\])?( = \S+)?, )+$IDENT(\[[^]]+\])?( = \S+)?;\$/ {
# print; exit
    static = \$1
    type = \$2
    \$1 = \"\"
    \$2 = \"\"

    split(\$0, array, \",\");
    for (i in array) {
        printf(\"%s %s%s;NEWLINELINE\", static, type, array[i]);
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
