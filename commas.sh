#!/bin/bash

# int32 auto_ntst = 15, auto_nmx = 200, auto_npr = 50, auto_ncol = 4;
IDENT="[[:alnum:]_]+"
BRACKETS='\[.*\]'

find src -iname "*.[ch]" | while read file; do

awk " /^[[:alnum:]_]+ ($IDENT( = \S+)?, )+$IDENT( = \S+)?;\$/ {
print; exit
    type = \$1
    \$1 = \"\"

    split(\$0, array, \",\");
    for (i in array) {
        printf(\"%s %s;NEWLINELINE\", type, array[i]);
    }
    getline
# }{
#     print
}" "$file" | tee "${file}.2"
# mv "${file}.2" "$file"

# sed -i ':a;N;$!ba; s/NEWLINELINE/\n/g' "$file"
# sed -i 's/,;$/;/' "$file"
# sed -i 's/;;$/;/' "$file"
done
