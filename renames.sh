#!/bin/bash

grep -E -R -l "^extern [[:alnum:]_]+ [[:alnum:]_]+;$"
