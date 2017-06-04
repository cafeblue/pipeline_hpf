#!/bin/awk -f
#print gensub(/(.+?)\.[0-9]+\t(.+)/, "\\1\t\\2\t\\2", "g")
{print $1"\t"$2"\t"$2}
