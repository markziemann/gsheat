#!/bin/bash

DAT=cybrid.xls.tabular.fmt
HUM=$DAT.hum
GMT=c2.cp.reactome.v6.1.symbols.gmt
OUT3COL=3col.txt
ORTHO=mouse2human.txt.sort

#human gene name conversion
head -1 $DAT > $HUM
join -1 1 -2 1 \
$ORTHO \
<(sort -k 1b,1 $DAT | sed 's/_/\t/' | cut -f1,3-) \
| tr ' ' '\t' | cut -f3- >> $HUM

MAXCOL=$(sed -n 5p $HUM | wc -w)
for COL in $(seq 2 $MAXCOL) ; do
  NAME=`head -1 $HUM | cut -f$COL`.rnk
  cut -f1,$COL $HUM | sort -k2g | cut -f1 > $NAME
done

while read line ; do
  GS=`echo $line | cut -d ' ' -f1`
  GENES=`echo $line | cut -d ' ' -f3- `

  for RNK in *rnk ; do
    SUM=`echo $GENES | tr ' ' '\n' \
    | grep -n -wFf - $RNK \
    | cut -d ':' -f1 | numsum`

    echo $GS $RNK $SUM
  done
done < $GMT > 3col.txt

Rscript gsheat.R

