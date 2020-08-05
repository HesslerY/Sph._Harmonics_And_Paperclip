#!/bin/bash
var=10
for i in `seq 1 $var` 
do
echo "scale = 2 ; $i/10" | bc
done
