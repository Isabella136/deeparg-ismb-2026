#!/bin/bash

eval "$(conda shell.bash hook)"
mamba activate cd-hit

cd-hit -i $db -o db90 -T 0 -M 0 -d 0
cd-hit -i $db -o db80 -T 0 -M 0 -d 0
cd-hit -i $db -o db70 -T 0 -M 0 -d 0
cd-hit -i $db -o db60 -T 0 -M 0 -d 0 -c 0.6 -n 4
cd-hit -i $db -o db50 -T 0 -M 0 -d 0 -c 0.5 -n 3
cd-hit -i $db -o db40 -T 0 -M 0 -d 0 -c 0.4 -n 2