#!/bin/bash
for ii in {0..40}; 
do
    python3 compas_sampler_macleod.py $ii &
done
