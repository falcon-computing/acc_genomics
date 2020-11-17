#!/bin/bash
data_loc="s3://fcs-data/testbench-data/bwa-sw/"
aws s3 sync --include "*" $data_loc ./
