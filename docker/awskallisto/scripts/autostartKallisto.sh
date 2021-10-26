#!/bin/bash
mkdir -p ~/.ncbi
echo '/repository/user/cache-disabled = "true"' > ~/.ncbi/user-settings.mkfg

while true
do
    rm -R /alignment/data/uploads/*
    rm -R /alignment/data/results/*
    python3 scripts/allAlign.py
done