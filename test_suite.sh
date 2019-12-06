#!/bin/sh
#
# Test script to run all benchmarking tests at the same time
# Disclaimer that some tests may take a long time to run due to
# memory measurement library. Please look at individual benchmark
# python files for instructions on how to comment this out in order
# improve script speed


echo Hash Table Benchmarking results: '\n'
python3 hash_bench_50kb.py
echo '\n'
python3 hash_bench_1mb.py
echo '\n'
python3 hash_bench_100mb.py

# echo '\n\n'
# echo FM Index Benchmarking results: '\n'
# python3 fm_bench_50kb.py
# echo '\n'
# python3 fm_bench_1mb.py
# echo '\n'
# python3 fm_bench_100mb.py


echo '\n\n'
echo Gk Array Benchmarking results: '\n'
python3 gk_bench_50kb.py
echo '\n'
python3 gk_bench_1mb.py
echo '\n'
python3 gk_bench_100mb.py