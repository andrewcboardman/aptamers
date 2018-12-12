#!/usr/bin/bash 

python gen_seq.py -i J_ones -o 10k_b1_1 -N 10000 -b 1
python gen_seq.py -i J_ones -o 10k_b1_2 -N 10000 -b 1
python gen_seq.py -i J_ones -o 10k_b1_3 -N 10000 -b 1

python gen_seq.py -i h_ones -o 10k_b1_1 -N 10000 -b 1
python gen_seq.py -i h_ones -o 10k_b1_2 -N 10000 -b 1
python gen_seq.py -i h_ones -o 10k_b1_3 -N 10000 -b 1

python gen_seq.py -i h_random -o 10k_b1_1 -N 10000 -b 1
python gen_seq.py -i h_random -o 10k_b1_2 -N 10000 -b 1
python gen_seq.py -i h_random -o 10k_b1_3 -N 10000 -b 1
python gen_seq.py -i h_random -o 10k_b0.1 -N 10000 -b 0.1


python gen_seq.py -i J_random -o 10k_b1_1 -N 10000 -b 1
python gen_seq.py -i J_random -o 10k_b1_2 -N 10000 -b 1
python gen_seq.py -i J_random -o 10k_b1_3 -N 10000 -b 1
python gen_seq.py -i J_random -o 10k_b0.1 -N 10000 -b 0.1




