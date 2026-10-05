#!/bin/bash
# usage: launch.sh U START STATE TAG [extra args]
cd /home/claude/s12push
U=$1; START=$2; STATE=$3; TAG=$4; shift 4
C=/tmp/claude-0/-home-claude/827af079-8944-5ba5-b8f5-20e5625f6a0c/scratchpad/evsp/s12/certificates
rm -rf runs/tmp*
nohup python3 -u runpool5.py --U $U --D 4000000 --N 48000 --rounds 80 --stride 24 --stride0 24 --slack0 0.003 --per_bin0 20 \
  --slack 0.003 --per_bin 20 --pool_top 100 --price_rounds 1 --purge 0.0005 --prune 0.01 --drop_rc 0.02 \
  --pool $C/s12_lower_3.9676.txt,$C/branch/s12_t3.98_corner_k0.txt,$C/branch/s12_t3.98_corner_k1.txt,$C/branch/s12_t3.98_corner_k2.txt \
  --start $START --state $STATE --out runs/$TAG.txt --log runs/$TAG.jsonl "$@" > runs/$TAG.out 2>&1 &
echo $! > runs/$TAG.pid
echo "launched $TAG pid $(cat runs/$TAG.pid)"
