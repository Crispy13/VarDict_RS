#!/usr/bin/env bash
# Auto-generated pairwise configs from generate_pairwise_configs.py
# 20 configs covering all 2-way parameter interactions
# Parameters: 27

PAIRWISE_CONFIGS=(
  "PW-000-default|-f 0.01 -r 2"
  "PW-001-pileup=on_freq=0.001_minr=1|-p -f 0.001 -r 1 -q 0 -Q 10 -k 0 -U --fisher -X 0 -m 3 -F 0x700 -3 -P 0 -o 0.5 --chimeric --deldupvar -K -x 150 -Y 600 -T 130 -M 25 -u -B 1"
  "PW-002-pileup=on_minr=5_qual=30|-p -f 0.01 -r 5 -q 30 -Q 30 -k 0 --fisher -X 5 -m 15 -F 0x100 -P 10 -o 2.0 --chimeric -K -Y 2000 -w 200 -W 50 -T 130 -UN -A 2 -L 500 -B 4"
  "PW-003-freq=0.05_minr=5_qual=0|-f 0.05 -r 5 -q 0 -Q 30 -U -X 5 -m 3 -F 0x100 -3 -o 2.0 --deldupvar -x 150 -Y 2000 -M 25 -UN -B 4"
  "PW-004-freq=0.001_minr=1_qual=30|-f 0.001 -r 1 -q 30 -U -X 0 -m 15 -3 -P 10 -o 0.5 -x 150 -Y 600 -u -B 1"
  "PW-005-pileup=on_freq=0.001_mapq=30|-p -f 0.001 -r 2 -Q 30 -k 0 -U --fisher -F 0x700 -3 -P 0 --deldupvar -K -T 130 -M 25"
  "PW-006-freq=0.05_qual=30_mapq=10|-f 0.05 -r 2 -q 30 -Q 10 -k 0 -m 3 -F 0x700 -P 0 -o 0.5 --chimeric -x 150 -w 200 -W 50 -M 25 -A 6 -L 2000"
  "PW-007-pileup=on_minr=1_mapq=10|-p -f 0.01 -r 1 -Q 10 -U --fisher -X 5 -o 2.0 --chimeric --deldupvar -K -x 150 -Y 600 -M 25 -u -B 4"
  "PW-008-pileup=on_freq=0.001_minr=5|-p -f 0.001 -r 5 -k 0 -X 0 -F 0x100 -o 0.5 --chimeric --deldupvar -K -Y 2000 -w 200 -W 50 -T 130 -M 25 -u -A 6 -L 2000 -B 1"
  "PW-009-freq=0.05_minr=1_fisher=on|-f 0.05 -r 1 --fisher -m 15 -F 0x700 -3 -P 10 --chimeric --deldupvar -K -Y 600 -w 200 -W 50 -T 130 -UN -A 2 -L 500"
  "PW-010-pileup=on_qual=0_fisher=on|-p -f 0.01 -r 2 -q 0 --fisher -X 5 -m 15 -F 0x700 -3 -P 0 -o 2.0 -w 200 -W 50 -T 130 -A 6 -L 2000 -B 4"
  "PW-011-pileup=on_freq=0.001_qual=0|-p -f 0.001 -r 2 -q 0 -Q 10 -k 0 -U -X 0 -P 10 --deldupvar -x 150 -Y 2000 -T 130 -M 25 -UN -B 4"
  "PW-012-pileup=on_freq=0.001_minr=5|-p -f 0.001 -r 5 -q 0 -Q 10 -k 0 -U -m 15 -F 0x100 -P 0 --deldupvar -M 25 -u"
  "PW-013-pileup=on_freq=0.001_minr=1|-p -f 0.001 -r 1 -q 30 -Q 30 -k 0 -U -X 0 -P 0 -o 2.0 --deldupvar -x 150 -M 25 -UN -B 1"
  "PW-014-pileup=on_qual=0_mapq=30|-p -f 0.01 -r 2 -q 0 -Q 30 -k 0 -U -X 0 -m 3 -P 10 -o 0.5 --deldupvar -Y 600 -B 1"
  "PW-015-pileup=on_freq=0.001_minr=1|-p -f 0.001 -r 1 -q 0 -Q 30 -k 0 -U -X 5 -m 3 -F 0x100 -P 0 --deldupvar -x 150 -Y 2000 -M 25"
  "PW-016-pileup=on_freq=0.001_qual=0|-p -f 0.001 -r 2 -q 0 -Q 30 -k 0 -U -F 0x100 -P 10 -o 2.0 --deldupvar -x 150 -M 25 -u"
  "PW-017-pileup=on_freq=0.001_minr=5|-p -f 0.001 -r 5 -q 0 -Q 10 -k 0 -U -P 0 -o 0.5 --deldupvar -x 150 -Y 600 -B 4"
  "PW-018-pileup=on_freq=0.001_qual=0|-p -f 0.001 -r 2 -q 0 -Q 30 -k 0 -U -X 0 -P 0 -o 0.5 --deldupvar -x 150 -M 25 -UN"
  "PW-019-pileup=on_freq=0.001_qual=0|-p -f 0.001 -r 2 -q 0 -Q 30 -k 0 -U -P 0 --deldupvar -x 150 -Y 2000 -M 25 -B 1"
)
