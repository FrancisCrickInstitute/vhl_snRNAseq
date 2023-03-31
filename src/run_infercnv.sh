#!/bin/bash
for nih_pid in N045 N059 N088 N090 N23 ; do
  for analysis_mode in subclusters ; do
    echo $nih_pid $analysis_mode
    . src/submit_infercnv.sh
  done
done
