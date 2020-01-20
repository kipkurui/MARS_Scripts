#!/usr/bin/env bash
for tf in Egr1 #$(ls /home/kipkurui/Project/Motif_Assessment/GIMME_MAP/)
    do
        cd /home/kipkurui/Project/Motif_Assessment/GIMME_MAP/"$tf"
        echo "$tf"
        gimme roc "$tf"_best.pwm "$tf"_test.fa "$tf"_bg.fa >"$tf"_roc.txt
    done