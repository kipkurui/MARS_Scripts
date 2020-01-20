#!/usr/bin/env bash
db="Data/Motifs/rap_unijasol_chipseq.meme Data/Motifs/unijasjol_chipseq.jasme Data/Motifs/unijasjol_chipseq.jolme Data/Motifs/unijasjol_chipseq.unime Data/Motifs/uniprobe_chipseq.meme Data/Motifs/zhao_unijasol_chipseq.meme"
for tf in $(cat Data/unijasjol_chipseq.txt);
do
	for cl in $(grep -i $tf encode);
	do
		#mkdir -p Results/Chip-seq/$tf
		centrimo -seqlen 500 -verbosity 1 -oc Results/Chip-seq/$tf/$cl/centrimo_out -bgfile Results/Chip-seq/$tf/$cl/background Results/Chip-seq/$tf/$cl/$cl.fa Results/Chip-seq/$tf/$cl/meme_out/meme.xml Data/Motifs/rap_unijasol_chipseq.meme Data/Motifs/jaspar_unijasjol_chipseq.meme Data/Motifs/Jolma_unijasjol_chipseq.meme Data/Motifs/Uniprobe_unijasjol_chipseq.meme Data/Motifs/zhao_unijasol_chipseq.meme Data/Motifs/hucomoco_chipseq.meme Data/Motifs/Kheradpour_encode_chipseq.meme Data/Motifs/zlab_encode_chipseq.meme
		
#meme-chip -oc Results/Chip-seq/$tf/$cl -db Data/Motifs/rap_unijasol_chipseq.meme -db Data/Motifs/unijasjol_chipseq.jasme -db Data/Motifs/unijasjol_chipseq.jolme -db Data/Motifs/unijasjol_chipseq.unime -db Data/Motifs/uniprobe_chipseq.meme -db Data/Motifs/Encode_chipseq.meme -db Data/Motifs/zhao_unijasol_chipseq.meme -nmeme 500 -meme-nmotifs 5 -dreme-m 5 Data/Chip-seq/$tf/$cl.fa
	done
done
