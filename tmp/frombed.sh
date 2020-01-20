for tf in $(cat /home/kipkurui/Project/Motif_Assessment/Data/unijasjol_chipseq.txt);
do
	for cl in $(grep -i $tf /tmp/encode);
	do
		mkdir -p $HOME/Project/Motif_Assessment/Data/Chip-seq/$tf
		cd $HOME/Project/Motif_Assessment/Data/Chip-seq/$tf
		cut -f 1-3 /home/kipkurui/ChIP-seq/wgEncodeAwgTfbs"$cl"UniPk.narrowPeak | bed-widen -width 500 > /tmp/$cl.bed
		sort -u -o $cl.bed /tmp/$cl.bed
		fastaFromBed -fi $HOME/Project/Data/genomes/hg19/hg19.fa -bed $cl.bed -fo $cl.fa
		#python /home/kipkurui/Dropbox/Pythontest/removemasked.py /tmp/$cl.fa $cl.fa
		
done
done
