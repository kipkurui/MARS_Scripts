#Above parameters can be changed to pick a different file
#Might need to make this into a python script and make it as generic as possible
function Getfasta {
				cut -f 1,2,3 $inputPEAK | bed-widen -width $siz >$cl-$tf-$siz.bed
				if [ -s "$cl-$tf-$siz.bed" ]; then
				
					#Extract the negative sequence
					python extractnegative.py $cl-$tf-$siz.bed $cl-$tf-$siz.negbed 500
					fastaFromBed -tab -fi $hg -bed $cl-$tf-$siz.bed -fo $cl-$tf-$siz.fa
					fastaFromBed -tab -fi $hg -bed $cl-$tf-$siz.negbed -fo $cl-$tf-$siz.negfa
					
					# Prepare the fasta seqences 
					cut -f 7 $inputPEAK >/tmp/f1
					cut -f 1 $cl-$tf-$siz.fa >/tmp/f2
					cut -f 2 $cl-$tf-$siz.fa >/tmp/f3
					paste /tmp/f2 /tmp/f1 /tmp/f3  >/tmp/$cl-$tf-$siz.fasta
					python removemasked.py /tmp/$cl-$tf-$siz.fasta $cl-$tf-$siz.fasta
					rm /tmp/f*

					#Do the same for the negative sequences
					cut -f 7 $inputPEAK >/tmp/f1
					cut -f 1 $cl-$tf-$siz.negfa >/tmp/f2
					cut -f 2 $cl-$tf-$siz.negfa >/tmp/f3
					paste /tmp/f2 /tmp/f1 /tmp/f3  >/tmp/$cl-$tf-$siz.negfasta
					python removemasked.py /tmp/$cl-$tf-$siz.negfasta $cl-$tf-$siz.negfasta
					rm /tmp/f*
					rm $cl-$tf-$siz.negfa
					rm $cl-$tf-$siz.fa
					head -500 $cl-$tf-$siz.fasta >$lab-$cl-$tf-$siz.posneg
					head -500 $cl-$tf-$siz.negfasta >>$lab-$cl-$tf-$siz.posneg
					
				#Cleanup
				rm *fasta
				rm *bed
				fi
	}

hg=$1 #$HOME/Project/Data/genomes/hg19/hg19.fa #change this to your hg path
inputPEAK=$2
cl=$3 #Hepg2
tf=$4 #Max
siz=$5 #100 #the length of the motifs
lab=$6 #lab or the source of the ChIP-seq data
Getfasta