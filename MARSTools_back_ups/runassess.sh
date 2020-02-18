#trying to make this a single run. 

###JASPAR

##BEEML
grep "MOTIF" Data/Motifs/zhao_unijasol_chipseq.meme |cut -f2 -d" " >/tmp/zhao.txt
for t in Egr1 #$(cat /home/kipkurui/Project/Motif_Assessment/Data/unijasjol_chipseq.txt);
do
	
	#for tf in $(grep -i "$t" /tmp/zhao.txt); #/home/kipkurui/Project/Motif_Assessment/Data/uniprobeforun.txt
	#do
	echo $tf >>/tmp/Egr1zhao.txt #/home/kipkurui/Project/Motif_Assessment/Results/Zhao_chipseq_scores.txt
	for score in gomeroccupancyscore energyscore maxoccupancyscore sumoccupancyscore sumlogoddsscore
	do
		for file in $(ls /home/kipkurui/Project/Motif_Assessment/Data/ENCODETFs/$t)
		do
			#mot=`python /home/kipkurui/Project/Motif_Assessment/convertuniprobe.py $tf`
			#echo $mot
			python /home/kipkurui/Dropbox/Pythontest/Example/Assess_motifs.py Data/Motifs/zhao_unijasol_chipseq.meme Data/ENCODETFs/$t/$file $score $tf >>/tmp/Egr1zhao.txt #/home/kipkurui/Project/Motif_Assessment/Results/Zhao_chipseq_scores.txt

		done
	done
done
done
: <<'END'


for tf in $(cat /home/kipkurui/Project/Motif_Assessment/Data/unijasjol_chipseq.txt);
do
	echo $tf >>/home/kipkurui/Project/Motif_Assessment/Results/Jaspar_chipseq_scores.txt
	for score in gomeroccupancyscore energyscore maxoccupancyscore sumoccupancyscore sumlogoddsscore
	do
		for file in $(ls /home/kipkurui/Project/Motif_Assessment/Data/ENCODETFs/$tf)
		do
			mot=`python /home/kipkurui/Project/Motif_Assessment/convertjaspar.py $tf`
			#echo $mot
			python /home/kipkurui/Dropbox/Pythontest/Example/Assess_motifs.py Data/Motifs/unijasjol_chipseq.jasme Data/ENCODETFs/$tf/$file $score $mot >>/home/kipkurui/Project/Motif_Assessment/Results/Jaspar_chipseq_scores.txt

		done
	done
done

###JOLMA


for t in $(cat /home/kipkurui/Project/Motif_Assessment/Data/unijasjol_chipseq.txt);
do
	for tf in $(grep -i "$t" /home/kipkurui/Project/Motif_Assessment/Data/Jolma2013);
	do
	echo $tf >>/home/kipkurui/Project/Motif_Assessment/Results/Jolma_chipseq_scores.txt

	for score in gomeroccupancyscore energyscore maxoccupancyscore sumoccupancyscore sumlogoddsscore
	do 
		for file in $(ls /home/kipkurui/Project/Motif_Assessment/Data/ENCODETFs/$t)
		do
			python /home/kipkurui/Dropbox/Pythontest/Example/Assess_motifs.py Data/Motifs/unijasjol_chipseq.jolme Data/ENCODETFs/$t/$file $score $tf >>/home/kipkurui/Project/Motif_Assessment/Results/Jolma_chipseq_scores.txt

		done
	done
done
done
#END
##UNIPROBE

for t in $(cat /home/kipkurui/Project/Motif_Assessment/Data/unijasjol_chipseq.txt);
do
	for tf in $(grep -i "$t" /home/kipkurui/Project/Motif_Assessment/Data/uniprobeforun.txt); #/home/kipkurui/Project/Motif_Assessment/Data/uniprobeforun.txt
	do
	echo $tf >>/home/kipkurui/Project/Motif_Assessment/Results/Uniprobe_chipseq_scores.txt
	for score in gomeroccupancyscore energyscore maxoccupancyscore sumoccupancyscore sumlogoddsscore
	do
		for file in $(ls /home/kipkurui/Project/Motif_Assessment/Data/ENCODETFs/$t)
		do
			mot=`python /home/kipkurui/Project/Motif_Assessment/convertuniprobe.py $tf`
			#echo $mot
			python /home/kipkurui/Dropbox/Pythontest/Example/Assess_motifs.py Data/Motifs/unijasjol_chipseq.unime Data/ENCODETFs/$t/$file $score $mot >>/home/kipkurui/Project/Motif_Assessment/Results/Uniprobe_chipseq_scores.txt

		done
	done
done
done



#run RAP
for t in $(cat /home/kipkurui/Project/Motif_Assessment/Data/unijasjol_chipseq.txt);
do
	for tf in $(grep -i "$t" /home/kipkurui/Project/Motif_Assessment/Data/rapmotifs.txt); #/home/kipkurui/Project/Motif_Assessment/Data/uniprobeforun.txt
	do
	echo $tf >>/home/kipkurui/Project/Motif_Assessment/Results/Rap_chipseq_scores.txt
	for score in gomeroccupancyscore energyscore maxoccupancyscore sumoccupancyscore sumlogoddsscore
	do
		for file in $(ls /home/kipkurui/Project/Motif_Assessment/Data/ENCODETFs/$t)
		do
			python /home/kipkurui/Dropbox/Pythontest/Example/Assess_motifs.py Data/Motifs/rap_unijasol_chipseq.meme Data/ENCODETFs/$t/$file $score $tf >>/home/kipkurui/Project/Motif_Assessment/Results/Rap_chipseq_scores.txt

		done
	done
done
done


#zlab
for t in $(cat /home/kipkurui/Project/Motif_Assessment/Data/unijasjol_chipseq.txt);
do
	for tf in $(grep -i "$t" /tmp/encode); #/home/kipkurui/Project/Motif_Assessment/Data/uniprobeforun.txt
	do
	echo $tf >>/home/kipkurui/Project/Motif_Assessment/Results/Encode_chipseq_scores.txt
	for score in gomeroccupancyscore energyscore maxoccupancyscore sumoccupancyscore sumlogoddsscore
	do
		for file in $(ls /home/kipkurui/Project/Motif_Assessment/Data/ENCODETFs/$t)
		do
			python /home/kipkurui/Dropbox/Pythontest/Example/Assess_motifs.py /home/kipkurui/Project/Motif_Assessment/Data/Motifs/Encode_chipseq.meme Data/ENCODETFs/$t/$file $score $tf >>/home/kipkurui/Project/Motif_Assessment/Results/Encode_chipseq_scores.txt

		done
	done
done
done
#END
for t in $(cat /home/kipkurui/Project/Motif_Assessment/Data/unijasjol_chipseq.txt);
do
	for tf in $(grep -i "$t" /home/kipkurui/Project/Motif_Assessment/Data/pourlist.txt |cut -f1 -d" "); #/home/kipkurui/Project/Motif_Assessment/Data/uniprobeforun.txt
	do
	echo $tf >>/home/kipkurui/Project/Motif_Assessment/Results/Pour_Encode_chipseq_scores.txt
	for score in gomeroccupancyscore energyscore maxoccupancyscore sumoccupancyscore sumlogoddsscore
	do
		for file in $(ls /home/kipkurui/Project/Motif_Assessment/Data/ENCODETFs/$t)
		do
			python /home/kipkurui/Dropbox/Pythontest/Example/Assess_motifs.py /home/kipkurui/Project/Motif_Assessment/Data/Motifs/Kheradpour_encode_chipseq.meme Data/ENCODETFs/$t/$file $score $tf >>/home/kipkurui/Project/Motif_Assessment/Results/Pour_Encode_chipseq_scores.txt

		done
	done
done
done


##HUMOCO

for t in $(cat /home/kipkurui/Project/Motif_Assessment/Data/unijasjol_chipseq.txt);
do
	for tf in $(grep -i "$t" /home/kipkurui/Project/Motif_Assessment/Data/hucomoco.txt); #/home/kipkurui/Project/Motif_Assessment/Data/uniprobeforun.txt
	do
	echo $tf >>/home/kipkurui/Project/Motif_Assessment/Results/Humoco_Encode_chipseq_scores.txt
	for score in gomeroccupancyscore energyscore maxoccupancyscore sumoccupancyscore sumlogoddsscore
	do
		for file in $(ls /home/kipkurui/Project/Motif_Assessment/Data/ENCODETFs/$t)
		do
			python /home/kipkurui/Dropbox/Pythontest/Example/Assess_motifs.py /home/kipkurui/Project/Motif_Assessment/Data/Motifs/hucomoco_chipseq.meme Data/ENCODETFs/$t/$file $score $tf >>/home/kipkurui/Project/Motif_Assessment/Results/Humoco_Encode_chipseq_scores.txt

		done
	done
done
END
done
