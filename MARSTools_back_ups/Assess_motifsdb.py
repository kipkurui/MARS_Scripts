"""
Assess_motif is a python module for assessing the motif Performance using
ChIP-seq and PBM data.
z
Requires:
    -> A motif file in meme format\n
    -> ChIP-seq file (for now) in tab delimited format Chr score sequence \n
    -> A motif scoring framework to use:\n
        -gomeroccupancyscore\n
        -sumoccupancyscore\n
        -maxoccuopancyscore\n
        -energyscore\n
    -> Output file
 Example :
     python Assess_motifs.py Max-jas.meme Haib-K562-Max-50.posneg gomeroccupancy score Haib-K562-Max-50
"""

#!/usr/bin/python
import sys
from math import log
from math import exp
import numpy as np
from sklearn import metrics
import os
import matplotlib.pyplot as pl

#
##TODO: Get all the functions working well and in generic forms and then
##separate out ChIP-seq and PBM main programs


def readpwm(motiffile):
    """
    Read the PWM and convert it into a dictionary

    This function has been replaced by Getmotif
    """
    found = 0
    nrows = 0
    n = 1
    areapwm = {}
    areapwm["A"] = []
    areapwm["C"] = []
    areapwm["G"] = []
    areapwm["T"] = []
    with open(motiffile, "r") as motiffile:
        for line in motiffile:
            words = line.split()
            if found == 0:
                if line.startswith("MOTIF"):
                    if len(words) < 3:
                        words.append("")
                    found = 1
                    continue
            if found == 1:
                if line.startswith("letter-probability"):
                    nrows = int((line.split("w="))[1].split()[0])
                    found = 2
                continue
            if found == 2:
                #check = 0
                if n < nrows + 1:
                    print words[0],
                    areapwm["A"].append(float(words[0]))
                    areapwm["C"].append(float(words[1]))
                    areapwm["G"].append(float(words[2]))
                    areapwm["T"].append(float(words[3]))
                    n += 1
        return areapwm, nrows


def getmotif(meme, motif="MOTIF"):
    """
    Extract a motif from meme file given a unique motif
    name and create dictionary for sequence scoring

    Default motif name is keyword MOTIF for single motif files. 
    """

    areapwm = {}
    areapwm["A"] = []
    areapwm["C"] = []
    areapwm["G"] = []
    areapwm["T"] = []
    flag = 0
    check = 0
    with open(meme, "r") as f1:
        for line in f1:
            if str(motif) in line:
                flag += 1
            if "letter-probability" in line and flag == 1:
                w = line.split(" ")[5]
                flag += 1
                continue
            if flag == 2 and int(check) < int(w):
                if line == "\n":
                    continue
                else:
                    words = line.split()
                    areapwm["A"].append(float(words[0]))
                    areapwm["C"].append(float(words[1]))
                    areapwm["G"].append(float(words[2]))
                    areapwm["T"].append(float(words[3]))
                    check += 1
        return areapwm, int(w)


def rcpwm(areapwm, pwmlen):
    """
    Takes as input the forward pwm and returns a reverse
    complement of the motif
    """

    rcareapwm = {}
    rcareapwm["A"] = []
    rcareapwm["C"] = []
    rcareapwm["G"] = []
    rcareapwm["T"] = []
    for i in range(pwmlen):
        rcareapwm["A"].append(areapwm["T"][pwmlen - i - 1])
        rcareapwm["C"].append(areapwm["G"][pwmlen - i - 1])
        rcareapwm["G"].append(areapwm["C"][pwmlen - i - 1])
        rcareapwm["T"].append(areapwm["A"][pwmlen - i - 1])
    return rcareapwm

#############################################
##Scoring functions
###########################################

##TODOO ##Convert them into classs objects


def gomeroccupancyscore(areapwm, pwmlen, seq):
    """
    Takes as input a PWM dictionary, and a sequences and
    computes the gomer score
    """

    value_gomer = 1
    areapwm_rc = rcpwm(areapwm, pwmlen)
    for i in range(pwmlen - 1, 1, -1):
        prod_gomer = 1
        prod_gomer_rc = 1
        for j in range(pwmlen):
            if j <= i:
                prod_gomer *= 0.25
                prod_gomer_rc *= 0.25
            elif (j + i) > len(seq):
                prod_gomer *= 0.25
                prod_gomer_rc *= 0.25
            else:
                s = seq[j + i]
                prod_gomer *= areapwm[s][j]
                prod_gomer_rc *= areapwm_rc[s][j]
        value_gomer *= (1 - prod_gomer) * (1 - prod_gomer_rc)
    for i in range(len(seq) - 1):
        prod_gomer = 1
        prod_gomer_rc = 1
        for j in range(pwmlen - 1):
            if (j + i) >= len(seq):
                prod_gomer *= 0.25
                prod_gomer_rc *= 0.25
            else:
                prod_gomer *= areapwm[seq[j + i]][j]
                prod_gomer_rc *= areapwm_rc[seq[j + i]][j]
        value_gomer *= (1 - prod_gomer) * (1 - prod_gomer_rc)
    value_gomer = 1 - value_gomer
    return value_gomer


def sumoccupancyscore(areapwm, pwmlen, seq):
    """
    Takes as input a PWM dictionary, and a sequences and
    computes the sum occupancy score.
    """

    value_gomer = 0
    areapwm_rc = rcpwm(areapwm, pwmlen)
    for i in range(len(seq) - 1):
        prod_gomer = 1
        prod_gomer_rc = 1
        for j in range(pwmlen - 1):
            if (j + i) >= len(seq):
                prod_gomer *= 0.25
                prod_gomer_rc *= 0.25
            elif seq[j + i] not in ["A", "C", "G", "T"]:
                prod_gomer *= 0.25
                prod_gomer_rc *= 0.25
            else:
                prod_gomer *= areapwm[seq[j + i]][j]
                prod_gomer_rc *= areapwm_rc[seq[j + i]][j]
        value_gomer += prod_gomer  # +(prod_gomer_rc)
    # value_gomer = 1 - value_gomer
    return value_gomer
    #log(0.17/0.25) / log(2) x 100


def sumlogoddsscore(areapwm, pwmlen, seq):
    """
    Takes as input a PWM dictionary, and a sequences and
    computes the sum of the log odds scores.
    
    This is the scoring approach that is used by meme

    We might  also need to check the effect of using Max score
    """

    value_gomer = 1
    areapwm_rc = rcpwm(areapwm, pwmlen)
    for i in range(len(seq) - 1):
        prod_gomer = 0
        prod_gomer_rc = 0
        for j in range(pwmlen - 1):
            if (j + i) >= len(seq):
                prod_gomer += 0.0
                prod_gomer_rc += 0.0
            elif seq[j + i] not in ["A", "C", "G", "T"]:
                prod_gomer += 0.0
                prod_gomer_rc += 0.0
            else:
                q = areapwm[seq[j + i]][j]
                if q == 0:
                    q = 0.000000000000000000000000000001  # make this as close to zero as possible
                else:
                    q = areapwm[seq[j + i]][j]
                prod_gomer += (np.log(q/0.25) / np.log(2)) * 100
                #prod_gomer_rc += (np.log(q/0.25) / np.log(2)) * 100
                #prod_gomer_rc += areapwm_rc[seq[j + i]][j]
        value_gomer += prod_gomer  #+(prod_gomer_rc)
    #value_gomer = 1 - value_gomer
    return value_gomer


def maxlogoddsscore(areapwm, pwmlen, seq):
    """
    Takes as input a PWM dictionary, and a sequences and
    computes the max of the log odds scores.
    
    This is the scoring approach that is used by meme

    We might  also need to check the effect of using Max score
    """

    value_gomer = 1
    prod = []
    areapwm_rc = rcpwm(areapwm, pwmlen)
    for i in range(len(seq) - 1):
        prod_gomer = 0
        prod_gomer_rc = 0
        for j in range(pwmlen - 1):
            if (j + i) >= len(seq):
                prod_gomer += 0.0
                prod_gomer_rc += 0.0
            elif seq[j + i] not in ["A", "C", "G", "T"]:
                prod_gomer += 0.0
                prod_gomer_rc += 0.0
            else:
                q = areapwm[seq[j + i]][j]
                if q == 0:
                    q = 0.000000000000000000000000000001  # make this as close to zero as possible
                else:
                    q = areapwm[seq[j + i]][j]
                prod_gomer += (np.log(q/0.25) / np.log(2)) * 100
                #prod_gomer_rc += (np.log(q/0.25) / np.log(2)) * 100
                #prod_gomer_rc += areapwm_rc[seq[j + i]][j]
          #+(prod_gomer_rc)
        prod = [prod_gomer]
    value_gomer = max(prod)
    return value_gomer


def maxoccupancyscore(areapwm, pwmlen, seq):
    """
    Takes as input a PWM dictionary and a sequences to
    compute the maximum  occupancy score.
    """
    prod = []
    areapwm_rc = rcpwm(areapwm, pwmlen)
    for i in range(len(seq) - 1):
        prod_gomer = 1
        prod_gomer_rc = 1
        #store site  scores
        for j in range(pwmlen - 1):
            if (j + i) >= len(seq):
                prod_gomer *= 0.25
                prod_gomer_rc *= 0.25
            else:
                prod_gomer *= areapwm[seq[j + i]][j]
                prod_gomer_rc *= areapwm_rc[seq[j + i]][j]
        prod = [prod_gomer + prod_gomer_rc]
    value_gomer = max(prod)
    return value_gomer


def amaoccupancyscore(areapwm, pwmlen, seq):
    """
    Score sequences using AMA scoring uses average occupancy scores

    This is still to be completed...
    """

    areapwm_rc = rcpwm(areapwm, pwmlen)
    for i in range(len(seq) - 1):
        prod_gomer = 1
        prod_gomer_rc = 1
        for j in range(pwmlen - 1):
            if (j + i) >= len(seq):
                prod_gomer *= 0.25
                prod_gomer_rc *= 0.25
            else:
                prod_gomer *= areapwm[seq[j + i]][j]
                prod_gomer_rc *= areapwm_rc[seq[j + i]][j]
        prod = [(prod_gomer) + (prod_gomer_rc)]
    value_gomer = sum(prod) / len(prod)
    return value_gomer


def energyscore(areapwm, pwmlen, seq):
    """
    Score sequences using the beeml energy scoring approach.

    Borrowed greatly from the work of Zhao and Stormo

    P(Si)=1/(1+e^Ei-u)

    Ei=sumsum(Si(b,k)e(b,k))

    Previous approaches seem to be using the the minimum sum of the
    energy contribution of each of the bases of a specific region.

    This is currently showing some promise but further testing is
    needed to ensure that I have a robust algorithm.
    """

    prod = []
    areapwm_rc = rcpwm(areapwm, pwmlen)
    for i in range(len(seq) - 1):
        prod_gomer = 0
        prod_gomer_rc = 0
        for j in range(pwmlen - 1):
            #fwrev = []
            if (j + i) >= len(seq):
                prod_gomer += 0.25
                prod_gomer_rc += 0.25
            else:
                prod_gomer += areapwm[seq[j + i]][j]
                prod_gomer_rc += areapwm_rc[seq[j + i]][j]
                #fwrev = [prod_gomer, prod_gomer_rc]
            prod.append(1 / (1 + (exp(prod_gomer))))
    value_gomer = min(prod)
    return value_gomer

######################################################################
##      Assessment metrics
######################################################################
#Convert all these assessment metrics into a class object and do the same
#for the rest of the  functions to respective categories of classes.


def computeauc(preds, cutoff, label):
    """
    Compute Area under ROC curve given a prediction file formatted as floows:
        seq_name                \t Score           \t Classification \n
        chr2:43019807-43019857	\t 4.251985e-06	\t 1 \n
        . \n
        . \n
        . \n
        chr2:43619807-43619857	\t 4.251985e-08	\t 0 \n
    log: Changed to be able to compute AUC for any scores
    """
    l = []
    for i in range(len(preds)):
        if i < cutoff:
            l.append(1)
        else:
            l.append(0)
    pr = np.array(preds)
    y = np.array(l)
    fpr, tpr, thresholds = metrics.roc_curve(y, pr, pos_label=label)
    AUC = metrics.auc(fpr, tpr)
    prin = "AUC score (sklearn): \t %.3f" % (AUC)
    return AUC
    
    ##TODO Need to find a way of reversing the labels when the Beeml is being used
    ##TODO to score the sequences since lower energy is better. Resolved \


def ROC_AUC(preds, cutoff, label):
    """
    This is the AUC computation adopted from Clarke 2003

    Want to compare the performance similarity with the sklearn
    algorithm (both do give same results)
    """

    from scipy.stats import stats
    from numpy import mean, array, hstack
    if label == 1:
        fg_vals = preds[:cutoff]
        bg_vals = preds[cutoff:]
    else:
        fg_vals = preds[cutoff:]
        bg_vals = preds[0:cutoff]
    if len(fg_vals) == 0 or len(bg_vals) == 0:
        return None

    fg_len = len(fg_vals)
    total_len = len(fg_vals) + len(bg_vals)

    if type(fg_vals) != type(array([])):
        fg_vals = array(fg_vals)
    if type(bg_vals) != type(array([])):
        bg_vals = array(bg_vals)

    fg_rank = stats.rankdata(fg_vals)
    #print fg_rank
    total_rank = stats.rankdata(hstack((fg_vals, bg_vals)))
    auc = (sum(total_rank[:fg_len]) - sum(fg_rank)) / (fg_len * (total_len - fg_len))
    #print "AUC score (Clarke): \t %.3f" % auc
    return auc
    #return (sum(total_rank[:fg_len]) - sum(fg_rank))/ (fg_len * (total_len - fg_len))


def MNCP(preds, cutoff, label):
    """
    This is the MNCP computation adopted from Clarke 2003

    MNCP is a rank based metric similar to AUC but its a plot of TP and all positives
    hence considered to be less affected by false positives.

    MNCP is the mean normalized
    """
    from scipy.stats import stats
    from numpy import mean, array, hstack
    if label == 1:
        fg_vals = preds[:cutoff]
        bg_vals = preds[cutoff:]
    else:
        fg_vals = preds[cutoff:]
        bg_vals = preds[:cutoff]
    fg_len = len(fg_vals)
    total_len = len(fg_vals) + len(bg_vals)

    if type(fg_vals) != type(array([])):
        fg_vals = array(fg_vals)
    if type(bg_vals) != type(array([])):
        bg_vals = array(bg_vals)

    fg_rank = stats.rankdata(fg_vals)
    total_rank = stats.rankdata(hstack((fg_vals, bg_vals)))
    slopes = []
    for i in range(len(fg_vals)):
        slope = ((fg_len - fg_rank[i] + 1) / fg_len) / ((total_len - total_rank[i] + 1) / total_len)
        slopes.append(slope)
    #print "MNCP score (Clarke): \t%.3f" % (mean(slopes))
    #print mean(slopes)
    return mean(slopes)


def computeaps(preds, mot):
    """
    Compute Area under precision curve given a prediction file formatted as follows:
        seq_name                \t Score           \t Classification \n
        chr2:43019807-43019857	\t 4.251985e-06	\t 1 \n
        . \n
        . \n
        . \n
        chr2:43619807-43619857	\t 4.251985e-08	\t 0 \n
    """
    a = []
    pred = []
    with open(preds) as debru:
        for line in debru:
            a += [float(line.split()[2])]
            pred += [float(line.split()[1])]
    """pr = np.array(pred)
    y = np.array(a)
    APS = aps(y, pr)
    prin = "APS score: \t %.3f" % (APS)
    return APS"""


def compute_spearman(obs, pred):
    """
    Compute Spearman correlation in the positive set of PBM data.

    Given pbm data, we need to compute the correlation of the predicted
    score and the intensity score
    """
    from scipy import stats

    #obs = obs[:getpbmpossitives(pbm)]
    #pred = pred[:getpbmpossitives(pbm)]
    spear = stats.spearmanr(obs[:500], pred[:500])[0]
    pr = ("Spearmansscore: \t%f" % spear)
    return spear


def compute_pearson(obs, pred):
    """
    Compute pearson correlation in the positive set of data.

    Given data, we need to compute the correlation of the predicted
    score and the existing score
    """
    from scipy import stats
    #print obs
    obs=[float(i) for i in obs]
    pred=[float(i) for i in pred]
    obs=[float(i)-np.mean(obs) for i in obs]
    pred=[float(i)-np.mean(pred) for i in pred]
    #obs = obs[:getpbmpossitives(pbm)]
    #pred = pred[:getpbmpossitives(pbm)]
    pear = stats.pearsonr(obs[:500], pred[:500])[0]
    pr = ("pearsonsscore: \t%f" % pear)
    return pear


##TODO: Find a way to speed up PBM scoring and analysis. Generally it
##is because I am using too
##TODO: much data for the negative set.

#######################################################
##          Assess motifs using PBM data
#########################################################


def getpbmpossitives(pbm):
    """
    Normalize Pbm data to extract positive and negative probes
    for AUC computation.

    Borrowed greatly from Cheng 2007

    The PBM intensity score are transformed to positive sequences by adding
    the minimum score+1 to each score then log transformed.

    The positive probes are generally, 4*MAD(mean absolute deviation) above
    the median of the normalized scores.

    The rest are considered as the negative probes but may revise this.


    """
    ##TODO: Include an option for optional printout of the output
    intensity = []
    sequences = []
    positive = []
    with open(pbm) as probes:
        for line in probes:
            intensity.append(float(line.split()[0]))
            sequences.append(line.split()[1])
        intensity = np.array(intensity)
        Min = min(intensity) - 1
        intensity += (-Min)
        transformed = np.log(intensity)
    Median = np.median(transformed) + (0.6745 * 4)
    for i in transformed:
        if i > Median:
            positive.append(i)
    return len(positive)
    #with open(output,"w") as out:
    #for j in range(len(positive)-1):
    #wr="%f\t%s\n"% (positive[j],sequences[j])
    #out.write(wr)


def scorepbm(pbm, scoref, motfile, motif="MOTIF"):
    """
    Given the PBM probe file, this function scores the sequences
    using the specified scoring function
    """
    pbmscore = []
    intensity = []
    pos = getpbmpossitives(pbm)
    os.system("head -%i %s >/tmp/pos.txt" % (pos, pbm))
    os.system("tail -%i %s >>/tmp/pos.txt" % (pos, pbm))
    #with open(outfile, "w") as out:
    with open("/tmp/pos.txt") as debru:
        for line in debru:
            details = line.split()
            seq = details[1]
            areapwm = getmotif(motfile, motif)[0]
            pwmlen = getmotif(motfile, motif)[1]
            seqscore = scoref(areapwm, pwmlen, seq[:36])
            pbmscore.append(seqscore)
            intensity.append(details[0])
            #wr = "%s\t%e\t%s\n" % (details[1], seqscore, seq)
            #out.write(wr)
    return intensity, pbmscore


def assesspbm():
    """
    This function runs the analysis of PBM data
    """


##TODO Change performance metric functional to a generic format
#that can take as input any data irrespective of the source

############################################################
## Assess motifs using ChIP-seq data
############################################################


def scorechipseq(chipseq, scoref, ID):
    """
    Give the ChIP-seq file, this function scores the sequences
    using the specified scoring function.
    """

    seqscore=[]
    chipscore=[]
    ou="/tmp/out-%s" % scoref
    with open(chipseq) as chip:
        with open(ou, "w") as out:
            for line in chip:
                details = line.split()
                areapwm = forasseessmotifs(ID)[0]
                pwmlen = len(areapwm['A'])
                out.write("%f \t %s \n" % (scoref(areapwm, pwmlen, details[2]),(details[1])))
                seqscore.append(scoref(areapwm, pwmlen, details[2]))
                chipscore.append(details[1])
            return chipscore, seqscore
  
##TODO Check if I can put an option to write the output to file.
    
#############################################################
# Execute the above functions using user specified details #
#############################################################

##TODO #Need to implement argparse at this point to handle the commandline #options to the program


################################################################################
#Read data from the databasese
################################################################################

import MySQLdb
mydb= MySQLdb.connect(host='localhost',
		      user = 'root',
		      passwd ='Ch11t0',
		      db = 'Caleb_db')
cur=mydb.cursor()
#command = cur.execute("SELECT * FROM JASPAR_test.MATRIX WHERE lower(NAME) LIKE '%max%' AND lower(NAME) NOT LIKE '%:%' ")
#IDS=cur.fetchall()


def getmotifdb(ID,tf):
    wrout=open("/tmp/"+tf,"a")
    statement="SELECT * FROM MATRIX_DATA WHERE ID='%i'" % ID
    command = cur.execute(statement)
    A=[]
    C=[]
    G=[]
    T=[]
    row=cur.fetchall()
    for i in row:
        if i[1]=="A":
            A.append(i[3])
        elif i[1]=="C":
            C.append(i[3])
        elif i[1]=="G":
            G.append(i[3])
        else:
            T.append(i[3])
    statement="SELECT * FROM MATRIX WHERE ID='%s'" % ID
    command = cur.execute(statement)
    det=cur.fetchall()
    mid=det[0][1]+"."+str(det[0][4])
    name=det[0][2]

    motif=("\nMOTIF %s %s\n\n"% (mid,name))
    wrout.write(motif)
    
    l=A[0]+C[0]+G[0]+T[0] #use this to convert to PWM from PFM
    header=("letter-probability matrix: alength= 4 w= %s nsites= 20 E= 0\n" % (str(len(A)-1)))
    wrout.write(header)
    for i in range(1,len(A)):
        out=("  %.6f\t  %.6f\t  %.6f\t  %.6f\t\n" %(A[i],float(C[i]),float(G[i]),float(T[i])))
        wrout.write(out)

def forasseessmotifs(ID):
    statement="SELECT * FROM MATRIX_DATA WHERE ID='%i'" % ID
    command = cur.execute(statement)
    areapwm = {}
    areapwm["A"] = []
    areapwm["C"] = []
    areapwm["G"] = []
    areapwm["T"] = []
    row=cur.fetchall()
    for i in row:
        if i[1]=="A":
            areapwm["A"].append(i[3])
        elif i[1]=="C":
            areapwm["C"].append(i[3])
        elif i[1]=="G":
            areapwm["G"].append(i[3])
        else:
            areapwm["T"].append(i[3])
    statement="SELECT * FROM MATRIX WHERE ID='%i'" % ID
    command = cur.execute(statement)
    det=cur.fetchall()
    mid=det[0][1]+"."+str(det[0][4])
    name=det[0][4]
    #print mid,
    return areapwm,mid

def getchip(ID):
    statement="SELECT CHIP_ID FROM CHIP_SEQ WHERE TF_ID='%s'" % ID
    command = cur.execute(statement)
    tfid=cur.fetchall()
    statement="SELECT AT_100 FROM CHIP_DATA WHERE CHIP_ID='%i'" % tfid[0]
    command = cur.execute(statement)
    tfid=cur.fetchall()
    tflist=[]
    for i in tfid:
        tflist.append(i[0])
    #print tflist
    return tflist
    
 
##############################################################################   
# Run the complete assess progrm as a function
#################################################################################

def runassess(tflist,sc,output,ID):
    """
    Use this fuction to summarize the data for each and
    every TF and scoring function
    """
    auc_max=[]
    spearman_max=[]
    mncp_max=[]
    enrich_max=[]
    TFs=[]
    pr=[]
    with open(output, "a") as out:
        for i in sc:
            auc = []
            spearman = []
            mncp = []
            enrich=[]
            pearson = []
            score = eval(i)#sys.argv[3])
            if i == "energyscore":  # cater for reverse labels in energy score sys.argv[3]
                label = 0
            else:
                label = 1
            for c in tflist:
                chipseq = "/home/kipkurui/Project/Motif_Assessment/Data/ChIP-seq/Derived/Posneg/%s" % c
                schip = scorechipseq(chipseq, score, ID)
                auc += [computeauc(schip[1], 500, label)]
                mncp += [MNCP(schip[1], 500, label)]
                #enrich+=[max_enrichment(schip[1], 500, label)[1]]
                spearman += [compute_spearman(schip[0], schip[1])]
                pearson += [compute_pearson(schip[0], schip[1])]
            motif=forasseessmotifs(ID)[1]
            pr = "%s \t %f \t %f \t %f \t %f \n" % (motif, np.mean(auc), np.mean(pearson), np.mean(spearman), np.mean(mncp))
            auc_max += [np.mean(auc)]
            spearman_max += [np.mean(spearman)]
            mncp_max += [np.mean(mncp)]
            out.write(pr)     

sc=["gomeroccupancyscore"] #other scoring approaches can also be specified
tf="Myc"
output="/home/kipkurui/Tests/%s.gomer" % tf
command = cur.execute("SELECT TF_ID FROM MATRIX WHERE lower(MOTIF_NAME) LIKE '%myc%' AND lower(MOTIF_NAME) NOT LIKE '%:%' ")
IDs=[]
a=cur.fetchall()
ID=a[0][0]
tflist=getchip(ID) 
command = cur.execute("SELECT ID FROM MATRIX WHERE TF_ID LIKE '%s'" % ID)
a=cur.fetchall()
for i in a:
    IDs.append(i[0])
for i in IDs:
    getmotifdb(i,tf)
    #runassess(tflist,sc,output,i)
    #what=forasseessmotifs(i)


