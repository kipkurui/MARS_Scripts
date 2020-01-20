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
import os
from math import exp

import numpy as np
from sklearn import metrics

import MARSTools.Assess_by_score
import MATOM.models
import MATOM.utils
from MATOM import fetch_motif
from MATOM.models import ChipSeq

##TODO: Get all the functions working well and in generic forms and then
##separate out ChIP-seq and PBM main programs

BASE_DIR = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))


def readpwm(motif_file):
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
    with open(motif_file, "r") as motif_file:
        for line in motif_file:
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


def get_motif(meme, motif="MOTIF"):
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
            #FIXME use the split to pik up right name
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


def rc_pwm(area_pwm, pwm_len):
    """
    Takes as input the forward pwm and returns a reverse
    complement of the motif
    """

    rcareapwm = {}
    rcareapwm["A"] = []
    rcareapwm["C"] = []
    rcareapwm["G"] = []
    rcareapwm["T"] = []
    for i in range(pwm_len):
        rcareapwm["A"].append(area_pwm["T"][pwm_len - i - 1])
        rcareapwm["C"].append(area_pwm["G"][pwm_len - i - 1])
        rcareapwm["G"].append(area_pwm["C"][pwm_len - i - 1])
        rcareapwm["T"].append(area_pwm["A"][pwm_len - i - 1])
    return rcareapwm

#############################################
##Scoring functions
###########################################

##TODOO ##Convert them into classs objects


def gomeroccupancyscore(pwm_dictionary, pwm_length, seq):
    """
    Takes as input a PWM dictionary, and a sequences and
    computes the gomer score
    """

    gomer_occupancy = 1
    pwm_dictionary_rc = rc_pwm(pwm_dictionary, pwm_length)
    for i in range(pwm_length - 1, 1, -1):
        prod_gomer = 1
        prod_gomer_rc = 1
        for j in range(pwm_length):
            if j <= i:
                prod_gomer *= 0.25
                prod_gomer_rc *= 0.25
            elif (j + i) > len(seq):
                prod_gomer *= 0.25
                prod_gomer_rc *= 0.25
            else:
                #print "got to else"
                s = seq[j + i]
                prod_gomer *= pwm_dictionary[s][j]
                prod_gomer_rc *= pwm_dictionary_rc[s][j]
        gomer_occupancy *= (1 - prod_gomer) * (1 - prod_gomer_rc)
    for i in range(len(seq) - 1):
        prod_gomer = 1
        prod_gomer_rc = 1
        for j in range(pwm_length - 1):
            if (j + i) >= len(seq):
                prod_gomer *= 0.25
                prod_gomer_rc *= 0.25
            else:
                prod_gomer *= pwm_dictionary[seq[j + i]][j]
                prod_gomer_rc *= pwm_dictionary_rc[seq[j + i]][j]
        gomer_occupancy *= (1 - prod_gomer) * (1 - prod_gomer_rc)
    gomer_occupancy = 1 - gomer_occupancy
    return gomer_occupancy


def sumoccupancyscore(pwm_dictionary, pwm_length, seq):
    """
    Takes as input a PWM dictionary, and a sequences and
    computes the sum occupancy score.
    """

    sum_occupancy = 0
    pwm_dictionary_rc = rc_pwm(pwm_dictionary, pwm_length)
    for i in range(len(seq) - 1):
        occupancy = 1
        occupancy_rc = 1
        for j in range(pwm_length - 1):
            if (j + i) >= len(seq):
                occupancy *= 0.25
                occupancy_rc *= 0.25
            elif seq[j + i] not in ["A", "C", "G", "T"]:
                occupancy *= 0.25
                occupancy_rc *= 0.25
            else:
                occupancy *= pwm_dictionary[seq[j + i]][j]
                occupancy_rc *= pwm_dictionary_rc[seq[j + i]][j]
        sum_occupancy += occupancy + occupancy_rc
    return sum_occupancy/2


def sumlogoddsscore(pwm_dictionary, pwm_length, seq):
    """
    Takes as input a PWM dictionary, and a sequences and
    computes the sum of the log odds scores.

    This is the scoring approach that is used by MEME Suite
    """

    sum_log_odds = 1
    pwm_dictionary_rc = rc_pwm(pwm_dictionary, pwm_length)
    for i in range(len(seq) - 1):
        log_odds = 0
        log_odds_rc = 0
        for j in range(pwm_length - 1):
            if (j + i) >= len(seq):
                log_odds += 0.0
                log_odds_rc += 0.0
            elif seq[j + i] not in ["A", "C", "G", "T"]:
                log_odds += 0.0
                log_odds_rc += 0.0
            else:
                q = pwm_dictionary[seq[j + i]][j]
                q_rc = pwm_dictionary_rc[seq[j + i]][j]
                if q == 0 or (q_rc == 0):
                    q = 0.000000000000000000000000000001  # make this as close to zero as possible
                    q_rc = 0.000000000000000000000000000001
                else:
                    q = pwm_dictionary[seq[j + i]][j]
                    q_rc = pwm_dictionary_rc[seq[j + i]][j]
                log_odds += (np.log(q/0.25) / np.log(2)) * 100
                log_odds_rc += (np.log(q_rc/0.25) / np.log(2)) * 100
        sum_log_odds += log_odds + log_odds_rc
    return sum_log_odds/2


def maxlogoddsscore(pwm_dictionary, pwm_length, seq):
    """
    Takes as input a PWM dictionary, and a sequences and
    computes the max of the log odds scores.

    """

    log_odds_list = []
    pwm_dictionary_rc = rc_pwm(pwm_dictionary, pwm_length)
    for i in range(len(seq) - 1):
        log_odds_score = 0
        log_odds_score_rc = 0
        for j in range(pwm_length - 1):
            if (j + i) >= len(seq):
                log_odds_score += 0.0
                log_odds_score_rc += 0.0
            elif seq[j + i] not in ["A", "C", "G", "T"]:
                log_odds_score += 0.0
                log_odds_score_rc += 0.0
            else:
                q = pwm_dictionary[seq[j + i]][j]
                q_rc = pwm_dictionary_rc[seq[j + i]][j]
                if q == 0 or q_rc == 0:
                    q = 0.000000000000000000000000000001  # make this as close to zero as possible
                    q_rc = 0.000000000000000000000000000001
                else:
                    q = pwm_dictionary[seq[j + i]][j]
                    q_rc = pwm_dictionary_rc[seq[j + i]][j]
            log_odds_score += (np.log(q/0.25) / np.log(2)) * 100
            log_odds_score_rc += (np.log(q_rc/0.25) / np.log(2)) * 100
        log_odds_list.append(log_odds_score)
        log_odds_list.append(log_odds_score)
    max_log_odds = max(log_odds_list)
    return max_log_odds


def maxoccupancyscore(pwm_dictionary, pwm_length, seq):
    """
    Takes as input a PWM dictionary, and a sequences and
    computes the sum occupancy score.
    """
    occupancy_list = []
    pwm_dictionary_rc = rc_pwm(pwm_dictionary, pwm_length)
    for i in range(len(seq) - 1):
        occupancy = 1
        occupancy_rc = 1
        for j in range(pwm_length - 1):
            if (j + i) >= len(seq):
                occupancy *= 0.25
                occupancy_rc *= 0.25
            elif seq[j + i] not in ["A", "C", "G", "T"]:
                occupancy *= 0.25
                occupancy_rc *= 0.25
            else:
                occupancy *= pwm_dictionary[seq[j + i]][j]
                occupancy_rc *= pwm_dictionary_rc[seq[j + i]][j]
        occupancy_list.append(occupancy)
        occupancy_list.append(occupancy_rc)
    max_occupancy = max(occupancy_list)
    return max_occupancy


def amaoccupancyscore(pwm_dictionary, pwm_length, seq):
    """
    Score sequences using AMA scoring uses average occupancy scores

    """
    occupancy_list = []
    pwm_dictionary_rc = rc_pwm(pwm_dictionary, pwm_length)
    for i in range(len(seq) - 1):
        occupancy = 1
        occupancy_rc = 1
        for j in range(pwm_length - 1):
            if (j + i) >= len(seq):
                occupancy *= 0.25
                occupancy_rc *= 0.25
            else:
                occupancy *= pwm_dictionary[seq[j + i]][j]
                occupancy_rc *= pwm_dictionary_rc[seq[j + i]][j]
        occupancy_list.append(occupancy + occupancy_rc)
    ama_occupancy = sum(occupancy_list) / len(occupancy_list)
    return ama_occupancy

##TODO convert the scoring function into class object to reduce with the and also increase the speed


def energyscore(pwm_dictionary, pwm_length, seq):
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

    energy_list = []
    pwm_dictionary_rc = rc_pwm(pwm_dictionary, pwm_length)
    for i in range(len(seq) - 1):
        energy = 0
        energy_rc = 0
        for j in range(pwm_length - 1):
            if (j + i) >= len(seq):
                energy += 0.25
                energy_rc += 0.25
            else:
                energy += pwm_dictionary[seq[j + i]][j]
                energy_rc += pwm_dictionary_rc[seq[j + i]][j]

            energy_list.append(1 / (1 + (exp(energy))))
            energy_list.append(1 / (1 + (exp(energy_rc))))
    energy_score = min(energy_list)
    return energy_score



# def gomeroccupancyscore(areapwm, pwmlen, seq):
#     """
#     Takes as input a PWM dictionary, and a sequences and
#     computes the gomer score
#     """
#
#     value_gomer = 1
#     area_pwm_rc = rc_pwm(areapwm, pwmlen)
#     for i in range(pwmlen - 1, 1, -1):
#         prod_gomer = 1
#         prod_gomer_rc = 1
#         for j in range(pwmlen):
#             if j <= i:
#                 prod_gomer *= 0.25
#                 prod_gomer_rc *= 0.25
#             elif (j + i) > len(seq):
#                 prod_gomer *= 0.25
#                 prod_gomer_rc *= 0.25
#             else:
#                 #print "got to else"
#                 s = seq[j + i]
#                 prod_gomer *= areapwm[s][j]
#                 prod_gomer_rc *= area_pwm_rc[s][j]
#         value_gomer *= (1 - prod_gomer) * (1 - prod_gomer_rc)
#     for i in range(len(seq) - 1):
#         prod_gomer = 1
#         prod_gomer_rc = 1
#         for j in range(pwmlen - 1):
#             if (j + i) >= len(seq):
#                 prod_gomer *= 0.25
#                 prod_gomer_rc *= 0.25
#             else:
#                 prod_gomer *= areapwm[seq[j + i]][j]
#                 prod_gomer_rc *= area_pwm_rc[seq[j + i]][j]
#         value_gomer *= (1 - prod_gomer) * (1 - prod_gomer_rc)
#     value_gomer = 1 - value_gomer
#     return value_gomer
#
#
# def sumoccupancyscore(areapwm, pwmlen, seq):
#     """
#     Takes as input a PWM dictionary, and a sequences and
#     computes the sum occupancy score.
#     """
#
#     value_gomer = 0
#     areapwm_rc = rc_pwm(areapwm, pwmlen)
#     for i in range(len(seq) - 1):
#         prod_gomer = 1
#         prod_gomer_rc = 1
#         for j in range(pwmlen - 1):
#             if (j + i) >= len(seq):
#                 prod_gomer *= 0.25
#                 prod_gomer_rc *= 0.25
#             elif seq[j + i] not in ["A", "C", "G", "T"]:
#                 prod_gomer *= 0.25
#                 prod_gomer_rc *= 0.25
#             else:
#                 prod_gomer *= areapwm[seq[j + i]][j]
#                 prod_gomer_rc *= areapwm_rc[seq[j + i]][j]
#         value_gomer += prod_gomer  # +(prod_gomer_rc)
#     # value_gomer = 1 - value_gomer
#     return value_gomer
#     #log(0.17/0.25) / log(2) x 100
#
#
# def sumlogoddsscore(areapwm, pwmlen, seq):
#     """
#     Takes as input a PWM dictionary, and a sequences and
#     computes the sum of the log odds scores.
#
#     This is the scoring approach that is used by meme
#
#     We might  also need to check the effect of using Max score
#     """
#
#     value_gomer = 1
#     areapwm_rc = rc_pwm(areapwm, pwmlen)
#     for i in range(len(seq) - 1):
#         prod_gomer = 0
#         prod_gomer_rc = 0
#         for j in range(pwmlen - 1):
#             if (j + i) >= len(seq):
#                 prod_gomer += 0.0
#                 prod_gomer_rc += 0.0
#             elif seq[j + i] not in ["A", "C", "G", "T"]:
#                 prod_gomer += 0.0
#                 prod_gomer_rc += 0.0
#             else:
#                 q = areapwm[seq[j + i]][j]
#                 if q == 0:
#                     q = 0.000000000000000000000000000001  # make this as close to zero as possible
#                 else:
#                     q = areapwm[seq[j + i]][j]
#                 prod_gomer += (np.log(q/0.25) / np.log(2)) * 100
#                 #prod_gomer_rc += (np.log(q/0.25) / np.log(2)) * 100
#                 #prod_gomer_rc += areapwm_rc[seq[j + i]][j]
#         value_gomer += prod_gomer  #+(prod_gomer_rc)
#     #value_gomer = 1 - value_gomer
#     return value_gomer
#
#
# def maxlogoddsscore(areapwm, pwmlen, seq):
#     """
#     Takes as input a PWM dictionary, and a sequences and
#     computes the max of the log odds scores.
#
#     """
#
#     value_gomer = 1
#     prod = []
#     areapwm_rc = rc_pwm(areapwm, pwmlen)
#     for i in range(len(seq) - 1):
#         prod_gomer = 0
#         prod_gomer_rc = 0
#         for j in range(pwmlen - 1):
#             if (j + i) >= len(seq):
#                 prod_gomer += 0.0
#                 prod_gomer_rc += 0.0
#             elif seq[j + i] not in ["A", "C", "G", "T"]:
#                 prod_gomer += 0.0
#                 prod_gomer_rc += 0.0
#             else:
#                 q = areapwm[seq[j + i]][j]
#                 if q == 0:
#                     q = 0.000000000000000000000000000001  # make this as close to zero as possible
#                 else:
#                     q = areapwm[seq[j + i]][j]
#                 prod_gomer += (np.log(q/0.25) / np.log(2)) * 100
#                 #prod_gomer_rc += (np.log(q/0.25) / np.log(2)) * 100
#                 #prod_gomer_rc += areapwm_rc[seq[j + i]][j]
#           #+(prod_gomer_rc)
#         prod = [prod_gomer]
#     value_gomer = max(prod)
#     return value_gomer
#
#
# def maxoccupancyscore(areapwm, pwmlen, seq):
#     """
#     Takes as input a PWM dictionary and a sequences to
#     compute the maximum  occupancy score.
#     """
#     prod = []
#     areapwm_rc = rc_pwm(areapwm, pwmlen)
#     for i in range(len(seq) - 1):
#         prod_gomer = 1
#         prod_gomer_rc = 1
#         #store site  scores
#         for j in range(pwmlen - 1):
#             if (j + i) >= len(seq):
#                 prod_gomer *= 0.25
#                 prod_gomer_rc *= 0.25
#             else:
#                 prod_gomer *= areapwm[seq[j + i]][j]
#                 prod_gomer_rc *= areapwm_rc[seq[j + i]][j]
#         prod = [prod_gomer + prod_gomer_rc]
#     value_gomer = max(prod)
#     return value_gomer
#
#
# def amaoccupancyscore(areapwm, pwmlen, seq):
#     """
#     Score sequences using AMA scoring uses average occupancy scores
#
#     This is still to be completed...
#     """
#
#     areapwm_rc = rc_pwm(areapwm, pwmlen)
#     for i in range(len(seq) - 1):
#         prod_gomer = 1
#         prod_gomer_rc = 1
#         for j in range(pwmlen - 1):
#             if (j + i) >= len(seq):
#                 prod_gomer *= 0.25
#                 prod_gomer_rc *= 0.25
#             else:
#                 prod_gomer *= areapwm[seq[j + i]][j]
#                 prod_gomer_rc *= areapwm_rc[seq[j + i]][j]
#         prod = [(prod_gomer) + (prod_gomer_rc)]
#     value_gomer = sum(prod) / len(prod)
#     return value_gomer
#
# ##TODO convert the scoring function into class object to reduce with the and also increase the speed
#
#
# def energyscore(areapwm, pwmlen, seq):
#     """
#     Score sequences using the beeml energy scoring approach.
#
#     Borrowed greatly from the work of Zhao and Stormo
#
#     P(Si)=1/(1+e^Ei-u)
#
#     Ei=sumsum(Si(b,k)e(b,k))
#
#     Previous approaches seem to be using the the minimum sum of the
#     energy contribution of each of the bases of a specific region.
#
#     This is currently showing some promise but further testing is
#     needed to ensure that I have a robust algorithm.
#     """
#
#     prod = []
#     areapwm_rc = rc_pwm(areapwm, pwmlen)
#     for i in range(len(seq) - 1):
#         prod_gomer = 0
#         prod_gomer_rc = 0
#         for j in range(pwmlen - 1):
#             #fwrev = []
#             if (j + i) >= len(seq):
#                 prod_gomer += 0.25
#                 prod_gomer_rc += 0.25
#             else:
#                 prod_gomer += areapwm[seq[j + i]][j]
#                 prod_gomer_rc += areapwm_rc[seq[j + i]][j]
#                 #fwrev = [prod_gomer, prod_gomer_rc]
#             prod.append(1 / (1 + (exp(prod_gomer))))
#     value_gomer = min(prod)
#     return value_gomer

######################################################################
##      Assessment metrics
######################################################################
#Convert all these assessment metrics into a class object and do the same
#for the rest of the  functions to respective categories of classes.

def compute_auc(preds, cutoff, label):
    """
    Compute Area under ROC curve given a prediction file formatted as floows:
        seq_name                \t Score           \t Classification \n
        chr2:43019807-43019857	\t 4.251985e-06	\t 1 \n
        . \n
        . \n
        . \n
        chr2:43619807-43619857	\t 4.251985e-08	\t 0 \n
    log: Changed to be able to compute auc for any scores
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
    auc = metrics.auc(fpr, tpr)
    return auc


def compute_roc_auc(preds, cutoff, label):
    """
    This is the AUC computation adopted from Clarke 2003

    Want to compare the performance similarity with the sklearn
    algorithm (both do give same results)
    """

    from scipy.stats import stats
    from numpy import array, hstack
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


def compute_mncp(preds, cutoff, label):
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
    mncp = mean(slopes)
    return mncp


#FIXME: Fails when the input data lacks scores for correlation. May have to fins a way around or ensure its entered....

def compute_spearman(observed, predicted, cutoff, label=0):
    """
    Compute Spearman correlation in the positive set of PBM data.

    Given pbm data, we need to compute the correlation of the predicted
    score and the intensity score
    :param label:
    """
    from scipy import stats

    spearmans_rank = stats.spearmanr(observed[:cutoff], predicted[:cutoff])[0]

    if label == 0:
        spearmans_rank *= -1

    return spearmans_rank


def compute_pearson(observed, predicted, cutoff, label=1):
    """
    Compute pearson correlation in the positive set of data.

    Given data, we need to compute the correlation of the predicted
    score and the existing score
    """
    from scipy import stats

    observed = [float(i) for i in observed]
    predicted = [float(i) for i in predicted]
    observed = [float(i)-np.mean(observed) for i in observed]
    predicted = [float(i)-np.mean(predicted) for i in predicted]

    pearson = stats.pearsonr(observed[:cutoff], predicted[:cutoff])[0]
    if label == 0:
        pearson *= -1
    return pearson


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


def scorepbm(pbm, scoref, ids, motif="MOTIF"):
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
            area_pwm = fetch_motif.for_assess_motifs(ids)
            pwmlen = len(area_pwm)
            seqscore = scoref(area_pwm, pwmlen, seq[:36])
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


def score_chipseq(chip_seq, score_function, ids, user_motif_details):
    """
    Give the ChIP-seq file, this function scores the sequences
    using the specified scoring function.
    """

    seq_score = []
    chip_score = []
    with open(chip_seq) as chip:
        for line in chip:
            details = line.split()
            if user_motif_details:
                area_pwm = user_motif_details[0]
            else:
                area_pwm = fetch_motif.for_assess_motifs(ids)
            pwm_length = len(area_pwm['A'])

            if len(details) == 2:
                # Use this to cater for the user input BED-file
                pos = 1
            else:
                pos = 2
            # TODO: Why do we need the line below?
            #out.write("%f \t %s \n" % (score_function(area_pwm, pwm_length, details[2]), (details[1])))
            seq_score.append(score_function(area_pwm, pwm_length, details[pos]))
            chip_score.append(details[1])
        return chip_score, seq_score

##TODO Check if I can put an option to write the output to file.
    
#############################################################
# Execute the above functions using user specified details #
#############################################################

##TODO #Need to implement argparse at this point to handle the commandline #options to the program


##############################################################################   
# Run the complete assess program as a function
#################################################################################


def run_assess(sc, output, user_motif_details, tflist, ids=None, tf='', seq_len=100):
    """
    Use this fuction to summarize the data for each and
    every TF and scoring function
    """
    import random
    auc_max = []
    spearman_max = []
    mncp_max = []
    enrich_max = []
    TFs = []
    pr = []

    #Only use 10 random sequences if there are more for visualization and speed (besides they add little value)
    if len(tflist) > 10:
        random.seed(10)
        tflist = random.sample(tflist, 10)
    with open(output, "a") as out:
        for i in sc:
            auc = []
            spearman = []
            mncp = []
            pearson = []
            score = eval(i)
            if i == "energyscore":
                label = 0
            else:
                label = 1
            if user_motif_details:
                motif = user_motif_details[1].strip()
            else:
                motif = ids.motif_id+"-"+ids.collection

            chipseq_path = "%s/Data/ChIP-seq/Derived/Posneg" % BASE_DIR
            write_raw = "MATOM/static/files/%s/%s_raw.%s" % (tf, tf, score_extensions[i])

            #Extract the cell line and lab name from the ChIP-seq file
            with open(write_raw, 'a') as raw_out:
                for b in tflist:
                    if b == '%s_tmp_chip_%s.posneg' % (tf, str(seq_len)):
                        "print is there tmp posneg?"
                        chipseq_path = '%s/MATOM/static/files/%s' % (BASE_DIR, tf)
                        cell_lab = b.split('.')[0]
                        chipseq = "%s/%s" % (chipseq_path, b)

                        schip = score_chipseq(chipseq, score, ids, user_motif_details)
                        cut_off = len(schip[1])/2  # use a flexible cut-off dictated by the sze of the input file
                        # finding a way to output the raw data before averaging)
                        au = compute_auc(schip[1], cut_off, label)
                        auc += [au]
                        mn = compute_mncp(schip[1], cut_off, label)
                        mncp += [mn]
                        #sp = compute_spearman(schip[0], schip[1], cut_off, label)
                        pe = 0
                        sp = 0
                        #spearman += [sp]
                        #pe = compute_pearson(schip[0], schip[1], cut_off, label)
                        #pearson += [pe]
                        pr = "%s\t%s\t%f\t%f\t%f\t%f\n" % (cell_lab, motif, au, pe, sp, mn)
                        raw_out.write(pr)
                    else:
                        cell_lab = b.split('/')[1].split('-')[0]  # b.split('/')[0] to avoid name duplicates, I had to use the ful names
                        chipseq = "%s/%s" % (chipseq_path, b)

                        schip = score_chipseq(chipseq, score, ids, user_motif_details)
                        cut_off = len(schip[1])/2  # use a flexible cut-off dictated by the sze of teh input file
                        # finding a way to output the raw data before averaging)
                        au = compute_auc(schip[1], cut_off, label)
                        auc += [au]
                        mn = compute_mncp(schip[1], cut_off, label)
                        mncp += [mn]
                        sp = compute_spearman(schip[0], schip[1], cut_off, label)
                        spearman += [sp]
                        pe = compute_pearson(schip[0], schip[1], cut_off, label)
                        pearson += [pe]
                        pr = "%s\t%s\t%.3f\t%.3f\t%.3f\t%.3f\n" % (cell_lab, motif, au, pe, sp, mn)
                        raw_out.write(pr)
                pr = "%s\t%s\t%.3f\t%.3f\t%.3f\t%.3f\n" % ('Average', motif, np.mean(auc), np.mean(mncp), np.mean(pearson),
                                                   np.mean(spearman))
                raw_out.write(pr)
            pr = "%s\t%.3f\t%.3f\t%.3f\t%.3f\n" % (motif, np.mean(auc), np.mean(mncp), np.mean(pearson), np.mean(spearman))
            #auc_max += [np.mean(auc)]
            #spearman_max += [np.mean(spearman)]
            #mncp_max += [np.mean(mncp)]
            out.write(pr)


def run_all(tf, sc, user_motif=None, tf_exists=True, uploaded_chip=False, tf_present=False, seq_len=100):
    user_motif_details = None
    score_option = [sc]
    output = "MATOM/static/files/%s/%s.%s" % (tf.lower(), tf, score_extensions[sc])
    raw_out = "MATOM/static/files/%s/%s_raw.%s" % (tf.lower(), tf, score_extensions[sc])
    pr = "%s\t%s\t%s\t%s\t%s\n" % ("Motif", "AUC", "MNCP", "Pearson", "Spearman")
    header = "%s\t%s\t%s\t%s\t%s\t%s\n" % ("Cell_lab", "Motif", "AUC", "MNCP", "Pearson", "Spearman")
    with open(output, "w") as out:
        with open(raw_out, "w") as write_raw:
            out.write(pr)
            write_raw.write(header)
    if tf_exists is True:
  # other scoring approaches can also be specified
        motifs = MATOM.models.get_tf(tf)[1]
        tf_class_id = MATOM.models.get_tf(tf)[0]
        # Here we do not care if the test Chip data exists, as long as
        if uploaded_chip is True:
            tf_list = "USER"
            print "Uploaded chip?"
        else:
            tf_list = MATOM.models.get_chip(tf_class_id)
        if user_motif:
            tf_names = get_tfnames(user_motif)
            #print tf_names
            for mot_name in tf_names:
                user_motif_details = MARSTools.Assess_by_score.get_motif_from_meme(user_motif, mot_name)
                #print user_motif_details[0] This printed the lists out
                run_assess(score_option, output, user_motif_details, tf_list)
                user_motif_details = None
        for mot in motifs:
            run_assess(score_option, output, user_motif_details, tf_list, mot, tf)
    else:

        if tf_present:
            tf_class_id = MATOM.models.get_tf(tf)[0]
            tf_names = get_tfnames(user_motif)
            if ChipSeq.objects.filter(tf_id=tf_class_id).exists():
                tf_list = MATOM.models.get_chip(tf_class_id)
            else:
                tf_list = []
            if uploaded_chip is True:
                tf_list.append('%s_tmp_chip_%s.posneg' % (tf, str(seq_len)))
            mot = None  # Since we are only using user motifs
            for mot_name in tf_names:
                #print mot_name  #TODO May use this for logging
                user_motif_details = MARSTools.Assess_by_score.get_motif_from_meme(user_motif, mot_name)
                run_assess(score_option, output, user_motif_details, tf_list, mot, tf)
        else:
            if uploaded_chip is True:
                tf_list = ['%s_tmp_chip_%s.posneg' % (tf, str(seq_len))]
                tf_names = get_tfnames(user_motif)
                mot = None  # Since we are only using user motifs
                for mot_name in tf_names:
                #print mot_name  #TODO May use this for logging
                    user_motif_details = MARSTools.Assess_by_score.get_motif_from_meme(user_motif, mot_name)
                    run_assess(score_option, output, user_motif_details, tf_list, mot, tf, seq_len)

score_extensions = {"gomeroccupancyscore": 'gomer', "energyscore": 'energy', "amaoccupancyscore": 'ama',
                    "maxoccupancyscore": 'maxoc', "sumlogoddsscore": 'sumlog', "maxlogoddsscore": 'maxlog',
                    "sumoccupancyscore": 'sumoc'}


def get_tfnames(user_motif):
    tf_names = []
    with open(user_motif) as motif_input:
        for line in motif_input:
            if line.startswith('MOTIF'):
                mot_name = line.split(" ")[1]
                tf_names.append(mot_name)
    return tf_names

