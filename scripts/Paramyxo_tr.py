#!/usr/bin/env python

"""
Author: Shohei Kojima @ RIKEN
Python3.7
"""

import os,sys,glob,gzip,string
import numpy as np

min_length = 10_000
max_length = 30_000
promoter_search_len = 40
leader_length = 55
motif_search_len = 120
min_orf_size = 450  # when 3000, 1000 aa
start_codon = 'ATG'
term_codons = {'TAG', 'TGA', 'TAA'}
allow_mismatch_ratio = 0.1
indel_cost = 1
mismatch_cost = 1  # A->G, G->A, C->T, T->C
mismatch_cost = 1  # A->T, A->C, G->T, G->C, T->A, T->G, C->A, C->G
AG_set={'A', 'G'}


def parse_fasta(path_to_file):
    tmp={}
    seq=[]
    with open(path_to_file) as infile:
        for line in infile:
            if '>' in line and seq:
                joined_seq = ''.join(seq)
                if min_length <= len(joined_seq) < max_length:
                    tmp[header]=joined_seq
                header=line.strip().replace('>', '')
                seq=[]
            elif '>' in line and not seq:
                header=line.strip().replace('>', '')
            else:
                seq.append(line.strip())
        joined_seq = ''.join(seq)
        if min_length <= len(joined_seq) < max_length:
            tmp[header]=joined_seq
    return tmp

def complement(seq):
    return seq.translate(str.maketrans('ATGCatgc', 'TACGtacg'))[::-1]

def bulk_fasta_complement(fa):
    for h in fa:
        fa[h]=complement(fa[h])
    return fa

def detect_promoter_perfect_match(seq):
    seq1 = seq[:promoter_search_len]
    seq2 = seq[-promoter_search_len:]
    seq2 = complement(seq2)
    n=0
    for c1,c2 in zip(seq1, seq2):
        if c1 == c2:
            n += 1
        else:
            break
    if n > 0:
        promoter=seq1[:n]
    else:
        promoter='NA'
    return promoter,n


def detect_promoter(seq):
    seq1 = seq[:promoter_search_len]
    seq2 = seq[-promoter_search_len:]
    seq2 = complement(seq2)
    seq1_len=len(seq1) + 1
    seq2_len=len(seq2) + 1
    dp=np.zeros((seq1_len, seq2_len))
    for x in np.arange(seq1_len):
        dp[x][0]= x * indel_cost
    for y in np.arange(seq2_len):
        dp[0][y]= y * indel_cost
    # calc. edit dist
    for x in np.arange(1, seq1_len):
        for y in np.arange(1, seq2_len):
            if seq1[x-1] == seq2[y-1]:
                c=0
            else:
                c=mismatch_cost
            dp[x][y]=min(dp[x-1][y]+indel_cost, dp[x][y-1]+indel_cost, dp[x-1][y-1]+c)  # in,del,sub
    for x in np.arange(1, max(seq1_len, seq2_len)):
        _min=np.min([np.min(dp[:, x]), np.min(dp[x, :])])
        penalty = max(2, int(np.ceil(x * allow_mismatch_ratio)))
        if _min > penalty:
            promoter1=seq1[:x]
            promoter2=seq2[:x]
            aseq1, aseq2, alignment = traceback(dp, x, seq1, seq2)
            promoter1=aseq1.replace('-', '')
            promoter2=aseq2.replace('-', '')
            if max(len(promoter1), len(promoter2)) <= 5 and alignment.count('|') <= 2:
                return '', '', 0, '', '', ''
            return promoter1, promoter2, max(len(promoter1), len(promoter2)), aseq1, aseq2, alignment
    aseq1, aseq2, alignment = traceback(dp, x, seq1, seq2)
    promoter1=aseq1.replace('-', '')
    promoter2=aseq2.replace('-', '')
    return promoter1, promoter2, max(len(promoter1), len(promoter2)), aseq1, aseq2, alignment


def traceback(dp, x, seq1, seq2):
    pos1=x
    pos2=x
    aseq1=[]
    aseq2=[]
    alignment=[]
    while (pos1 + pos2) > 0:
        if dp[pos1][pos2] == (dp[pos1 -1][pos2] + indel_cost):
            alignment.append(' ')
            aseq1.append(seq1[pos1 -1])
            aseq2.append('-')
            pos1 -= 1
        elif dp[pos1][pos2] == (dp[pos1][pos2 -1] + indel_cost):
            alignment.append(' ')
            aseq1.append('-')
            aseq2.append(seq2[pos2 -1])
            pos2 -= 1
        elif dp[pos1][pos2] == dp[pos1 -1][pos2 -1]:
            alignment.append('|')
            aseq1.append(seq1[pos1-1])
            aseq2.append(seq2[pos2-1])
            pos1 -= 1
            pos2 -= 1
        else:
            alignment.append(' ')
            aseq1.append(seq1[pos1-1])
            aseq2.append(seq2[pos2-1])
            pos1 -= 1
            pos2 -= 1
    trim=0
    for c in alignment:
        if c != '|':
            trim += 1
        else:
            break
    aseq1= ''.join(aseq1[::-1])
    aseq2= ''.join(aseq2[::-1])
    alignment= ''.join(alignment[::-1])
    return aseq1[:-trim], aseq2[:-trim], alignment[:-trim]


def test_ACC(seq):
    if seq[:3] == 'ACC' and seq[-3:] == 'GGT':
        return 'True'
    return 'False'

def detect_N_orf(cseq):
    seq=complement(cseq)
    orf_pos=[]
    for frame in [0, 1, 2]:
        is_orf = False
        for i in range(frame + leader_length, len(seq), 3):
            if seq[i:i+3] == 'ATG':
                if is_orf is False:
                    start = i
                    is_orf = True
            elif seq[i:i+3] in term_codons:
                if is_orf:
                    if i - start > min_orf_size:
                        orf_pos.append((start, i))
                    is_orf = False
    if len(orf_pos) >= 1:
        N_start = 100000000
        for s,e in orf_pos:
            if s < N_start:
                N_start = s
        N_start = str(N_start)
    else:
        N_start = 'NA'
    return orf_pos, N_start


def find_NNNNNC(seq):
    seqlen = len(seq)
    matched_pos = []
    for i in range(5, seqlen - 13 + 1):
        if seq[i] == 'C' and seq[i + 6] == 'C' and seq[i + 12] == 'C':
            matched_pos.append(i - seqlen)
    return matched_pos


def test_6(seq):
    if len(seq) % 6 == 0:
        return 'True'
    else:
        return 'False'


def find_CGNNNN(seq):
    seqlen = len(seq)
    matched_pos = []
    for i in range(0, seqlen - 18 + 1):
        if seq[i:i+2] == 'CG' and seq[i + 6: i + 8] == 'CG' and seq[i + 12: i + 14] == 'CG':
            matched_pos.append(i - seqlen)
    return matched_pos


def motif_formatter(seq, pos_set, poss, back):
    tmp=[]
    for p in poss:
        if p in pos_set:
            if p + 18 - back == 0:
                motif=seq[p - back:]
            else:
                motif=seq[p - back:p + 18 - back]
            tmp2=[]
            for n,c in enumerate(motif):
                if (n > 0) and (n % 6 == 0):
                    tmp2.append(' ')
                tmp2.append(c)
            tmp.append(''.join(tmp2))
        else:
            tmp.append('')
    return '\t'.join(tmp).replace('T', 'U')


def motif_formatter2(seq, _from, _to):
    if _to == 0:
        tmp=seq[_from:]
    else:
        tmp=seq[_from:_to]
    tmp2=[]
    for n,c in enumerate(tmp):
        if (n > 0) and (n % 6 == 0):
            tmp2.append(' ')
        tmp2.append(c)
    return ''.join(tmp2).replace('T', 'U')




f='../220803_1/Paramyxoviridae_220803_nucl.fa'
fa=parse_fasta(f)
fa=bulk_fasta_complement(fa)

# pick up motif pos
c_pos=set()
cg_pos=set()
for h in fa:
    if not 'complete genome' in h:
        continue
    seq = fa[h]
    NNNNNC = find_NNNNNC(seq[-motif_search_len:])
    c_pos |= set(NNNNNC)
    CGNNNN = find_CGNNNN(seq[-motif_search_len:])
    cg_pos |= set(CGNNNN)
c_pos=sorted(list(c_pos))
cg_pos=sorted(list(cg_pos))
c_pos_header = '\t'.join([ 'NNNNNC:%d' % pos for pos in c_pos ])
cg_pos_header = '\t'.join([ 'CGNNNN:%d' % pos for pos in cg_pos ])


out=['fasta_header\tseq_length\tis_rule_of_6\tstarts_with_ACC..GGT\tpromoter_length\tpromoter_leader\tpromoter_trailer\talignment_leader\talignment_trailer\talignment\tnum_ORF_longer_than_150aa\tstart_of_first_ORF(0-based)\t5UTR_length\tle_num_NNNNNC_found\tle_NNNNNC_positions\tle_num_CGNNNN_found\tle_CGNNNN_positions\tseq[72:96]\t%s\t%s\n' % (c_pos_header, cg_pos_header)]
for h in fa:
    if not 'complete genome' in h:
        continue
    seq = fa[h]
    is6 = test_6(seq)
    promoter1, promoter2, promoter_len, aseq1, aseq2, alignment = detect_promoter(seq)
    is_ACC = test_ACC(seq)
    orf_pos, N_start = detect_N_orf(seq)
    NNNNNC = find_NNNNNC(seq[-motif_search_len:])
    CGNNNN = find_CGNNNN(seq[-motif_search_len:])
    if N_start == 'NA':
        UTR_length = 'NA'
    else:
        UTR_length = int(N_start)
    out.append('%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\n' % (
        h,
        len(seq),
        is6,
        is_ACC,
        promoter_len,
        promoter1,
        promoter2,
        aseq1,
        aseq2,
        alignment,
        len(orf_pos),
        N_start,
        UTR_length,
        len(NNNNNC),
        ';'.join([ str(abs(pos)) for pos in NNNNNC ]),
        len(CGNNNN),
        ';'.join([ str(abs(pos)) for pos in CGNNNN ]),
        motif_formatter2(seq, -96, -72),
        motif_formatter(seq, NNNNNC, c_pos, 5),
        motif_formatter(seq, CGNNNN, cg_pos, 0),
    ))

with open('result_221226.wide.tsv', 'w') as outfile:
    outfile.write(''.join(out))
