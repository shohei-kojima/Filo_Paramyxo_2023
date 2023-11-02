#!/usr/bin/env python

"""
Author: Shohei Kojima @ RIKEN
Python3.7
"""

import os,sys,glob,gzip
import numpy as np

min_length = 17_000
max_length = 21_000
promoter_length = 200



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

def complement_rna(seq):
    return seq.translate(str.maketrans('AUGCaugc', 'UACGuacg'))[::-1]

def bulk_fasta_complement(fa):
    for h in fa:
        fa[h]=complement(fa[h])
    return fa

def test_6(seq):
    if len(seq) % 6 == 0:
        return 'True'
    return 'False'


def find_UNNNNN(seq):
    seqlen = len(seq)
    appended = set()
    matched_pos = []
    for i in range(0, promoter_length - 18 + 1):
        if i in appended:
            continue
        if seq[i] == 'A':
            for j in range(1, int(promoter_length / 6) + 1):
                index = i + (6 * j)
                if seq[index] == 'A':
                    continue
                elif j >= 3:
                    matched_pos.append([i, j])
                    for k in range(j):
                        appended.add(i + (6 * k))
                    break
                else:
                    break
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




f='Filoviridae_230726_nucl.fa'
fa=parse_fasta(f)
print(len(fa))
exit()

# pick up motif pos
out=[]
tmp = ['fasta_header', 'seq_length', 'is_rule_of_6', 'strand', 'motif_start_0based', 'motif_end_1based', 'number_of_U', 'seq', 'seq_separated']
out.append('\t'.join([ str(i) for i in tmp ]) + '\n')
for h in fa:
    if not 'complete genome' in h:
        continue
    seq = fa[h].replace('T', 'U')
    is6 = test_6(seq)
    matched_pos = find_UNNNNN(seq)
    for start, rep in matched_pos:
        end = start + (rep * 6)
        motif_seq = seq[start:end]
        motif_seq_sep = ':'.join([ motif_seq[i*6:(i+1)*6] for i in range(rep) ])
        tmp = [h, len(seq), is6, 'genome', start, end, rep, complement_rna(motif_seq), complement_rna(motif_seq_sep)]
        out.append('\t'.join([ str(i) for i in tmp ]) + '\n')
    seq = complement(fa[h]).replace('T', 'U')
    matched_pos = find_UNNNNN(seq)
    for start, rep in matched_pos:
        end = start + (rep * 6)
        motif_seq = seq[start:end]
        motif_seq_sep = ':'.join([ motif_seq[i*6:(i+1)*6] for i in range(rep) ])
        tmp = [h, len(seq), is6, 'antigenome', start, end, rep, complement_rna(motif_seq), complement_rna(motif_seq_sep)]
        out.append('\t'.join([ str(i) for i in tmp ]) + '\n')


with open('result_Filoviridae.230802.tsv', 'w') as outfile:
    outfile.write(''.join(out))

