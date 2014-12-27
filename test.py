#!/usr/bin/python
# -*- coding: utf8 -*-
#
# yes, the code is horrible. a complete rewrite wouldn't be a bad idea
#
# example usage:
#  mkdir resultsr05/
#  python test.py ./motif -t0 -Mpam30 -T0.5 -i1000 -m50 -r0.5 -cdata/motifs -Oresultsr05/ `ls -Sr data/PS* | head -5 | tr '\n' ' '` > res.txt
#
import sys
import os
import re
import subprocess

def usage():
    print sys.argv[0], "[motif-binary] [motif-options] [-c prosite_motifs] [-O out_prefix] files ..."

# parse arguments
options = ""
files = []
cur_opt = motifs_fn =  None
out_pref = ""
exe = sys.argv[1]
for arg in sys.argv[2:]:
    if cur_opt is not None:
        if cur_opt == '-c':
            motifs_fn = arg
        elif cur_opt == '-O':
            out_pref = arg[2:]
        else:
            if len(arg.split(' ')) > 1:
                arg = '"' + arg + '"'
            options = " ".join([options, cur_opt, arg])
        cur_opt = None
    elif arg[0] == '-':
        if len(arg) <= 2:
            cur_opt = arg
        elif arg[0:2] == '-c':
            motifs_fn = arg[2:]
        elif arg[0:2] == '-O':
            out_pref = arg[2:]
        else:
            if len(arg.split(' ')) > 1:
                arg = arg[0:2] + '"' + arg[2:] + '"'
            options = " ".join([options, arg])
    else:
         # positional argument (file)
        files.append(arg)

#print "./motif", options, "<>", files, "<>", motifs_fn
if len(files) == 0:
    usage()
    exit()
if motifs_fn is None:
    print "[-] Warning! PROSITE motifs not given"

# make sure that the given input files are accessible
for f in files:
    h = open(f, 'r')
    h.close()

# read motifs
motifs = {}
if motifs_fn:
    h = open(motifs_fn, 'r')
    for line in h.readlines():
        line = line.strip()
        m = re.match('([^:]+): (.*)', line)
        if not m.groups() or len(m.groups()) != 2:
            print "[-] Motif file in wrong format!"
            h.close()
            exit()
        motifs[m.group(1)] = m.group(2)
    h.close()

# converts PROSITE motif to corresponding regular expression
def motif_to_regex(motif):
    motif = re.sub(r"{([^}]+)}", r"[^\1]", motif)
    motif = re.sub(r"\(([^\)]+)\)", r"{\1}", motif)
    motif = re.sub(r"x", r".", motif)
    motif = re.sub(r"-", r"", motif)
    return motif

# run_test returns an instance of this class.
# dict would also work...
class test_result:
    def __init__(self):
        self.consensuses = []            # consensus words
        self.consensuses_by_score = []   # ordered by score
        self.motifs = []                 # consensus words containing the correct motif.
        self.motifs_by_score = []
        self.cmotif = None               # correct motif
        self.best_consensus = []

# returns instance of test_result
# TODO: a complete rewrite
def run_test(options, seqfile, psmotif):
    motif = motif_to_regex(psmotif)
    cmd = " ".join([exe, options, seqfile])

    out_file = out_pref + os.path.basename(seqfile) + "_stdout.txt"
    res_file = out_pref + os.path.basename(seqfile) + "_res.txt"
    out = open(out_file, 'w')
    res = open(res_file, 'w')
    out.write(cmd + "\n" + psmotif + "\n\n")
    res.write(cmd + "\n" + psmotif + "\n\n")

    print cmd + ";\nmotif: " + motif
    process = subprocess.Popen(cmd, shell=True, stdout=subprocess.PIPE)

    # process inpout from the program
    found = []
    scores = []
    ret = test_result()
    while True:
        line = process.stdout.readline()
        if not line: break;
        out.write(line)
        line = line.rstrip()

        # iteration line
        m = re.match(r"Iteration (\d+):[^>]+> (.*)$", line)
        if m and len(m.groups()) == 2:
            iteration = m.group(1)
            consensus = m.group(2)
            ret.consensuses.append(consensus)
            
            if len(re.findall(motif, consensus)) > 0:
                ret.motifs.append(consensus)
                found.append((line, consensus))

                if ret.cmotif is None or len(consensus) < len(ret.best_consensus):
                    ret.best_consensus = consensus
                    m = re.search(motif, consensus)
                    ret.cmotif = m.group(0)
                continue

        # sorted by score line
        m = re.match(r"(\d+)\. (\w*)", line)
        if m and len(m.groups()) == 2:
            rank = m.group(1)
            consensus = m.group(2)
            ret.consensuses_by_score.append(consensus)

            if len(re.findall(motif, consensus)) > 0:
                ret.motifs_by_score.append((rank, consensus))
                scores.append((line, consensus))
                continue

    if len(found) == 0:
        print "Correct motif not found in file:", seqfile
        return None
    print "Correct motif found in {}/{} iterations ({})".format(len(ret.motifs), len(ret.consensuses), out_file)
    res.write("Correct motif found in {}/{} iterations ({})\n".format(len(ret.motifs), len(ret.consensuses), out_file))


    for (line, consensus) in found:
        res.write("{}, len: {}\n".format(line, len(consensus)))
    res.write("\n\nSorted by score:\n")
    for (line, consensus) in scores:
        res.write("{}, len: {}\n".format(line, len(consensus)))
     
    out.close()
    res.close()
    return ret

motif = None
matches = dist10 = dist5 = 0
cons_count = cons_len_sum = 0
for f in files:
    # find motif
    if motifs_fn:
        try:
            motif = motifs[os.path.basename(f)]
        except KeyError:
            print "[-] Missing correct motif for sequence set:", f
    print "-" * 70
    ret = run_test(options, f, motif)
    if ret and ret.best_consensus is not None and ret.cmotif is not None:
        print "len: {}/{}".format(len(ret.best_consensus), len(ret.cmotif))
        matches += 1
        d = len(ret.best_consensus) - len(ret.cmotif)
        if d <= 10: dist10 += 1
        if d <= 5: dist5 += 1

        cons_count += len(ret.consensuses)
        for cons in ret.consensuses:
            cons_len_sum += len(cons)

print "-" * 70
print "cons_count: {}, cons_len_sum: {}".format(cons_count, cons_len_sum)
print "matches [%d/%d] %.2f" % (matches, len(files), 1.0*matches/len(files))
print "dist10 [%d/%d] %.2f" % (dist10, len(files), 1.0*dist10/len(files))
print "dist5 [%d/%d] %.2f" % (dist5, len(files), 1.0*dist5/len(files))
