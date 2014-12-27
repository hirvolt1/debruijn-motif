#!/usr/bin/python
# -*- coding: utf8 -*-
#
# yes, the code is horrible. a complete rewrite wouldn't be a bad idea
#
# example usage:
#  mkdir resultsr05/
#  python test.py ./motif -t0 -Mpam30 -T0.5 -i1000 -m50 -r0.5 -cdata/motifs 
#  -Oresultsr05/ `ls -Sr data/PS* | head -5 | tr '\n' ' '` > res.txt

# Minor changes by Kalle Karhu on fall 2012 (kalle.karhu@aalto.fi Minor changes
# by Kalle Karhu on fall 2012 (kalle.karhu@aalto.fi)
import sys
import os
import re
import subprocess

#smaller result class to cut down the huuuge size of the original structure
#below
class MiniResult:
  def __init__(self): 
    self.top_GM_rank = 0
    self.top_score_rank_prefiltered = 0
    self.motif_name = ""
    self.is_a_miss = False

  def copy_from_test_result(self, test_result):
    self.top_GM_rank = test_result.top_GM_rank
    self.top_score_rank_prefiltered = test_result.top_score_rank_prefiltered
    self.motif_name = test_result.motif_name
    self.is_a_miss = test_result.is_a_miss


class ResultFolder:
  def __init__(self):
    self.results = []
    self.hits_at_rank = []
    self.hits_at_GM_rank = []
    self.score_prefiltered_missed = 0
    self.gm_would_have_missed = 0
    self.hits_total = 0
    #only using positions 1-5 out of 0-10 of this array, better safe though
    for rank in range(0,11):
      self.hits_at_rank.append(0)
      self.hits_at_GM_rank.append(0)
    
  #end of def
  
#end of class

class LineVariables:
  def __init__(self):
    self.iteration = 0
    self.weight = 0
    self.gm = 0
    self.ngm = 0
    self.consensus = "WRONG NUMBER OF LINE VARIABLES"

def getResultDetails(line):
  result = LineVariables()
  m = re.match(r"Iteration (\d+):[^:]+: (\d+)[^:]+: (\d+)[^:]+: (\d+\.\d+)[^>]+> (.*)$", line)
  if m and len(m.groups()) == 5:
    result.iteration = int(m.group(1))
    result.weight = int(m.group(2))
    result.gm = int(m.group(3))
    result.ngm = float(m.group(4))
    result.consensus = m.group(5)
    #print "read line: {} {} {} {} {}".format(result.iteration, result.weight, result.gm, result.ngm, result.consensus)
  return result
#end def

def getResultDetailsFromScoreLine(line):
  result = LineVariables()
  m = re.match(r"(\d+)\.(.*)$", line)
  #if m and len(m.groups()) == 5:
  result.iteration = int(m.group(1))
  #  result.weight = m.group(2)
  #  result.gm = m.group(3)
  #  result.ngm = m.group(4)
  #  result.consensus = m.group(5)
  #print "read line: {}".format(result.iteration) #, result.weight, result.gm, result.ngm, result.consensus
  return result

class RetAndResultFolder():
  def __init__(self):
    self.ret = []
    self.result_folder = ResultFolder()

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
        self.cons_contains_motif = []    # by consensuses
        self.gms = []                     # by consensuses
        self.highest_gm = 0             
        self.weights = []                # by consensuses
        self.matching_motif_lengths = [] # by consensuses
        self.top_GM_rank = 0
        self.top_score_rank_prefiltered = 0
        self.motif_name = ""
        self.is_a_miss = False

#    def print_roc_data_to(writehandle):
#        size = len(self.consensuses)
#        i = 0
#        while(i < size): #if these are all lists and don't work like vectors, this will be horrible
#            writehandle.write(len(self.consensuses[i]
#      
#            i = i + 1
        
        



def usage():
    print sys.argv[0], "[motif-binary] [motif-options] [-c prosite_motifs] [-O out_prefix] files ..."


# converts PROSITE motif to corresponding regular expression

def motif_to_regex(motif):
    motif = re.sub(r"{([^}]+)}", r"[^\1]", motif)
    motif = re.sub(r"\(([^\)]+)\)", r"{\1}", motif)
    motif = re.sub(r"x", r".", motif)
    motif = re.sub(r"-", r"", motif)
    return motif

# returns instance of test_result
# TODO: a complete rewrite


def run_test(options, seqfile, psmotif):
    motif = motif_to_regex(psmotif)
    cmd = " ".join([exe, options, seqfile])
    roc_file = out_pref + os.path.basename(seqfile) + "_for_rocs.txt"
    motif_cw_file = out_pref + os.path.basename(seqfile) + "_motifs_and_consensus.txt"
    out_file = out_pref + os.path.basename(seqfile) + "_stdout.txt"
    res_file = out_pref + os.path.basename(seqfile) + "_res.txt"
    not_in_top5_filename = out_pref + os.path.basename(seqfile) + "_missed.txt"
    out = open(out_file, 'w')
    res = open(res_file, 'w')
    roc = open(roc_file, 'w')
    mcw = open(motif_cw_file, 'w')
    not_in_top5 = open(not_in_top5_filename, 'w')

    out.write(cmd + "\n" + psmotif + "\n\n")
    res.write(cmd + "\n" + psmotif + "\n\n")
    not_in_top5.write(cmd + "\n" + psmotif + "\n\n")

    print cmd + ";\nmotif: " + motif
    process = subprocess.Popen(cmd, shell=True, stdout=subprocess.PIPE)

    # process inpout from the program
    found = []
    scores = []

    ret = test_result()
    motif_found = 0
    while True:
        line = process.stdout.readline()
        if not line: break;
        out.write(line)
        line = line.rstrip()

        # iteration line
        m = re.match(r"Iteration (\d+):[^:]+: (\d+)[^:]+: (\d+)[^:]+: (\d+\.\d+)[^>]+> (.*)$", line)

        if m and len(m.groups()) == 5:
            iteration = m.group(1)
            weight = m.group(2)
            gm = m.group(3)
            ngm = m.group(4)
            consensus = m.group(5)
            if iteration == "1":
                ret.highest_gm = gm
            
            ret.consensuses.append(consensus)
            ret.gms.append(gm)      #sorted by consensuses
            ret.weights.append(weight) #same
            motif_length = 0
            matching_motif = 0
            #roc.write(len(consensus) + ";" + gm + ";" + ret.highest_gm + ";" + weight + ";")
            roc.write("{};{};{};{};{};".format(len(consensus), gm, ret.highest_gm, weight, ngm))
            
            if len(re.findall(motif, consensus)) > 0:
                if motif_found == 0:
                    motif_found = 1
                    matching_motif = 1
                    motif_length = len(motif)
                    mcw.write("m= {} - {} =cw\n".format(motif, consensus))
                    
                #roc.write(matching_motif + ";" + motif_length  + "\n")
                roc.write("{};{}\n".format(matching_motif, motif_length))
                ret.matching_motif_lengths.append(motif_length)
                ret.cons_contains_motif.append(matching_motif)


                ret.motifs.append(consensus)
                found.append((line, consensus))
                
                if ret.cmotif is None or len(consensus) < len(ret.best_consensus):
                    ret.best_consensus = consensus
                    m = re.search(motif, consensus)
                    ret.cmotif = m.group(0)
                continue
            else:
                ret.matching_motif_lengths.append(motif_length)
                ret.cons_contains_motif.append(matching_motif)
                #roc.write(matching_motif + ";" + motif_length  + "\n")
                roc.write("{};{}\n".format(matching_motif, motif_length))

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

    #one result will be gotten per round, 
    #...that is, per calling this function

    for (line, consensus) in found:
        #print "{}".format(line)
        line_variables = getResultDetails(line)
        res.write("{}, len: {}\n".format(line, len(consensus)))
        if ret.top_GM_rank == 0:# or line_variables.iteration < ret.top_GM_rank:
            ret.top_GM_rank = line_variables.iteration
            print "gm rank: |{}|".format(line_variables.iteration)
            
        
    res.write("\n\nSorted by score:\n")
    for (line, consensus) in scores:
        #print "{}".format(line)
        line_variables = getResultDetailsFromScoreLine(line)
        res.write("{}, len: {}\n".format(line, len(consensus)))
        if ret.top_score_rank_prefiltered == 0:# or line_variables.iteration < ret.top_score_rank_prefiltered:
            ret.top_score_rank_prefiltered = line_variables.iteration
            print "score rank: |{}|".format(line_variables.iteration)

     
    out.close()
    res.close()
    roc.close()
    mcw.close()
    
    ret.motif_name = os.path.basename(seqfile) 
    if ret.top_score_rank_prefiltered > 5 \
        or ret.top_score_rank_prefiltered == 0:
      ret.is_a_miss = True
      print "considered to be a miss due to score {}. Type of score: {}".format(ret.top_score_rank_prefiltered, type(ret.top_score_rank_prefiltered))

    return ret






# MAIN :)
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

motif = None
matches = dist10 = dist5 = 0
cons_count = cons_len_sum = 0

result_folder = ResultFolder()

for f in files:
    # find motif
    if motifs_fn:
        try:
            motif = motifs[os.path.basename(f)]
        except KeyError:
            print "[-] Missing correct motif for sequence set:", f
    print "-" * 70
    ret = run_test(options, f, motif) 
    mini_result = MiniResult()
    if not ret == None:
      mini_result.copy_from_test_result(ret)

    result_folder.results.append(mini_result)
    if mini_result.top_GM_rank > 5 or mini_result.top_GM_rank == 0:
      result_folder.gm_would_have_missed += 1
    else:
      print "attempting to increment hits_at_GM_rank at {}, type {}".format(mini_result.top_GM_rank, type(mini_result.top_GM_rank))
      result_folder.hits_at_GM_rank[mini_result.top_GM_rank] += 1
    if mini_result.is_a_miss:
      result_folder.score_prefiltered_missed += 1
    else:
      print "attempting to increment hits_at_rank at {}, type {}".format(mini_result.top_score_rank_prefiltered, type(mini_result.top_score_rank_prefiltered))
      result_folder.hits_at_rank[mini_result.top_score_rank_prefiltered] += 1
      result_folder.hits_total += 1

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

print "-" * 70
print "score missed total of {} motifs".format(result_folder.score_prefiltered_missed)
print "GM alone would have missed {} motifs".format(result_folder.gm_would_have_missed)
print "missed motifs: [motif name] [rank in GM] [rank in score]"
for result in result_folder.results:
  if result.is_a_miss:
    print "{}\t{}\t{}".format(result.motif_name, result.top_GM_rank, result.top_score_rank_prefiltered)
print "nonmissed motifs: [motif name] [rank in GM] [rank in score]"
for result in result_folder.results:
  if not (result.is_a_miss):
    print "{}\t{}\t{}".format(result.motif_name, result.top_GM_rank, result.top_score_rank_prefiltered)

print "counts of found motifs at top 5 ranks by score"

for i in range(1,6):
  print"{}. : {}".format(i, result_folder.hits_at_rank[i])
print "counts of found motifs at top 5 ranks by GM"

for i in range(1,6):
  print"{}. : {}".format(i, result_folder.hits_at_GM_rank[i])
