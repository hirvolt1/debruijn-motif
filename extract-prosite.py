#!/usr/bin/python2
# -*- coding: utf-8 -*-
#
# a python script to extract PROSITE motifs and associated sequences from
# Swiss-Prot protein database.
#
# - tommi hirvola

import argparse
import gzip
import os
import sys
import re
import xml.sax

# <setup>

# specify PROSITE entries that we are interested in.
# the function receives PROSITE entry name, type and accession number.
# the function should return true on entries we want to extract.
excluded_ACs = set([13, 15, 16, 17, 29, 38, 40, 43, 44, 107])
def filter_prosite_entry(name, etype, accession_number):
    if etype != "PATTERN": # no rules or matrices
        return False
    
    # The first data set described in:
    # An Approximate de Bruijn Graph Approach to
    # Multiple Local Alignment and Motif Discovery
    # in Protein Sequences
    # R. Patwardhan et al.
    AC = int(accession_number[2:]) # PS12345 -> 12345
    if AC < 10 or AC > 119:
        return False
    if AC in excluded_ACs:
        return False
    return True

# dr[0] = ac_nb, dr[1] = entry_name, dr[2] = character_flag
#
# character flags:
# - T  For a true positive.
# - P  For a 'potential' hit; a sequence that belongs to the set under
#      consideration, but which was not picked up because the region(s) that are
#      used as a 'fingerprint' (pattern or profile) is not yet available in the
#      database (partial sequence).
# - N  For a false negative; a sequence which belongs to the set under
#      consideration, but which has not been picked up by the pattern or
#      profile.
# - ?  For an unknown; a sequence which possibly could belong to the set under
#      consideration.
# - F  For a false positive; a sequence which does not belong to the set in
#      consideration.
def filter_database_reference(ac, name, flag):
    return flag == 'T' # only true positives

# </setup>

# parse arguments
parser = argparse.ArgumentParser()
parser.add_argument('-p', '--patterns', default="prosite.dat",
        help="PROSITE pattern data file (default: prosite.dat)")
parser.add_argument('-d', '--database', default="uniprot_sprot.xml.gz",
        help="UniProtKB/Swiss-Prot protein sequence database in XML format"
             " (default: uniprot_sprot.xml.gz)")
parser.add_argument('-o', '--output-dir', default="",
        help="Output directory (default: current directory)")
args = parser.parse_args()

# open PROSITE pattern file
fpat = None
try:
    fpat = open(args.patterns)
except IOError as e:
    print e
    print "[-] Can't open file '{0}'".format(args.patterns)
    exit()

# open the protein database
fdb = None
try:
    if args.database[-3:] == '.gz':
        fdb = gzip.open(args.database, 'rb')
    elif args.database[-4:] == '.xml':
        fdb = open(args.database, 'r')
    else:
        print "[-] Was expecting Swiss-Prot database in .gz or .xml format"
        fpat.close()
        parser.print_help()
        exit()
except IOError as e:
    print e
    print "[-] Can't open file '{0}'".format(args.database)
    fpat.close()
    exit()

# prepare output directory
odir = ""
if args.output_dir and len(args.output_dir) > 0:
    odir = os.path.normpath(args.output_dir)
    if not os.path.exists(odir):
        print "[+] Creating output directory: {0}".format(odir)
        os.makedirs(odir)
    odir = odir + "/"

# progress bar class
class ProgressBar:
    cur = 0
    max = 0
    prev_progress = None
    width = None

    def set_max(self, max): self.max = max
    def set_cur(self, cur): self.cur = cur

    def __init__(self, cur=0, max=100):
        self.cur = cur
        self.max = max
        try:
            self.width = int(os.popen('stty size', 'r').read().split()[1]) - 7
        except:
            self.width = 20
    
    def set_bar_width(self, w): self.width = w
    def print_bar(self):
        progress = (100 * self.cur) / self.max
        if progress != self.prev_progress:
            w = (progress * self.width)/100
            sys.stdout.write('\r[{0}] {1}%'
                    .format('#' * w + ' ' * (self.width - w), progress))
            sys.stdout.flush()
            self.prev_progress = progress

    def clear(self): sys.stdout.write('\r' + ' ' * (self.width + 7) + '\r')

# PROSITE entry class
class PSEntry:
    def parse_ID(self, line): # entry name and type
        m = re.match('ID\s+([^;]*);\s+([^.]*)', line)
        self.name = m.group(1)
        self.etype = m.group(2)
    def parse_AC(self, line): # accession number
        m = re.match('AC\s+([^;]+)', line)
        self.ac = m.group(1)
    def parse_PA(self, line): # pattern
        m = re.match('PA\s+([^;.]+)', line)
        self.pattern += m.group(1)
    def parse_DR(self, line): # database reference
        m = re.match('DR\s+(.*)$', line)
        line = m.group(1)
        m = re.findall('\s*([^,\s]+)\s*,\s*([^,\s]+)\s*,\s*(.)\s*;', line)
        m = tuple(filter(
                lambda x: filter_database_reference(x[0], x[1], x[2]), m))
        self.drs = self.drs.union(m)
        self.sequences += len(m)

    def create_file(self, prefix=None):
        self.f = open(prefix + self.ac, "w")
        return self.f
    def close_file(self):
        self.f.close()
        self.f = None

    def write_seq(self, seq):
        self.f.write(seq if seq[-1] == '\n' else seq + '\n')
        self.written += 1

    def __init__(self):
        self.name = self.etype = self.ac = self.f = None
        self.pattern = ""
        self.sequences = self.written = 0
        self.drs = set()
    def __del__(self):
        if self.f:
            self.close_file()

# parse PROSITE data lines
print "[+] Parsing {0} ...".format(args.patterns)
psentries = []
entry = PSEntry()
handles = {} # maps prot db reference ACs to lists of associated ps entries
progress = ProgressBar(0, os.path.getsize(args.patterns))
for line in fpat:
    ltype = line[0:2]
    if ltype == 'ID':
        entry.parse_ID(line)
    elif ltype == 'AC':
        entry.parse_AC(line)
    elif ltype == 'PA':
        entry.parse_PA(line)
    elif ltype == 'DR':
        entry.parse_DR(line)
    elif ltype == '//': # entry terminator
        if filter_prosite_entry(entry.name, entry.etype, entry.ac):
            entry.create_file(odir)
            for dr in entry.drs:
                dr_ac = dr[0]
                if dr_ac in handles:
                    hlist = handles[dr_ac]
                    hlist.append(entry)
                    handles[dr_ac] = hlist
                else:
                    handles[dr_ac] = [entry]
            psentries.append(entry)
        entry = PSEntry()

        # print progress
        progress.set_cur(fpat.tell())
        progress.print_bar()
progress.clear()

print "[+] {0} PROSITE entries selected".format(len(psentries))

# write motifs
print "[+] Writing motifs to file {0}".format(odir + "motifs")
try:
    f = open(odir + "motifs", "w")
    for pse in psentries:
        f.write("{0}: {1}".format(pse.ac, pse.pattern))
        if pse.pattern[-1] != '\n':
            f.write("\n")
    f.close()
except IOError as e:
    print e

# SAX parser for XML Swiss-Prot protein database
progress = ProgressBar(0, len(handles))
class SP_xml_parser(xml.sax.ContentHandler):
    is_ac = False
    is_seq = False

    cur_ac = ""
    cur_seq = ""

    extract = False
    this_acs = []
    extracted = 0   # for progress bar
    k = 0

    def startElement(self, name, attrs):
        if name == "accession":
            self.is_ac  = True
        elif name == "sequence" and "length" in attrs:
            self.is_seq = True

    # the actual content may be splitted into two or more calls, so we have to
    # concatenate the contents passed to this function and handle the whole
    # content in endElement.
    def characters(self, content):
        if self.is_ac:
            self.cur_ac += content
        elif self.is_seq:
            self.cur_seq += content

    def endElement(self, name):
        if name == "accession":
            # handle <accession>content</accession>.
            # entry may have more than one accession number associated with it.
            self.this_acs.append(self.cur_ac)
            self.cur_ac = ""
            self.is_ac = False
        elif self.is_seq and name == "sequence":
            for dr_ac in self.this_acs:
                if dr_ac not in handles:
                    continue

                # loop through PS entries that reference this sequence
                for pse in handles[dr_ac]:
                    pse.write_seq(self.cur_seq)

                # important! the prot database contains duplicate entries with
                # the same AC and sequence. do not output them multiple times
                handles.pop(dr_ac)

                # increment extracted count for progress bar
                self.extracted += 1
                progress.set_cur(self.extracted)
                progress.print_bar()

                # stops if all sequences were extracted
                if len(handles) == 0:
                    raise xml.sax.SAXException('Done')
            self.cur_seq = ""
            self.is_seq = False
        elif name == "entry":
            self.this_acs = []

# parse
print "[+] Extracting {0} protein sequences ...".format(len(handles))
progress.print_bar()
parser = xml.sax.make_parser()
parser.setContentHandler(SP_xml_parser())
try:
    parser.parse(fdb)
    progress.clear()
    print "[+] Done"
except xml.sax.SAXException as e:
    progress.clear()
    print "[+] Parsing stopped:", str(e)

# check results
for pse in psentries:
    if pse.sequences != pse.written:
        print "[-] Expected {0} sequences for {1} but {2} found".format(pse.sequences, pse.ac, pse.written)
if len(handles) > 0:
    print "[-] {0} sequences missing in total".format(len(handles))
    for key, entries in handles.items():
        print "==>", key
