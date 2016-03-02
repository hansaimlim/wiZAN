#!/usr/bin/python
from Bio import SeqIO
from Bio.SeqRecord import SeqRecord
from Bio.Seq import Seq
import sys
#Truncate FASTA sequence up to 3 segments per file

def OutFasta(frags,outfile):
	o_handle = open(outfile, "w")
	SeqIO.write(frags, o_handle, "fasta")
	o_handle.close()
	return


seq_frags=[]
if len(sys.argv)==5:
	#take one segment
	seqfile = sys.argv[1]
	start = int(sys.argv[2])
	end = int(sys.argv[3])
	out = sys.argv[4]
	
	inseq = SeqIO.read(open(seqfile), "fasta")
	frag = inseq.seq[start:end]
	rec = SeqRecord(frag, 'fragment %i:%i' % (start,end),'','')
	seq_frags.append(rec)
	
	OutFasta(seq_frags,out)
elif len(sys.argv)==7:
	#take two segments
	seqfile = sys.argv[1]
	start1 = int(sys.argv[2])
	end1 = int(sys.argv[3])
	start2 = int(sys.argv[4])
	end2 = int(sys.argv[5])
	out = sys.argv[6]

	inseq = SeqIO.read(open(seqfile), "fasta")
	frag1 = inseq.seq[start1:end1]
	frag2 = inseq.seq[start2:end2]
	fr = (str(frag1), str(frag2))
	frags = ''.join(fr)
	rec = SeqRecord(Seq(frags), 'fragment %i:%i-%i:%i' % (start1,end1,start2,end2),'','')


	seq_frags.append(rec)
	OutFasta(seq_frags,out)
elif len(sys.argv)==9:
	#take three segments
	seqfile = sys.argv[1]
	start1 = int(sys.argv[2])
	end1 = int(sys.argv[3])
	start2 = int(sys.argv[4])
	end2 = int(sys.argv[5])
	start3 = int(sys.argv[6])
	end3 = int(sys.argv[7])
	out = sys.argv[8]
	
	inseq = SeqIO.read(open(seqfile), "fasta")
	frag1 = inseq.seq[start1:end1]
	frag2 = inseq.seq[start2:end2]
	frag3 = inseq.seq[start3:end3]
	fr = (str(frag1), str(frag2), str(frag3))
	print fr
#	frags = ''.join(fr)
#	rec = SeqRecord(Seq(frags), 'fragment %i:%i-%i:%i-%i:%i' % (start1,end1,start2,end2,start3,end3),'','')
#	seq_frags.append(rec)
#	OutFasta(seq_frags,out)
else:
	print "Usage: python <truncate_seq.py> <input fasta> (<start> <end>){1-3}"
	print "You may put up to 3 segments in order"
	exit(1)

