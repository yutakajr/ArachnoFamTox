from .util import warn 
import subprocess
import sys 
from pathlib import Path
from Bio import SearchIO 

def run_hmmscan(fasta,args):
    warn("stdout",f"running hmmscan with evalue {args.eHMM}")
    subprocess.run([
                    "hmmscan",
                    "--cpu", str(args.cpu),
                    "-E" , str(args.eHMM),
                    "--domE", str(args.eHMM),
                    "--domtblout" , Path(args.out, "out.hmmer.domtab"),
                    Path("db", "Arachnida.hmm"),
                    fasta],
                    stdout=subprocess.DEVNULL)
    return str(Path(args.out, "out.hmmer.domtab")) 

def run_rpsblast(fasta,args):
    warn("stdout",f"running RPS-BLAST with evalue {args.ePSSM}")
    subprocess.run([
                    "rpsblast", 
                    "-db", Path("db", "ArachnidaToxinsV3.pssm"), #here
                    "-query", fasta, 
                    "-out",  Path(args.out, "out.pssm"),
                    "-evalue", str(args.ePSSM),
                    "-outfmt", "6 qseqid sseqid evalue bitscore ppos pident qlen qstart qend",
                    "-num_threads",  str(args.cpu)])
    return str(Path(args.out, "out.pssm")) 

def run_blastp(fasta,args):
    warn("stdout",f"running BLASTp with evalue {args.eBLASTP}")
    subprocess.run([
                    "blastp",
                    "-db", Path("db", "toxprot"),
                    "-query", fasta, 
                    "-out",  Path(args.out, "out.blastp.toxprot"),
                    "-evalue", str(args.eBLASTP),
                    "-outfmt", "6 qseqid sseqid pident length mismatch gapopen qstart qend sstart send evalue bitscore qcovs ppos", 
                    "-num_threads",  str(args.cpu)])
    return str(Path(args.out, "out.blastp.toxprot")) 

def parse_hmmscan(domtab):
    output = open(domtab + '.parsed', 'w')
    with open(domtab, 'r') as f:
        for record in SearchIO.parse(f, 'hmmscan3-domtab'):
            hits = record.hits
            hsps = record.hsps
            num_hits = len(hits)
            num_hsps = len(hsps)

            if num_hits > 0: 
                for i in range(0,num_hits):
                    if (hsps[i].evalue <= 0.1) or (hsps[i].evalue == 0):
                        output.write(str(record.id) + '\t' + str(hits[i].id) + '\t' + str(hsps[i].evalue)  
                            + '\t' + str(hsps[i].bitscore) + '\t' + str(hsps[i].acc_avg) +  
                            '\t' +  str(record.seq_len) + '\t' + str(hits[i].description) + '\t' + 
                            str(hsps[i].env_start) + '\t' + str(hsps[i].env_end) + '\n')
            else:
                if (hsps.evalue <= 0.1) or (hsps.evalue == 0):
                    output.write(str(record.id) + '\t' + str(hits.id) + '\t' + str(hsps.evalue) + 
                        '\t' + str(hsps.bitscore) + '\t' + str(hsps.acc_avg) + 
                        '\t' +  str(record.seq_len) + '\t' + str(hits.description) + '\t' + 
                        str(hsps.env_start) + '\t' + str(hsps.env_end) + '\n' )
    f.close()
    return domtab + ".parsed"
