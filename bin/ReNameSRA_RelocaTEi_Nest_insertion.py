#!/opt/Python/2.7.3/bin/python
import sys
from collections import defaultdict
import numpy as np
import re
import os
import argparse
import glob
import time
from Bio import SeqIO
import fnmatch

def usage():
    test="name"
    message='''
python ReNameSRA_RelocaTEi.py --input Japonica_fastq
python ReNameSRA_RelocaTEi_Nest_insertion.py --input Run_folder --repeat ping

Run relocaTE_Nested_insertionFinder.py in all the RelocaTE_i runs in Run_folder for nested Ping/Pong insertions.
1. merge and mapped flanking_seq of Ping/Pong junction reads
2. identify Ping/Pong insertions on Ping/Pong elements
3. clean temporary files  

    '''
    print message


def runjob(script, lines):
    cmd = 'perl /rhome/cjinfeng/BigData/software/bin/qsub-pbs.pl --maxjob 80 --lines %s --interval 120 --resource nodes=1:ppn=1,walltime=10:00:00,mem=10G --convert no %s' %(lines, script)
    #print cmd 
    os.system(cmd)

def fasta_id(fastafile):
    fastaid = defaultdict(str)
    for record in SeqIO.parse(fastafile,"fasta"):
        fastaid[record.id] = 1
    return fastaid


def readtable(infile):
    data = defaultdict(str)
    with open (infile, 'r') as filehd:
        for line in filehd:
            line = line.rstrip()
            if len(line) > 2: 
                unit = re.split(r' ',line)
                data[unit[0]] = line
                #print unit[0], line
    return data

def rerun_shell(outfile, infile, topdir):
    ofile = open(outfile, 'w')
    r = re.compile(r'step_(\d+)')
    with open (infile, 'r') as filehd:
        for line in filehd:
            line = line.rstrip()
            step = 8
            if r.search(line):
                step = int(r.search(line).groups(0)[0])
            if step >= 4:
                print >> ofile, line
    ofile.close()

def main():
    parser = argparse.ArgumentParser()
    parser.add_argument('-i', '--input')
    parser.add_argument('-o', '--output')
    parser.add_argument('-g', '--genome')
    parser.add_argument('-r', '--repeat')
    parser.add_argument('-v', dest='verbose', action='store_true')
    args = parser.parse_args()
    try:
        len(args.input) > 0
    except:
        usage()
        sys.exit(2)

    if not args.output:
        args.output = '%s' %(os.path.abspath(args.input))

    if not args.genome:
        args.genome = '/rhome/cjinfeng/BigData/00.RD/RelocaTE_i/Simulation/Reference/ping.fa'
    if not args.repeat:
        args.repeat = 'ping'

    #-t ../input/mping_UNK.fa -g /rhome/cjinfeng/HEG4_cjinfeng/seqlib/MSU_r7.fa -d ../input/FC52_7 -e HEG4 -o mPing_HEG4_UNK -r 1 -p 1 -a 1   
    #RelocaTE = 'python /rhome/cjinfeng/software/tools/RelocaTE_1.0.3_i/RelocaTE/scripts/relocaTE.py'
    #RelocaTE = 'python /rhome/cjinfeng/BigData/00.RD/RelocaTE2/scripts/relocaTE.py'
    RelocaTE = 'python /rhome/cjinfeng/BigData/00.RD/RelocaTE2_mPing/scripts/relocaTE.py'
    Reference= os.path.abspath(args.genome)
    project = os.path.split(args.output)[1]
    cpu = 16
    if not os.path.exists(project):
        os.mkdir(project)
    print project
    read_dirs = glob.glob('%s/*Ping' %(os.path.abspath(args.input)))
    shell = []
    script= '/rhome/cjinfeng/Rice/Rice_population_sequence/Rice_3000/analysis/Ping_silence/bin'
    repeat= args.repeat
    ofile = open('%s.run.sh' %(args.output), 'w')
    for read_dir in sorted(read_dirs):
        outdir = '%s/%s' %(os.path.abspath(args.output), os.path.split(read_dir)[1])
        #existingTE  = '%s.mPing.RepeatMasker.out' %(Reference)
        # relocate will not run if there is result exists
        #if not os.path.exists(outdir):
        if fnmatch.fnmatch(read_dir, '*Ping'):
            shell.append('cat %s/repeat/flanking_seq/*_1.te_repeat.flankingReads.fq.matched > %s/repeat/flanking_seq/sample_1.te_repeat.flankingReads.fq.matched' %(outdir, outdir))
            shell.append('cat %s/repeat/flanking_seq/*_2.te_repeat.flankingReads.fq.matched > %s/repeat/flanking_seq/sample_2.te_repeat.flankingReads.fq.matched' %(outdir, outdir))
            shell.append('cat %s/repeat/flanking_seq/*.flankingReads.unPaired.info > %s/repeat/flanking_seq/ALL.flankingReads.unPaired.info' %(outdir, outdir))
            read1 = '%s/repeat/flanking_seq/sample_1.te_repeat.flankingReads.fq.matched' %(outdir)
            read2 = '%s/repeat/flanking_seq/sample_2.te_repeat.flankingReads.fq.matched' %(outdir)
            shell.append('/rhome/cjinfeng/BigData/00.RD/RelocaTE2/tools/bwa-0.6.2/bwa aln -l 20 %s %s > %s.sai' %(Reference, read1, read1))
            shell.append('/rhome/cjinfeng/BigData/00.RD/RelocaTE2/tools/bwa-0.6.2/bwa aln -l 20 %s %s > %s.sai' %(Reference, read2, read2))
            prefix = '%s/repeat/bwa_aln/sample' %(outdir)
            shell.append('/rhome/cjinfeng/BigData/00.RD/RelocaTE2/tools/bwa-0.6.2/bwa sampe %s %s.sai %s.sai %s %s > %s.sam' %(Reference, read1, read2, read1, read2, prefix))
            shell.append('/opt/tyler/bin/samtools view -bS -o %s.raw.bam %s.sam' %(prefix, prefix))
            shell.append('/opt/tyler/bin/samtools sort %s.raw.bam %s.sort' %(prefix, prefix))
            shell.append('java -jar /opt/picard/1.81/MarkDuplicates.jar ASSUME_SORTED=TRUE REMOVE_DUPLICATES=TRUE VALIDATION_STRINGENCY=LENIENT INPUT=%s.sort.bam OUTPUT=%s.bam METRICS_FILE=%s.dupli' %(prefix, prefix, prefix))
            shell.append('/opt/tyler/bin/samtools index %s.bam' %(prefix))
            all_info = '%s/repeat/flanking_seq/ALL.flankingReads.unPaired.info' %(outdir)
            shell.append('rm %s %s %s' %(read1, read2, all_info))
            shell.append('rm %s.sam %s.raw.bam %s.sort.bam %s.sai %s.sai' %(prefix, prefix, prefix, read1, read2))
            shell.append('python %s/relocaTE_Nested_insertionFinder.py %s/repeat/bwa_aln/sample.bam %s MSU_r7.fa repeat %s/regex.txt HEG4 100 MSU_r7.fa.mPing.RepeatMasker.out 3 0 500' %(script, outdir, repeat, outdir))
            
            #relocaTE = '%s --te_fasta %s --genome_fasta %s --fq_dir %s --outdir %s --reference_ins %s' %(RelocaTE, Repeat, Reference, read_dir, outdir, existingTE)
            #os.system('cp /rhome/cjinfeng/Rice/Rice_population_sequence/Rice_3000/CAAS/existingTE.bed %s/repeat/' %(outdir))
            #os.system('rm -R %s/repeat/results' %(outdir))
            #os.system('rm %s/repeat/bwa_aln/MSU*' %(outdir))
            #rerun_shell('%s/run_these_jobs_rerun.sh' %(outdir), '%s/run_these_jobs.sh' %(outdir), outdir)
            #shell    = 'bash %s/run_these_jobs_rerun.sh > %s/run.log 2> %s/run.log2' %(outdir, outdir, outdir)
            #os.system(relocaTE)
            #print >> ofile, relocaTE
            print >> ofile, '\n'.join(shell)
    ofile.close()
    runjob('%s.run.sh' %(args.output), 13)
 
if __name__ == '__main__':
    main()

