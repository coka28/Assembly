# python main -- calls java assembler with command line args

import subprocess as sp
from random import randint,random
from math import tan,pi,log
from time import time,sleep
import traceback as tb
import os

def A1(length=100,lower=8,upper=14,coverage=20,err_rate=0,correction=False):
    assemble(shotgun(randseq(length),(lower,upper),coverage,err_rate,err_rate,print_sequence=False),
             print_reads=False,print_merges=True,subchance=err_rate,gapchance=err_rate,threads=4,error_correction_lvl=int(correction))

def A2():
    truth = input('Type true sequence : \t')
    lower = int(input('min read length : \t'))
    upper = int(input('max read length : \t'))
    cov = int(input('coverage :        \t'))
    err = int(input('error rate :      \t'))
    correct = (1 if input('error correction (y/n) : \t').lower().startswith('y') else 0)
    reads = shotgun(truth,(lower,upper),cov,err,err,print_sequence=False)
    assemble(reads,truth=truth,print_reads=False,print_merges=True,subchance=err,gapchance=err,threads=4,error_correction_lvl=correct)

def A3():
    truth = input('Type true sequence : \t')
    reads = input('Type read set, seperate with comma :\n').replace(' ','').split(',')
    assemble(reads,truth=truth,print_reads=False,print_merges=True,subchance=0.005,gapchance=0.005,threads=4,error_correction_lvl=0)

def B():
    threadNr = int(input('How many threads? \t'))
    if threadNr <1: threadNr = 1
    assemble(shotgun(randseq(3500),(10,20),13,0,0,print_sequence=False),
             print_reads=False,print_merges=True,subchance=0,gapchance=0,threads=threadNr,error_correction_lvl=0)

def assemble(
    reads, # ............................ set of reads
    truth=None, # ....................... true sequence
    print_reads=False, # ................. print reads before assembly
    print_merges=True, # ................ print merges during assembly
    print_init_stats=True, # ............ print nr of reads and avrg length of actually used reads
    force_result=True, # ................ force algorithm to find an assembly, even if a good merge is not found
    print_alignment=True, # ............. print needleman-wunsch alignment of truth and assembly
    print_progress_bar=True, # .......... print progress bar for every round
    error_correction_lvl=1, # ........... run error correction on reads before assembly this many times
    correct_all=True, # ................. correct all reads when error correcting
    silent_mode=False, # ................ don't print anything (used for benchmarking error rate)
    subchance=1e-10, # .................. value used to create sub penalty ( = -ln(subchance) )
    gapchance=1e-10, # .................. value used to create gap penalty ( = -ln(gapchance) )
    confidence_threshold=0.5, # ......... internally used value .. >0.5 should be enough for anything
    threads=4, # ........................ nr of threads to run in java process
    line_width=100, # ................... line width of alignment print-out
    return_error_rate=False): # ......... return error rate as float (used for benchmarking)

    if silent_mode:
        # define passive print method for silent mode
        def print_(*args,**kwargs):
            return
    else: print_ = print

    # print reads to file for inter-process-communication
    with open('reads.txt','w+') as read_file:
        read_file.write(','.join(reads))

    # print reads to console
    if print_reads:
        for r in reads:
            print(r)

    # print truth to file for inter-process-communication .. only used internally for error count
    if isinstance(truth,str):
        with open('truth.txt','w+') as truth_file:
            truth_file.write(truth)
        print_('True sequence:',truth,sep='\n')

    error_rate = ""
    all_output = ""

    # create subprocess using java with arguments
    worker = sp.Popen(['javaw','-jar','Assembly.jar',
                       str(subchance),str(gapchance),str(threads),
                       '1' if print_merges else '0',
                       '1' if print_init_stats else '0',
                       '1' if force_result else '0',
                       '1' if print_alignment else '0',
                       '1' if print_progress_bar else '0',
                       str(int(error_correction_lvl)),
                       str(int(correct_all)),
                       str(confidence_threshold),
                       str(line_width)],
                      stdout=sp.PIPE,stderr=sp.PIPE) # create PIPE to current python stdout

    try:
        # while worker is still running
        while worker.poll()==None:
            output = str(worker.stdout.readline(),'utf8')
            all_output += output
            print_((output if not (output in ['~\r\n',':\r\n','.\r\n']) else output.rstrip('\r\n')),end='')
            if output.find('Estimated minimum runtime:')>-1:
                print_('\r\n| 0 %'+' '*88+'100 % |',end='')
            if output.find('Error rate')>-1:
                error_rate = output
    except KeyboardInterrupt:
        worker.kill()
        raise KeyboardInterrupt('aborted assembly')
        return

    # while output is not empty
    while (line:=(str(worker.stdout.readline(),'utf8'))):
        line = line.rstrip('\r\n')
        all_output += line
        print_(line,end=('' if line in [':','.','~'] else '\r\n'))
        if line.find('Error rate')>-1:
            error_rate = line

    # flush last output (necessary?)
    output = str(worker.stdout.read(),'utf8').rstrip('\r\n')
    all_output += output
    print_()
    print_(output.replace(':\r\n',':'))
    if output.find('Error rate')>-1:
        error_rate = output

    # extract error rate from strings
    epos = error_rate.find(':')
    while epos<len(error_rate) and epos!=-1:
        if error_rate[epos].isnumeric(): break
        else: epos += 1

    try:
        os.remove('reads.txt')
        os.remove('truth.txt')
    except: pass

    if return_error_rate:
        try: error_rate = float(error_rate[epos:error_rate.find('%')])
        except: raise Exception('!')
        return error_rate

# random sequence
def randseq(length,alphabet='ATGC'):
    seq = ""
    for i in range(length):
        seq += alphabet[randint(0,len(alphabet)-1)]

    with open('truth.txt','w+') as truth_file:
        truth_file.write(seq)

    return seq
# random variable with gaussian distribution (kind of?), SD depends on arg n
def bigrand(n):
    r = random()-1/2
    r *= 2
    r_abs = abs(r)
    r_sig = -1 if r<0 else 1
    r_abs **= 1/(2**n) 
    r = r_abs
    r /= 2*r_sig
    r *= pi
    return tan(r)

# shotgun function .. 1st arg : true sequence; 2nd arg : tuple of (min,max) length;
#                     3rd arg : coverage of bare read set (without eliminating reads contained in other reads)
#                     4th arg : chance of substitution in a read; 5th arg : -"- of gap in a read
#                     6th arg : print the sequence (if it was created with randseq method)
def shotgun(sequence,length_interval,coverage,subchance=0,gapchance=0,print_sequence=True):
    if not isinstance(coverage,int) or coverage<=0:
        print('coverage needs to be a positive integer. Rounding abs value of coverage.')

    if print_sequence: print(sequence)
    alphabet = list(set(sequence))
    
    coverage = int(round(abs(coverage)))
    middle_length = (length_interval[1]+length_interval[0])/2
    # this is a in a = 9.4 * e**0.7x --> find x
    def findx(middle_length):
        two_standard_deviations = middle_length/2
        x = two_standard_deviations/9.4
        x = log(x)
        x /= 0.7
        return x
    
    reads = []

    for i in range(coverage):
        seq = sequence
        j = 0
        while j<len(seq):
            if random()<subchance:
                seq = seq[:j] + [a for a in alphabet if a!=seq[j]][randint(0,len(alphabet)-2)] + seq[j+1:]
            elif random()<gapchance*(1+subchance):
                if random()<0.5: seq = seq[:j] + seq[j+1:]
                else: seq[:j] + alphabet[randint(0,len(alphabet)-1)] + seq[j+1:]
            j += 1
        k = 0
        while len(seq)>0:
            x = findx(middle_length)
            cut_at = min(max(int(middle_length+bigrand(x)),length_interval[0]),length_interval[1])
            reads.append(seq[:cut_at])
            k += 1
            seq = seq[cut_at:]
            if len(seq) < length_interval[0]:
                reads[-1] += seq
                # remove spike at higher end of length distribution
                for i in range(1,k):
                    x = findx((len(reads[-i])+length_interval[0])/2)
                    cut_at = min(max(int(middle_length+bigrand(x)),length_interval[0]),len(reads[-i]))
                    shift = reads[-i][:-cut_at]
                    reads[-i-1] += shift
                    reads[-i] = reads[-i][-cut_at:]
                break
    return reads
