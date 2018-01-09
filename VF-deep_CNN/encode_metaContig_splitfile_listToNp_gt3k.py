def encodeSeq(seq) : 
    seq_code = list()
    for pos in range(len(seq)) :
        letter = seq[pos]
        if letter in ['A', 'a'] :
            code = [1,0,0,0]
        elif letter in ['C', 'c'] :
            code = [0,1,0,0]
        elif letter in ['G', 'g'] :
            code = [0,0,1,0]
        elif letter in ['T', 't'] :
            code = [0,0,0,1]
        else :
            code = [1/4, 1/4, 1/4, 1/4]
        seq_code.append(code)
    return seq_code 


import os, sys, glob, path
import random
import csv
from Bio.Seq import Seq
import numpy as np
import pickle



    
#reader = csv.reader(open("/staging/fs3/renj/DL/source/hostData.csv", "r"), delimiter=",")
#x = list(reader)
#result = np.array(x)
#print(result[1:])
#faList = [x for x in os.listdir(fileDir) if ".fa" in x ]



#import os, sys, glob, path
#os.environ['THEANO_FLAGS'] = "floatX=float32,openmp=True" 
#os.environ['OMP_NUM_THREADS'] = str(multiprocessing.cpu_count())


import optparse
prog_base = os.path.split(sys.argv[0])[1]

parser = optparse.OptionParser()
parser.add_option("-i", "--fileName", action = "store", type = "string", dest = "fileName",
									help = "fileName")
parser.add_option("-l", "--contigLength", action = "store", type = int, dest = "contigLength",
									help = "contigLength")
(options, args) = parser.parse_args()
if (options.fileName is None or options.contigLength is None ) :
	sys.stderr.write(prog_base + ": ERROR: missing required command-line argument")
	parser.print_help()
	sys.exit(0)


#contigLength = 1000
contigLength = options.contigLength
contigLengthk = contigLength/1000
if contigLengthk.is_integer() :
    contigLengthk = int(contigLengthk)

contigType = "phage"

fileName = options.fileName
NCBIName = os.path.splitext((os.path.basename(fileName)))[0]


#fileDir = "/auto/cmb-12/fs3/renj/metagenome/marine_virus_POV/rawData/assemble/POVcontigs_newbler"
fileDir = os.path.dirname(fileName)
outDir0 = fileDir
outDir = os.path.join(outDir0, "encode")
if not os.path.exists(outDir):
    os.makedirs(outDir)

#codeFileName = contigType+"#"+NCBIName+"#"+str(contigLengthk)+"k_code.npy"
#nameFileName = contigType+"#"+NCBIName+"#"+str(contigLengthk)+"k_name.npy"
#print(codeFileName)

fileCount = 0
with open(fileName, 'r') as faLines :
    code = list()
    name = []
    head = ''
    lineNum = 0
    seqCat = ''
    flag = 0
    for line in faLines :
        #print(line)
        if flag == 0 and line[0] == '>' :
            lineNum += 1
            head = line.strip()
            continue
        elif line[0] != '>' :
            seqCat = seqCat + line.strip()
            flag += 1
            lineNum += 1
        elif flag > 0 and line[0] == '>' :
            lineNum += 1
            print("seqCatLen="+str(len(seqCat)))
            if len(seqCat) < 3000 :
                print("seq < 3kb")
            else :
                pos = 0
                posEnd = pos + contigLength
                #pos = 0
                #posEnd = pos + contigLength
                while posEnd <= len(seqCat) :
                    #print(pos)
                    contigName = contigType+"#"+NCBIName+"#"+str(contigLengthk)+"k#"+head+"#"+str(pos)+"#"+str(posEnd)
                    seq = seqCat[pos:posEnd]

                    countN = seq.count("N")
                    if countN/len(seq) <= 0.3 : 
                        name.append(contigName)
                        seq_code = encodeSeq(seq)
                        code.append(seq_code)
                        seqR = Seq(seq).reverse_complement()
                        seqR_code = encodeSeq(seqR)
                        code.append(seqR_code)
                    else :
                        print("NNN")

                    pos = posEnd
                    posEnd = pos + contigLength

                    if len(name) > 0 and len(name) % 1000 == 0 :
                        print("lineNum="+str(lineNum)+",contigNum="+str(len(name)))
                        fileCount += 1
                        codeFileName = contigType+"#"+NCBIName+"#"+str(contigLengthk)+"k_num"+str(fileCount)+"_seq"+str(int(len(code)/2))+"_code.npy"
                        nameFileName = contigType+"#"+NCBIName+"#"+str(contigLengthk)+"k_num"+str(fileCount)+"_seq"+str(len(name))+"_name.npy"     
                        print(codeFileName)
                        np.save( os.path.join(outDir, codeFileName), np.array(code) )
                        np.save( os.path.join(outDir, nameFileName), np.array(name) )
                        #pickle.dump(code, open( os.path.join(outDir, codeFileName), "wb" ))
                        #pickle.dump(name, open( os.path.join(outDir, nameFileName), "wb" ))

                        code = []
                        name = []

            flag = 0
            seqCat = ''
            head = line.strip()
          
    if flag > 0 :
        lineNum += 1
        print("seqCatLen="+str(len(seqCat)))
        if len(seqCat) < 3000 :
            print("seq < 3kb")
        else :
            pos = 0
            posEnd = pos + contigLength
            #pos = 0
            #posEnd = pos + contigLength
            while posEnd <= len(seqCat) :
                #print(pos)
                contigName = contigType+"#"+NCBIName+"#"+str(contigLengthk)+"k#"+head+"#"+str(pos)+"#"+str(posEnd)
                seq = seqCat[pos:posEnd]

                countN = seq.count("N")
                if countN/len(seq) <= 0.3 : 
                    name.append(contigName)
                    seq_code = encodeSeq(seq)
                    code.append(seq_code)
                    seqR = Seq(seq).reverse_complement()
                    seqR_code = encodeSeq(seqR)
                    code.append(seqR_code)
                else :
                    print("NNN")

                pos = posEnd
                posEnd = pos + contigLength

                if len(name) > 0 and len(name) % 1000 == 0 :
                    print("lineNum="+str(lineNum)+",contigNum="+str(len(name)))
                    fileCount += 1
                    codeFileName = contigType+"#"+NCBIName+"#"+str(contigLengthk)+"k_num"+str(fileCount)+"_seq"+str(int(len(code)/2))+"_code.npy"
                    nameFileName = contigType+"#"+NCBIName+"#"+str(contigLengthk)+"k_num"+str(fileCount)+"_seq"+str(len(name))+"_name.npy"
                    print(codeFileName)
                    np.save( os.path.join(outDir, codeFileName), np.array(code) )
                    np.save( os.path.join(outDir, nameFileName), np.array(name) )
                    #pickle.dump(code, open( os.path.join(outDir, codeFileName), "wb"))
                    #pickle.dump(name, open( os.path.join(outDir, nameFileName), "wb"))

                    code = []
                    name = []

if len(code) > 0 :
    codeFileName = contigType+"#"+NCBIName+"#"+str(contigLengthk)+"k_num"+str(fileCount+1)+"_seq"+str(int(len(code)/2))+"_code.npy"
    nameFileName = contigType+"#"+NCBIName+"#"+str(contigLengthk)+"k_num"+str(fileCount+1)+"_seq"+str(len(name))+"_name.npy"
    print(codeFileName)
    np.save( os.path.join(outDir, codeFileName), np.array(code) )
    np.save( os.path.join(outDir, nameFileName), np.array(name) )
    #pickle.dump(code, open( os.path.join(outDir, codeFileName), "wb" ))
    #pickle.dump(name, open( os.path.join(outDir, nameFileName), "wb" ))

           
