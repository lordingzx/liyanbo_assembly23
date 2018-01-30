import sys
import os
import tools
def read_blasr_2(f):
    line = f.readline().strip()
    
    print "[info] read blasr"
    
    bestLink = "" 
    highScore = 0
    while line:
        word = line.split()
         
        #contigId = word[0].split("_")[1].strip()
        readDirection = int(word[2].strip()) 
        contigDirection = int(word[3].strip())

        #readId = word[1].split("_")[1].strip()
        
        ''' 
        if readDirection != 0:
            line = f.readline().strip()
            continue
        if contigDirection != 0:
            line = f.readline().strip()
            continue
        '''
 
        
        sim = float(word[5].strip())         
        score = float(word[4])
        
        readStart = int(word[6].strip())
        readEnd = int(word[7].strip())
        readLen = int(word[8].strip())
        contigStart = int(word[9].strip())
        contigEnd = int(word[10].strip())
        contigLen = int(word[11].strip())

        
        if contigStart < 2000 and contigEnd > contigLen-2000:
            print line
        
        if score < highScore:
            highScore = score
            bestLink = word[0]
        line = f.readline().strip()

        
    print "score:", bestLink
            
    return bestLink 

def read_blasr(f):
    line = f.readline().strip()
    
    print "[info] read blasr"
    link = {} 
    while line:
        word = line.split()
        contigId = word[0].split("_")[1].strip()
        readDirection = int(word[2].strip()) 
        contigDirection = int(word[3].strip())

        readId = word[1].split("_")[1].strip()

        if readDirection != 0:
            line = f.readline().strip()
            continue
        if contigDirection != 0:
            line = f.readline().strip()
            continue
        if not contigId.endswith("tail"):
            line = f.readline().strip()
            continue
        if not readId.endswith("head"):
            line = f.readline().strip()
            continue
        
            
        #sim = float(word[5].strip())         
        #score = float(word[4])
        

        readStart = int(word[6].strip())
        readEnd = int(word[7].strip())
        readLen = int(word[8].strip())
        contigStart = int(word[9].strip())
        contigEnd = int(word[10].strip())
        contigLen = int(word[11].strip())
        
        if readStart != 0:
            line = f.readline().strip()
            continue
        if contigEnd != contigLen:
            line = f.readline().strip()
            continue
        if contigId not in link:
            link[contigId] = []
        link[contigId].append((readId, readEnd))
        line = f.readline().strip()
    for key in link:
        print key, link[key]
        
    return link 
   

def read_scaffolds(f):
    scaffolds = {}
    count = 0
    for line in filter(lambda x: len(x)>0, map(str.strip, f)):
        a = line.split('|')
        if len(a) == 2:
            count = count + 1
            scaffolds[count] = ""
        else:
            scaffolds[count] = scaffolds[count] + line
    contigs = {}

    for key in scaffolds:
        print key, len(scaffolds[key])
        key_h = str(key) + "head"
        contigs[key_h] = scaffolds[key][0:10000]

        key_t = str(key) + "tail"
        contigs[key_t] = scaffolds[key][-10000:]
        contigs["R"+key_t] = tools.reverse_string( contigs[key_h] )
        contigs["R"+key_h] = tools.reverse_string( contigs[key_t] )

    print len(contigs) 
    return scaffolds, contigs







            
if __name__ == '__main__':
    if len(sys.argv) != 2:
        sys.exit("Usage: %s scaffolds.filename" % sys.argv[0])
       
    
    with open(sys.argv[1], "r") as f:
        scaffolds, contigs = read_scaffolds(f)
    '''  
    fileName = "head_tail.fasta"
    fout = open(fileName, "w")
    for key in contigs:
        print key, len(contigs[key])
        fout.write(">path_%s_%s\n" % (key, len(contigs[key])))
        fout.write("%s\n" % contigs[key])
    
    #runPara = " -m 1 -bestn 10 -nproc 12 -minMatch 15" 
    outputFile = " -out head_tail.blasr" 
    s = "/home/liyanbo/ARCS23_v2/src/blasr " + fileName + " " + fileName + outputFile
    print s
    #os.system(s)
    '''
     
    with open("head_tail.blasr", "r") as f:
        link = read_blasr(f)
    
    #foutRunBlasr = open("run_blasr.sh", "w")
    
    for key in link:
        '''   
        fileName = key+"_ChooseNext.fasta"
        
        fout = open(fileName, "w")
        for (Id, pos) in link[key]: 
            #string = contigs[key][-4000:] + contigs[Id][pos:pos+1000]
           
            string = contigs[key][8000-pos:] + contigs[Id][pos:pos+2000]
            fout.write(">path_%s_%s_%s\n" % (key, Id, len(string)))
            fout.write("%s\n" % string)
          
        
        #runPara = " -m 1 -bestn 10 -nproc 12 -minMatch 15" 
        outputFile = " -out " + key + "_ChooseNext.blasr" 
        s = "/home/liyanbo/ARCS23_v2/src/blasr " + fileName + " " + "re_pacbio_read.fasta" + outputFile
        print s

        foutRunBlasr.write(s) 
        foutRunBlasr.write("\n")
        #os.system(s)
       '''
        outputFileName = key + "_ChooseNext.blasr" 
      
        with open(outputFileName, "r") as f:
            one_link = read_blasr_2(f)
        
        #break
    
    

    
