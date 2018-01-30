import contig
import tools
import graph

def set_start(longContigsId, startContigsNumber):
    
    startContigsId = set()
    if len(longContigsId)/2 <= startContigsNumber:
        startContigsId = longContigsId
    else:
        count = 0
        cc = 0
        while count < startContigsNumber:    
            node = list(longContigsId)[cc]
            cc += 1
            if node not in startContigsId:
                startContigsId.add(node)
                startContigsId.add(tools.reverse_id(node))
                count += 1  
        #print startContigsNumber
        #print len(startContigsId)
        assert len(startContigsId) == 2*startContigsNumber
    return startContigsId    
  
def read_parameter(f):

    for line in filter(lambda x: len(x)>0, map(str.strip, f)):
        a = line.split(':')
        if a[0] == 'K':
            K = int(a[1])
        elif a[0] == 'startContigsNumber':
            startContigsNumber = int(a[1])
        elif a[0] == 'sim':
            sim = int(a[1])
        elif a[0] == 'maxDeep':
            maxDeep = int(a[1])
        elif a[0] == 'startLen':
            startLen = int(a[1])
        elif a[0] == 'minMatch':
            minMatch = int(a[1])
        elif a[0] == 'scoreThreshold':
            scoreThreshold = int(a[1])
        elif a[0] == 'supportNumberThreshold':
            supportNumberThreshold = int(a[1])
        elif a[0] == 'coverRatio':
            coverRatio = float(a[1])

    return K, startContigsNumber, sim, maxDeep, startLen, minMatch, scoreThreshold, supportNumberThreshold, coverRatio 

def read_SPAdes_graph(f, startContigsNumber, startLen, maxDeep, K):
    
    contigGraph = {}
    reverseContigGraph = {}
    contigs = {}
    longContigsId = set()
    for line in filter(lambda x: len(x)>0, map(str.strip, f)):
        
        a = line.split(':')
        if len(a) == 1 and not a[0].startswith(">EDGE"):
            contigs[ left_id ].string += a[0]
            continue
        line = line[:-1]
        a = line.split(':')
        left = a[0].split('_')
        left_id = left[1]
        if a[0][-1] == '\'':
            left_id = "R" + left_id
            a[0] = a[0][:-1]
        
        left = a[0].split('_')
        contigGraph[left_id] = []    
        contigs[ left_id ] = contig.Contig(left_id, "", float(left[5]))        
        if int(left[3]) > startLen:
            longContigsId.add(left_id)    
        
        if len(a) == 1:
            continue
        adj = a[1].split(',')
        for right in adj:
            right_id = right.split('_')[1]
            if right[-1] == '\'' or right[-2] == '\'':
                right_id = "R" + right_id 
            contigGraph[left_id].append(right_id)
            if right_id not in reverseContigGraph:
                reverseContigGraph[ right_id ] = []
            reverseContigGraph[ right_id ].append(left_id)
  
       
    #print "[info] longContigsId number: ", len(longContigsId)/2
    startContigsId = set_start(longContigsId, startContigsNumber)
    print "[info] contigs number: ", len(contigs)
    print "[info] longContigsId number: ", len(longContigsId)/2
    print "[info] choose ", startContigsNumber, " contigs as start"
    print "[info] choose ", startContigsId, " as start"

     
    for Id in contigs:
        assert tools.reverse_id(Id) in contigs
    
    fout = open("initial_contigs.fasta", "w")
    for key in contigs:
        fout.write(">contigs_%s\n" % key)
        fout.write("%s\n" % contigs[key].string)
    edges = {}

    graph.ADJGraph.init_class_var(contigs, K)
    g = graph.ADJGraph(contigGraph, reverseContigGraph, edges, maxDeep, 2*startLen)
    return g, contigs, startContigsId

def read_blasr(fileName, simThreshold, startLen):
    f = open(fileName, "r")
    line = f.readline().strip()
    detail = {}  
    # key: contigId, value: s, e, l, readId, s, e, l, sim, score
    number = {}  # key: contig, value: read number
    length = {}  # key: contig, value: length
    print "[info] read blasr"
    
    bestContigs = {} #key: readId value: (contigId, score) 
    while line:
        word = line.split()
        contigId = word[0].split("_")[1].strip()
        readDirection = int(word[2].strip()) 
        contigDirection = int(word[3].strip())
        readId = int(word[1].split("_")[1].strip())
        sim = float(word[5].strip())         
        if (sim < simThreshold):                 ## 2, high sim
            line = f.readline().strip();
            continue 
        score = float(word[4])
        readStart = int(word[6].strip())
        readEnd = int(word[7].strip())
        readLen = int(word[8].strip())
        contigStart = int(word[9].strip())
        contigEnd = int(word[10].strip())
        contigLen = int(word[11].strip())
        # startLen = 2000
        #debug 2-15  startLen*0.2 not contigLen*0.2
        #(cause contigLen is different)
        if (contigStart > contigLen*0.25 or contigStart >= startLen - 200):  # debug 20x 0.25, 5x 10x both 0.2 
            line = f.readline().strip();
            continue 
        if (contigEnd < contigLen*0.8 or contigEnd <= startLen):
            line = f.readline().strip();
            continue  
        if (contigEnd < contigLen*0.9 and readEnd <= readLen*0.9):  
            line = f.readline().strip();
            continue 
        
        if contigId not in detail:
            detail[contigId]= []

        detail[contigId].append( (contigStart, contigEnd, contigLen, readId, readStart, readEnd, readLen, sim, score) )

        if readId not in bestContigs:
            bestContigs[ readId ] = []
            bestContigs[ readId ].append((contigId, score))
        elif score < bestContigs[ readId ][0][-1]:
            bestContigs[ readId ] = []
            bestContigs[ readId ].append((contigId, score))
        elif score == bestContigs[ readId ][0][-1]:
            bestContigs[ readId ].append((contigId, score))

        if contigId not in number:
            number[contigId] = 0
        number[contigId] += 1

        if contigId not in length:
            length[ contigId ] = contigEnd
        else:
            if contigEnd > length[ contigId ]:   # map longest length
                length[ contigId ] = contigEnd
        line = f.readline().strip()
    f.close()
    return detail, number, length, bestContigs
   

