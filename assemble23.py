import sys
import os
from multiprocessing import Pool
import contig
import graph
import path
import read

def write_final_paths(filename, paths):

    fout = open(filename, "w")
    for key in paths:
        print "debug:", key

        string = paths[key].get_string()
        fout.write(">path_%s_%s\n" % (key, len(string)))
        fout.write("%s\n" % string)



def get_local(Id, g, read3FileName, minMatch, startLen, paths1, scoreThreshold, supportNumberThreshold):

    localGraphId = 0
    firstId = Id 
    paths1[firstId] = []
    checkLoopPath = []
    length = 0
    while 1:
        if localGraphId == 100:
            print "too many graph"
            break

        print "local Id", localGraphId
        paths = g.get_paths_start_one_node(Id)

        #print "debug paths:", paths
        fileName = firstId + "_" + str(localGraphId)
        print fileName
        path.write_paths(fileName, paths)
        run_blasr(fileName, read3FileName, minMatch)
        pathId, mapLength = find_best_path(fileName, sim, startLen, scoreThreshold, supportNumberThreshold)
        if pathId == -1 and localGraphId == 0:
            paths1[firstId].append(firstId)
            break
        elif pathId == -1:
            break
        #print "debug paths:", paths
        print "debug choose pathId:", pathId
        mapPath, mapString = paths[pathId].get_map_path(mapLength)
        if len(mapPath) <= 1:
            break
        if localGraphId == 0:
            paths1[firstId].extend(mapPath)
            length = len(mapString) + path.Path.Contigs[firstId].get_len() - startLen 
        else:
            paths1[firstId].extend(mapPath[1:])
            length = length + len(mapString) - startLen
        
        
        print "debug"
        print localGraphId, "final path", mapPath, "len: ", len(mapString)
        print firstId, "final path", "len: ", length

        Id = mapPath[-1]
        if Id == firstId or mapPath in checkLoopPath:
            print "Case 3:", Id, " have loop"
            break
        checkLoopPath.append(mapPath) 
        if Id not in g.graph or len(g.graph[Id])==0:
            print "Case 4:", Id, " have no output"
            break
        newId = Id + "_l"
        g.addNode(Id, newId)
        path.Path.Contigs[newId] = contig.Contig(newId, mapString[-startLen:], path.Path.Contigs[Id].copyNum) 
        localGraphId += 1
        Id = newId
    print firstId, "final path", paths1[firstId]
    p = path.Path(paths1[firstId][0], paths1[firstId])
    p.write_seq()

def run_blasr(FileName, read3FileName, minMatch):

    pathFileName = FileName + "_path.fasta"
    blasrFileName = FileName + "_path.blasr"
    inputFile = pathFileName + " " + read3FileName 
    runPara = " -m 1 -bestn 20 -nproc 12 -minMatch " + str(minMatch) 
    outputFile = " -out "  + blasrFileName 
    s = "/home/liyanbo/ARCS23_v2/src/blasr " + inputFile + runPara + outputFile
    #os.system(s)

def check_repeat(read2bestContig, length, detail, startLen, scoreThreshold, supportNumberThreshold):
  
    diffPath = []
    supportNumber = {}
    pathScore = {}
    pathId = set()
    for read in read2bestContig:
        print read, " : ", read2bestContig[read]
        onePath = set()
        for (c, s) in read2bestContig[read]:
            onePath.add(c)
            if c in supportNumber:
                supportNumber[c] += 1
            else:
                supportNumber[c] = 1
            if c in pathScore:
                pathScore[c] = min(pathScore[c], s)
            else:
                pathScore[c] = s
        flag = 1
        for d in diffPath:
            if len(d.intersection(onePath)) != 0:
                flag = 0
                break
        if flag:
            diffPath.append(onePath)      
            for c in onePath:
                pathId.add(c)
    print len(diffPath), " different paths"
    print pathId
    print diffPath
    
    sortedScore = sorted(pathScore.items(), lambda x,y: cmp(x[1], y[1]))
    
    # (1) support number not enough
    print "supportNumber:", supportNumber[ sortedScore[0][0] ] 
    if supportNumber[ sortedScore[0][0] ] <= supportNumberThreshold: # best support
        print "[info] Case: support Number small"
        return 1, -1

    # (1) only one path have best read support,  path's length small
    if len(diffPath) == 1 or length[sortedScore[0][0]] < startLen + 300: # 3100 3080, because short, no distinct
        return 0, sortedScore[0][0]
    
    supportRead = {}
    supportScore = {}
    for d in pathId:
        supportRead[d] = set()
        for (s, e, l, readId, s, e, l, sim, score) in detail[d]:
            print d, detail[d]
            supportRead[d].add(readId)
            supportScore[(d,readId)] = score
    '''
    for key in supportRead:
        print key, supportRead[key]

    for key in supportScore:
        print key, supportScore[key]
    '''

    highestScore = sortedScore[0][1]
    highScoreId = []
    for (Id, s) in sortedScore:
        if Id in pathId:
            if s > highestScore + scoreThreshold:
                break
            highScoreId.append(Id)
    
    print "high score Id:", highScoreId
    for Id1 in highScoreId:
        for Id2 in highScoreId:
            inSamePath = 0
            if Id1 != Id2:
                for onePath in diffPath:
                    if Id1 in onePath and Id2 in onePath:
                        inSamePath = 1
                        break
                if inSamePath == 0:
                    #print Id1, Id2
                    sameReadSupport = supportRead[Id1].intersection(supportRead[Id2])
                    if len(sameReadSupport) == 0:
                        print "aa"
                        return 1, -1
                    isDiffScore = 0
                    for r in sameReadSupport:
                        #print r
                        if supportScore[(Id1,r)] != supportScore[(Id2,r)]:
                            isDiffScore = 1
                            break
                    if isDiffScore == 0:
                        #print "cc"    
                        return 1, -1

    return 0, sortedScore[0][0]

    '''
    diffPathMap = {} # key: Id, value: pathId
    diffPathSupportNumber = {} # key: Id, value: support read 
    diffPathLen = {}
    diffPathScore = {} 

    Id = 0
    for d in diffPath: 
        diffPathScore[Id] = 0
        for dd in d:
            if pathScore[dd] < diffPathScore[Id]: 
                diffPathScore[Id] = pathScore[dd] 
                diffPathMap[Id] = dd         
                diffPathLen[Id] = length[dd]
                diffPathSupportNumber[Id] = supportNumber[dd]
        Id += 1        
    '''
    '''
    Id = 0
    for d in diffPath:
        diffPathMap[Id] = d        
        representPath = list(d)[0]
        diffPathLen[Id] = length[representPath]
        diffPathSupportRead[Id] = set()
        diffPathScore[Id] = 0
        for read in read2bestContig:
            for (c, s) in read2bestContig[read]:
                if c == representPath:
                    diffPathSupportRead[Id].add(read)
                    diffPathScore[Id] = min(diffPathScore[Id], s)
        Id += 1
    '''    

    '''
    sortedDiffLen = sorted(diffPathLen.items(), lambda x,y: cmp(y[1], x[1]))
    print ("[debug] different paths' length", sortedDiffLen)
    
    l = sortedDiffLen[0][1]
    sameLenPath = []
    for (Id, ll) in sortedDiffLen:
        if ll >= l-5:
            sameLenPath.append(Id)
        else:
            break
    print "same len number: ", len(sameLenPath)
    
    print diffPathScore
    sortedDiffScore = sorted(diffPathScore.items(), lambda x,y: cmp(x[1], y[1]))
    print sortedDiffScore
    print sortedDiffScore[0][1]
    print sortedDiffScore[1][1]
    
    
    highScoreId = sortedDiffScore[0][0]
   
    # (1) support number not enough
    if diffPathSupportNumber[ highScoreId ] <= supportNumberThreshold:
        print "[info] Case: support Number small"
        return 1, -1
   
    # (2) only one long path or path's length small
    if len(sameLenPath) == 1 or l < startLen + 300: # 3100 3080, because short, no distinct  
        if sortedDiffLen[0][0] != highScoreId:
            print "notice: longest not highest score"
        return 0, diffPathMap[highScoreId]
    

      #score similar, repeat
    if (abs(sortedDiffScore[0][1]-sortedDiffScore[1][1]) < scoreThreshold):
        # same read support and diff socre not repeat
        return 1, -1
     
    return 0, diffPathMap[highScoreId]
    '''




def find_best_path(fileName, sim, startLen, scoreThreshold, supportNumberThreshold): 
    blasrFileName = fileName + "_path.blasr"
    print "read blasr file"
    detail, number, length, read2bestContig = read.read_blasr(blasrFileName, sim, startLen)
    print "finish read"
    bestPathId = -1
    if len(detail) == 0:
        print "[info] Case 1: no read support"
        return -1, -1
    isRepeat, bestPathId = check_repeat(read2bestContig, length, detail, startLen, scoreThreshold, supportNumberThreshold)
    if isRepeat:
        print "[info] Case 2: maybe long repeat, break"
        return -1, -1
    print bestPathId
    return bestPathId, length[bestPathId]  

 

if __name__ == '__main__':
        
    if len(sys.argv) != 6:
        sys.exit("Usage: %s input_type condened_before/spades_graph ref.fasta(F.tul/*.fasta) read3_file_name parameter_file_name" % sys.argv[0])
    with open(sys.argv[5], "r") as f:   
        K, startContigsNumber, sim, maxDeep, startLen, minMatch, scoreThreshold, supportNumberThreshold, coverRatio = read.read_parameter(f)
    
    graphType = int(sys.argv[1]) 
    read3FileName = sys.argv[4]
     
    if graphType == 0:
        with open(sys.argv[2], "r") as f:
            contigGraph, reverseContigGraph, contigs, long_nodes = ReadGraph(f, startContigsNumber, startLen)
    
    print "[info] parameter: K, startContigsNumber, sim, maxDeep, startLen, minMatch, coverRatio", K, startContigsNumber, sim, maxDeep, startLen, minMatch, coverRatio  
    if graphType == 1:
        with open(sys.argv[2], "r") as f:
            g, contigs, startContigsId = read.read_SPAdes_graph(f, startContigsNumber, startLen, maxDeep, K)
    
    print "StartLen 1: ", startLen
    path.Path.init_class_var(K, startLen, contigs, coverRatio)
    print "StartLen 2: ", path.Path.StartLen
    ''' 
    pool = Pool(16)
    paths1 = {}
    #startContigsId = ["238"]
    
    for Id in startContigsId: 
        if Id not in g.graph or len(g.graph[Id]) == 0:
            contigIdSeq = []
            contigIdSeq.append(Id)
            p = path.Path(Id, contigIdSeq)
            p.write_seq()
            print "[info] ", Id, " do not have output"
            continue 
        #get_local(Id, g, read3FileName, minMatch, startLen, paths1)
        pool.apply_async(get_local, args=(Id, g, read3FileName, minMatch, startLen, paths1, scoreThreshold, supportNumberThreshold)) 
    pool.close()
    pool.join()
    #print paths1
     
    paths = path.read_path()
     
    
    write_final_paths("all.fasta", paths)
    '''
     
    Id = "R450"
    paths1 = {}
    get_local(Id, g, read3FileName, minMatch, startLen, paths1, scoreThreshold, supportNumberThreshold)
    paths = {}
    for key in paths1:
       paths[key] = path.Path(key, paths1[key])
    write_final_paths("R450.fasta", paths)
    
    
