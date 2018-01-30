import os
import tools
import graph

def write_paths(fileName, paths):
    print "StartLen:", Path.StartLen
    minLen = 2 * Path.StartLen
    for p in paths:
        string = p.get_local_path_string()
        minLen = min(len(string), minLen)
   
    print "minLen: ", minLen
    fileName = fileName + "_path.fasta"
    fout = open(fileName, "w")
    for p in paths:
        fout.write(">path_%s_%s\n" % (p.Id, minLen))
        fout.write("%s\n" % p.string[:minLen])

def diff_seq(seq1, seq2):
    if len(seq1) != len(seq2):
        return 10000
    count = 0
    for i in range(len(seq1)):
        if seq1[i] != seq2[i]:
            count += 1
    return count        


'''
def link_Rcontig_contig(paths):
    newPaths = {}
    for Id in paths:
        if Id.startswith("R"):
'''

def read_path():
    s = "cat *final_path > all_contigs_id" 
    os.system(s)
    fileName = "all_contigs_id"
    f = open(fileName, "r")
    
    paths = {}
    for line in filter(lambda x: len(x)>0, map(str.strip, f)):
        a = line.split()
        paths[a[0]] = Path(a[0],a)
    
    print "debug number of paths:", len(paths)
    paths2 = {}    
    for key in paths:
        if key.startswith("R"):
            if key[1:] not in paths:
                print key[1:], "not finish"
            assert key[1:] in paths
        if not key.startswith("R"):
            reverseKey = "R" + key
            paths2[key] = paths[key].merge_reverse(paths[reverseKey])
    print "debug number of paths2:", len(paths2)
    '''
    for key in paths2:
        print key, paths[key].contigIdSeq
    '''
    assert 2*len(paths2) == len(paths)
    paths3 = remove_cover(paths2)

    print "debug number of paths3:", len(paths3)
     
    for key in paths3:
        print key, paths3[key].contigIdSeq
     
    paths4 = link(paths3)  
    
    
    print "debug number of paths4:", len(paths4)
    for key in paths4:
        print key, paths4[key].contigIdSeq
    
    #return paths3
    return paths4 

def consensus_link(linkInfo):
    visited = set()
    linkF = {}
    for (l, r) in linkInfo:
        lR = tools.reverse_id(l) 
        rR = tools.reverse_id(r) 
        if lR in visited:
            linkF[(rR, lR)] = linkInfo[(l,r)]
            visited.add(rR)
        elif rR in visited:
            linkF[(rR, lR)] = linkInfo[(l,r)]
            visited.add(lR)
        else:
            linkF[(l, r)] = linkInfo[(l,r)]
            visited.add(l)
            visited.add(r)
    return linkF



def filter_link(linkinfo):
    visitedl = {}
    visitedr = {}

    filteredge = set()
    for (l, r) in linkinfo:
        if l not in visitedl:
            visitedl[l] = r
        elif r != visitedl[l] and linkinfo[(l, visitedl[l])][1] == linkinfo[(l,r)][1]:
            filteredge.add((l, visitedl[l]))
            filteredge.add((l, r))
        elif r != visitedl[l] and linkinfo[(l, visitedl[l])][1] > linkinfo[(l,r)][1]:
            filteredge.add((l, r))
        elif r != visitedl[l] and linkinfo[(l, visitedl[l])][1] < linkinfo[(l,r)][1]:
            filteredge.add((l, visitedl[l]))
            visitedl[l] = r
      
        if r not in visitedr:
            visitedr[r] = l
        elif l != visitedr[r] and linkinfo[(visitedr[r],r)][1] == linkinfo[(l,r)][1]:
            filteredge.add((visitedr[r],r))
            filteredge.add((l, r))
        elif l != visitedr[r] and linkinfo[(visitedr[r],r)][1] > linkinfo[(l,r)][1]:
            filteredge.add((l, r))
        elif l != visitedr[r] and linkinfo[(visitedr[r],r)][1] < linkinfo[(l,r)][1]:
            filteredge.add((visitedr[r],r))
            visitedr[r] = l
    linkf= {}
    for key in linkinfo:
        if key not in filteredge:
            linkf[key] = linkinfo[key]
    print "after filter:", linkf        
    return linkf



def link(paths):

    linkInfo = {} 
    visitedPair = set()
    for l in paths:
        for r in paths:
            if l == r:
                continue
            if (r, l) in visitedPair:
                continue
            visitedPair.add( (l,r) )
            pre, after, linkLen = paths[l].get_link(paths[r]) 
            if linkLen[0] != 0:
                linkInfo[(pre,after)] = linkLen
    print linkInfo            
    linkF1 = consensus_link(linkInfo)
    linkF2 = filter_link(linkF1)
    
    finalGraph = graph.ADJGraph({},{},linkF2,0,0)
    p = finalGraph.link_one_in_one_out()
    finalPaths = {}
    visitedId = set()
   
    for pp in p:
        
        Id = "merge"
        print Id, pp
        if pp[0] not in paths:
            r = tools.reverse_id(pp[0])
            assert r in paths
            paths[pp[0]] = Path(pp[0], tools.reverse_seq(paths[r].contigIdSeq))
            visitedId.add(r)
        seq = paths[pp[0]].contigIdSeq
        print "0:", paths[pp[0]].contigIdSeq
        Id = Id + pp[0]
        visitedId.add(pp[0])
        for i in range(1,len(pp)):
            assert (pp[i-1], pp[i]) in linkF2

            if pp[i] not in paths:
                r = tools.reverse_id(pp[i])
                assert r in paths
                paths[pp[i]] = Path(pp[i], tools.reverse_seq(paths[r].contigIdSeq))
                visitedId.add(r)
            Id = Id + "_" + pp[i]
            print i, paths[pp[i]].contigIdSeq
            visitedId.add(pp[i])
            '''
            print pp[i]
            print len(paths[pp[i]].contigIdSeq)
            print linkInfo[(pp[i-1], pp[i])]
            '''
            seq.extend(paths[pp[i]].contigIdSeq[linkF2[ (pp[i-1], pp[i]) ][0]:])
            print seq
        finalPaths[Id] = Path(Id, seq)

    for key in paths:
        if key not in visitedId:
            finalPaths[key] = paths[key]
    return finalPaths



def remove_cover(paths):
    print "cover Ratio", Path.CoverRatio
    beCovered = set()
    for l in paths:
        for r in paths:
            if r in beCovered or l in beCovered:
                continue
            if l != r and (paths[r].seqLen < paths[l].seqLen
                            or (paths[r].seqLen == paths[l].seqLen and paths[r].get_len() <= paths[l].get_len())): 
                if paths[l].cover(paths[r]):
                    beCovered.add(r)
    pathsNoBeCovered = {}
    for key in paths:
        if key not in beCovered:
            pathsNoBeCovered[key] = paths[key]
    return pathsNoBeCovered        

    
class Path(object):
    K = 127
    StartLen = 3000
    Contigs = {}
    CoverRatio = 0.05
    @classmethod
    def init_class_var(cls, k, startLen, contigs, coverRatio):
        cls.K = k
        cls.StartLen = startLen
        cls.Contigs = contigs
        cls.CoverRatio = coverRatio

    def __init__(self, Id, contigIdSeq):
        self.Id = Id
        self.contigIdSeq = contigIdSeq
        self.seqLen = len(contigIdSeq)
        self.string = ""
        self.len = 0

    def merge_reverse(self, r):
        assert "R" + self.Id == r.Id
        #self.print_seq()
        #r.print_seq()
        res = tools.reverse_seq(r.contigIdSeq)
        res.extend(self.contigIdSeq[1:]) 
        p = Path(self.Id, res)
        return p

    def diff(self, seq):
        assert self.seqLen == len(seq)
        count = diff_seq(self.contigIdSeq, seq)
        '''
        for i in range(self.seqLen):
            if self.contigIdSeq[i] != seq[i]:
                count += 1     
        '''
        return count


    def cover(self, short):      
        for i in range(0,  self.seqLen - short.seqLen+1): 
            if short.contigIdSeq == self.contigIdSeq[i : short.seqLen + i]:            
                print "cover case 1: ",self.Id, " cover ", short.Id
                return 1
            if tools.reverse_seq(short.contigIdSeq) == self.contigIdSeq[i : short.seqLen + i]:
                print "cover case 2: ",self.Id, " cover ", short.Id
                return 1
            #if short.seqLen > 10 and short.diff(self.contigIdSeq[i : short.seqLen + i]) < 4:             

            if short.diff(self.contigIdSeq[i : short.seqLen + i]) / float(short.seqLen)  < Path.CoverRatio:              
                print "cover case 3: ", self.Id, "cover", short.Id
                print short.diff(self.contigIdSeq[i : short.seqLen + i]) / float(short.seqLen) 
                print self.contigIdSeq
                print short.contigIdSeq
                return 1
            if short.diff(tools.reverse_seq(self.contigIdSeq[i : short.seqLen + i])) / float(short.seqLen) < Path.CoverRatio:
                print "cover case 4: ", self.Id, "cover", short.Id
                print short.diff(tools.reverse_seq(self.contigIdSeq[i : short.seqLen + i])) / float(short.seqLen) 
                print self.contigIdSeq
                print tools.reverse_seq(short.contigIdSeq)
                return 1
        return 0
    
    def get_link(self, other):
        length = max(self.seqLen, other.seqLen)
        
        for i in range(0, length):
            if self.seqLen - i <= other.seqLen and self.seqLen-i > 0:
                
                cc = diff_seq(self.contigIdSeq[i:self.seqLen], other.contigIdSeq[0:self.seqLen-i])
                if cc/float(self.seqLen-i)<0.05: # tail and end
                    ##############################
                    #self    1 2 3 4 5
                    #other         4 5 6 7 8 9
                    #############################
                    print "case1: ", self.Id, " links ", other.Id, " length ", self.seqLen-i, " diff ", cc 
                    path1 = Path("path1", self.contigIdSeq[i:self.seqLen]) 
                    path2 = Path("path2", other.contigIdSeq[0:self.seqLen-i])
                    print path1.get_len()
                    print path2.get_len()
                    len1 = path1.get_len()
                    len2 = path2.get_len()
                    #if abs(len1 - len2) > 0.01*len1 or min(len1, len2) < 1000:
                    if abs(len1 - len2) > 0.01*len1 or length < 2:
                        return -1, -1, (0,0)
                    return self.Id, other.Id, (self.seqLen-i, min(len1, len2))
            
            if other.seqLen - i <= self.seqLen and other.seqLen-i > 0:
                cc = diff_seq(other.contigIdSeq[i:other.seqLen], self.contigIdSeq[0:other.seqLen-i])
                if cc/float(other.seqLen-i)<0.05: # tail and end
                    ##############################
                    #self          4 5 6 7 8 9
                    #other   1 2 3 4 5       
                    #############################
                    print "case2: ", other.Id, " links ", self.Id, " length ", other.seqLen-i, " diff ", cc 
                    path1 = Path("path1", other.contigIdSeq[i:other.seqLen]) 
                    path2 = Path("path2", self.contigIdSeq[0:other.seqLen-i])
                    print path1.get_len()
                    print path2.get_len()
                    len1 = path1.get_len()
                    len2 = path2.get_len()
                    #if abs(len1 - len2) > 0.01*len1 or min(len1, len2) < 1000:
                    if abs(len1 - len2) > 0.01*len1 or length < 2:
                        return -1, -1, (0,0)

                    
                    return other.Id, self.Id, (other.seqLen-i, min(len1, len2))

        for i in range(min(self.seqLen, other.seqLen),0,-1):
            cc = diff_seq(tools.reverse_seq(self.contigIdSeq[self.seqLen-i:self.seqLen]) , other.contigIdSeq[other.seqLen-i:other.seqLen])
            if cc/float(i)<0.05: # both tail  
                    ##############################
                    #self    1  2  3   4  5  6
                    #other   9  8  7  R6 R5 R4
                    #############################
                    print "case3:", self.Id, "links", "R" + other.Id, "length", i, " diff ", cc
                    path1 = Path("path1", self.contigIdSeq[self.seqLen-i:self.seqLen]) 
                    path2 = Path("path2", other.contigIdSeq[other.seqLen-i:other.seqLen])
                    print path1.get_len()
                    print path2.get_len()
                    len1 = path1.get_len()
                    len2 = path2.get_len()
                    #if abs(len1 - len2) > 0.01*len1 or min(len1, len2) < 1000:
                    if abs(len1 - len2) > 0.01*len1 or length < 2:
                        return -1, -1, (0,0)

                    return self.Id, "R"+other.Id, (i, min(len1, len2))

            cc = diff_seq(tools.reverse_seq(self.contigIdSeq[0:i]) , other.contigIdSeq[0:i])
            if cc/float(i)<0.05: # both head
                    ##############################
                    #self    4   5  6   7   8  9
                    #other  R6  R5 R4  R3  R2 R1
                    #############################
                    print "case4:", "R"+self.Id, "links", other.Id, "length", i, " diff ", cc

                    path1 = Path("path1", self.contigIdSeq[0:i]) 
                    path2 = Path("path2", other.contigIdSeq[0:i])
                    print path1.get_len()
                    print path2.get_len()

                    len1 = path1.get_len()
                    len2 = path2.get_len()
                    #if abs(len1 - len2) > 0.01*len1 or min(len1, len2) < 1000:
                    if abs(len1 - len2) > 0.01*len1 or length < 2:
                        return -1, -1, (0,0)
                    return "R"+self.Id, other.Id, (i, min(len1, len2))
        return -1, -1, (0,0)
  
    
    
    def print_seq(self):
        print self.Id, "seq:", self.contigIdSeq
            


    def write_seq(self):
        fileName = self.contigIdSeq[0] + "_final_path"
        fout = open(fileName, "w")
        for Id in self.contigIdSeq:
            fout.write("%s " % (Id))
        fout.write("\n")    

    
    def get_len(self):
        self.get_string()
        if self.len != len(self.string):
            self.len = len(self.string)
        return self.len    
     

    def get_string(self): # whole string
        if len(self.contigIdSeq) == 0:
            return ""
        if len(self.string) > 0:
            return self.string
        firstId = self.contigIdSeq[0]
        #print Path.Contigs
        #string = Path.Contigs[firstId].string[-Path.StartLen:]
        #debug
        #two type string , final and mid
        string = Path.Contigs[firstId].string
        for Id in self.contigIdSeq[1:]:
            #print Id
            string = string + Path.Contigs[Id].string[Path.K:]
        self.string = string
        return string
    
    
    def get_local_path_string(self):
        if len(self.contigIdSeq) == 0:
            return ""
        if len(self.string) > 0:
            return self.string
        firstId = self.contigIdSeq[0]
        string = Path.Contigs[firstId].string[-Path.StartLen:]
        
        for Id in self.contigIdSeq[1:]:
            #print Id
            string = string + Path.Contigs[Id].string[Path.K:]
        self.string = string
        return string
    
    def get_map_path(self, mapLength):

        mapPath = []
        firstId = self.contigIdSeq[0]
        mapPath.append(firstId)
        mapString = Path.Contigs[firstId].string[-Path.StartLen:]
        pathLength = Path.StartLen
        
        for Id in self.contigIdSeq[1:]:
            assert Path.Contigs[Id].get_len() > Path.K
            pathLength += (Path.Contigs[Id].len - Path.K)
            if pathLength > mapLength:
                break
            mapString = mapString + Path.Contigs[Id].string[Path.K:]
            mapPath.append(Id)

        if len(mapPath) == 1:
            Id = self.contigIdSeq[1]
            mapString = mapString + Path.Contigs[Id].string[Path.K:]
            mapPath.append(Id)

        return mapPath, mapString
