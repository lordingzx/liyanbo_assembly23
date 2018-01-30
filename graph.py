import sys
import copy
import path
import tools
sys.setrecursionlimit(100000)

class ADJGraph(object): 
    
    Len= {}
    K = 0
    @classmethod
    def init_class_var(cls, contigs, k):
        for c in contigs:
            cls.Len[c] = contigs[c].get_len()
        cls.K = k

    
    def __init__(self, graph, graphR, edges, maxDeep, maxLen):
        self.graph = graph
        self.graphR = graphR
        self.maxDeep = maxDeep 
        self.maxLen = maxLen
        self.edges = edges
        if len(edges) != 0 and len(graph) == 0:
            self.edges_2_graph()

     
    def edges_2_graph(self):
        for (l, r) in self.edges:
            
            if l not in self.graph:
                self.graph[l] = []
            self.graph[l].append(r)
            if r not in self.graphR:
                self.graphR[r] = []
            self.graphR[r].append(l)   

            '''
            if l not in self.graph:
                self.graph[l] = []
                self.graph[l].append(r)
            elif len(self.graph[l]) == 1:
                if (self.edges[(l,r)] == self.edges[(l,self.edges[l][0])]):
                    self.graph[l] = []
                elif (self.edges[(l,r)] > self.edges[(l,self.edges[l][0])]):
                    self.graph[l].pop()
                    self.graph[l].append(r)

            if r not in self.graphR:
                self.graphR[r] = []
                self.graphR[r].append(l)   
            elif len(self.graphR[r]) == 1:
                if (self.edges[(l,r)] == self.edges[(self.edgesR[r][0], r)]):
                    self.graphR[r] = []
                elif (self.edges[(l,r)] > self.edges[(self.edgesR[r][0], r)]):
                    self.graphR[r].pop()
                    self.graphR[r].append(l)
            '''
        for key in self.graph:
            print key, self.graph[key]
        
        for key in self.graphR:
            print key, self.graphR[key]

 

    def link_one_in_one_out(self):
        print "bb"
        linkId = {}
        for (l, r) in self.edges:
            lR = tools.reverse_id(l)
            rR = tools.reverse_id(r)
            if (len(self.graph[l]) == 1 and        # l only one output 
                    len(self.graphR[r]) == 1 and lR not in self.graphR  # r only one input
                    and rR not in self.graph):
                assert l not in linkId
                linkId[l] = r
        paths = []
        visited = set()
        print "cc"
        for key in linkId:
            if key in visited:
                continue
            path = []
            path.append(key)
            visited.add(key)
            path.append(linkId[key])
            key = linkId[key]
            visited.add(key)
            while key in linkId:
                key = linkId[key]
                path.append(key)
                if key in visited:
                    break
                visited.add(key)

            paths.append(path)    
        
        print "aa", paths
        return paths        
 

      

    
    def calc_edgeNumber(self):
        edgeNumber = 0 
        for Id in self.graph:
            edgeNumber += len(self.graph[Id])
        print "[info] edge number: ", edgeNumber

    '''
    def divide_by_nodes(self, startContigsId):
        p = Pool()
    '''  
    
    def get_paths_start_one_node(self, Id): 
        print ("[info] start with Id: %s" % (Id))
        print ("[info] max deep: %s" % (self.maxDeep))
        paths = []
        self.dfs2(Id, "", paths, 0, 0) # dfs should change to length
        print "get paths start one node:", paths
        pathsList = []
        pathId = 0
        for p in paths:
            pathsList.append(path.Path(pathId, p.split(" ")))
            pathId += 1
        '''    
        for p in pathsList:
            p.print_seq()
        '''
        return pathsList


    def dfs(self, Id, temp, paths, deep): # by deep
        if deep + 1 == self.maxDeep:
            paths.append(temp + str(Id))
            return
        if Id in self.graph and len(self.graph[Id])>0:
            for adj in self.graph[Id]:
                self.dfs(adj, temp + str(Id) + " ", paths, deep+1)
        else:
            paths.append(temp + str(Id))
        return

    def dfs2(self, Id, temp, paths, lens, deep): # by length
        if deep + 1 == self.maxDeep: # set deep unper threshold, cause only control by length paths number maybe too much
            paths.append(temp + str(Id))
            return

        if lens > 0 and lens + ADJGraph.Len[Id] - ADJGraph.K >= self.maxLen:
            paths.append(temp + str(Id))
            return
        if Id in self.graph and len(self.graph[Id])>0:
            for adj in self.graph[Id]:
                if lens ==  0:
                    self.dfs2(adj, temp + str(Id) + " ", paths, self.maxLen/2, deep+1)
                else:
                    self.dfs2(adj, temp + str(Id) + " ", paths, lens + ADJGraph.Len[Id]-ADJGraph.K, deep+1)
        else:
            paths.append(temp + str(Id))
        return


    def addNode(self, Id, newId):
        self.graph[newId] = self.graph[Id]
        
'''
    def get_local_graph(self, Id, maxDeep):

    def find_path(self): # all paths from start to end
    
class SubGraph(ADJGraph):
    def __init__(self):
        self.start = set()
        self.end = set()
        self.graph = {}
        self.graphR = {}
'''
