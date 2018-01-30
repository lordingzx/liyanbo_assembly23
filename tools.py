import sys

def reverse_id(node):
    if node.startswith("R"):
        r_node = node[1:]
    else:
        r_node = "R"+node
    return r_node

def reverse_string(string):
    ans = ""
    for i in range(len(string)-1, -1, -1):
        if string[i] == 'A':
            ans += 'T'
            continue
        elif string[i] == 'T':
            ans += 'A'
            continue
        elif string[i] == 'C':
            ans += 'G'
            continue
        elif string[i] == 'G':
            ans += 'C'
            continue
        else:
            print string[i]
            os._exit("not ATCG")
    return ans    

def reverse_seq(contigSeq):
    res = []
    for i in range(0, len(contigSeq)):
        key = reverse_id(contigSeq[i])
        res.insert(0,key)
    return res    
