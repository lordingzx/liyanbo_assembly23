import tools

class Contig(object):
    def __init__(self, Id, string, copyNum):
        self.Id = Id
        self.copyNum = copyNum
        self.string = string
        self.len = len(string) 

    def print_string(self):
        print('%s: %s' % (self.Id, self.string))

    def reverse_contig(self):
        return Contig(reverse_id(self.Id), reverse_string(self.string), self.copy_num)

    def get_len(self):
        if self.len != len(self.string):
            self.len = len(self.string)
        return self.len

    def print_contig(self):
        print self.Id, self.copyNum, self.len
        #print self.string
        
