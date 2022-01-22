
class MultipleAlign(object):
    """
    
    Alinha sequências multiplas.
    
    """
    def __init__(self, seqs, alignseq):
        self.seqs = seqs
        self.alignpars = alignseq
        
    def numSeqs(self):
        return len(self.seqs)
    
    def printMat (mat):
        for i in range(0, len(mat)):
            print(mat[i])
            
            
    def addSeqAlignment (self, alignment, seq):
            """
            
            Adiciona uma nova sequência a um alinhamento existente (fazendo o consenso do alinhamento).
            
            """
            res = [ ]
            for i in range(len(alignment.listseqs)  +1):
                res.append("")
            cons =  MySeq(alignment.consensus(),alignment.tipo)
            self.alignpars.needlemanWunsch(cons,    seq)
            align2 = self.alignpars.recoverAlignment()
            orig = 0
            for i in range(len(align2)): # cada coluna
                if align2[0,i]=='-':
                    for k in range(len(alignment.listseqs)):
                        res[k] += "-"
                else:
                    for k in range(len(alignment.listseqs)):
                        res[k] += alignment[k,orig]
                    orig+=1
            res[len(alignment.listseqs)] = align2.listseqs[1]
            return MyAlign(res, alignment.tipo)   
        
    def alignConsensus(self):
            """
            
            Implementa o algoritmo completo de alinhamento.
            
            """
            self.alignpars.needlemanWunsch(self.seqs[0], self.seqs[1])
            res = self.alignpars.recoverAlignment()
            for i in range(2, len(self.seqs)):
                res = self.addSeqAlignment(res, self.seqs[i])
            return res