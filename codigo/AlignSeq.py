class AlignSeq:
    """
    
    Esta classe permite fazer o alinhamento entre duas sequências.
    
    """
    def __init__(self, sm, g):
            self.g = g
            self.sm = sm
            self.S = None
            self.T = None
            self.seq1 = None
            self.seq2 = None
            
    def printMat (mat):
        for i in range(0, len(mat)):
            print(mat[i])

    def scorePos (self, c1, c2):
        """
        
        Marca o score de uma posição.
        
        """
        if c1 == "-" or c2=="-":
            return self.g
        else:
            return self.sm[c1,c2]
            
    def scoreAlin (self, alin):
        """
        
        Marca o score de um alinhamento.
        
        """
        res = 0
        for i in range(len(alin)):
            res += self.scorePos (alin[0][i], alin[1][i])
        return res

    def needlemanWunsch (self, seq1, seq2, ties = False):
            """
            
            Permite realizar o alinhamento global de duas sequencias.
            
            """
            if (seq1.tipo != seq2.tipo): return None
            self.S = [[0]]
            self.T = [[0]]
            self.seq1 = seq1
            self.seq2 = seq2
            for j in range(1, len(seq2)+1):
                self.S[0].append(self.g * j)
                if ties: self.T[0].append([3])
                else: self.T[0].append(3)
            for i in range(1, len(seq1)+1):
                self.S.append([self.g * i])
                if ties: self.T.append([[2]])
                else: self.T.append([2])
            for i in range(0, len(seq1)):
                for j in range(len(seq2)):
                    s1 = self.S[i][j] + self.scorePos (seq1[i], seq2[j])
                    s2 = self.S[i][j+1] + self.g
                    s3 = self.S[i+1][j] + self.g
                    self.S[i+1].append(max(s1, s2, s3))
                    if ties:
                        self.T[i+1].append(max3t_with_ties(s1, s2, s3))
                    else:
                        self.T[i+1].append(max3t(s1, s2, s3))
            return self.S[len(seq1)][len(seq2)]

    def recoverAlignment (self):
            """
            
            Pega na matriz construída na última função, em conjunto com as duas sequencias e implementa o processo de recuperação de alinhamento ótimo global.
            
            """
            res = ["", ""]
            i = len(self.seq1)
            j = len(self.seq2)
            while i>0 or j>0:
                if self.T[i][j]==1:
                    res[0] = self.seq1[i-1] + res[0]
                    res[1] = self.seq2[j-1] + res[1]
                    i -= 1
                    j -= 1
                elif self.T[i][j] == 3:
                    res[0] = "-" + res[0]
                    res[1] = self.seq2[j-1] + res[1] 
                    j -= 1
                else:
                    res[0] = self.seq1[i-1] + res[0]
                    res[1] = "-" + res[1]
                    i -= 1
            return MyAlign(res, self.seq1.tipo)

    def recoverAlignment_with_ties (self):
            """
            
            Semelhante ao anterior mas tem em conta os empates.
            
            """
            i = len(self.seq1)
            j = len(self.seq2)  
            alins = [["", "", i,j]]
            res = []
            while alins:
                al = alins.pop(0)
                i = al[2]
                j = al[3]
                if i==0 and j==0:
                    res.append(al[:2])
                else:
                    for t in self.T[i][j]:
                        p = []
                        if t==1:
                            p.append(self.seq1[i-1] + al[0])
                            p.append(self.seq2[j-1] + al[1])
                            p.append(i-1)
                            p.append(j-1)
                        elif t == 3:
                            p.append("-" + al[0])
                            p.append(self.seq2[j-1] + al[1])
                            p.append(i)
                            p.append(j-1)
                        else:
                            p.append(self.seq1[i-1] + al[0])
                            p.append("-" + al[1])
                            p.append(i-1)
                            p.append(j)
                        alins.append(p)
            return res

    def smithWaterman (self, seq1, seq2, ties = False):
            """
            
            Algoritmo bem conhecido para a realização de alinhamento de sequencias local. Em vez de olhar a sequência total, o algoritmo de Smith-Waterman compara segmentos de todos os comprimentos possíveis e otimiza a medida de similaridade.
            
            """
            if (seq1.tipo != seq2.tipo): return None
            self.S = [[0]]
            self.T = [[0]]
            self.seq1 = seq1
            self.seq2 = seq2
            maxscore = 0
            for j in range(1, len(seq2)+1):
                self.S[0].append(0)
                if ties: self.T[0].append([0])
                else: self.T[0].append(0)
            for i in range(1, len(seq1)+1):
                self.S.append([0])
                if ties: self.T.append([[0]])
                else: self.T.append([0])
            for i in range(0, len(seq1)):
                for j in range(len(seq2)):
                    s1 = self.S[i][j] + self.scorePos(seq1[i], seq2[j]) 
                    s2 = self.S[i][j+1] + self.g
                    s3 = self.S[i+1][j] + self.g
                    b = max(s1, s2, s3)
                    if b <= 0:
                        self.S[i+1].append(0)
                        self.T[i+1].append(0)
                    else:
                        self.S[i+1].append(b)
                        if ties:
                            self.T[i+1].append(max3t_with_ties(s1, s2, s3))
                        else:
                            self.T[i+1].append(max3t(s1, s2, s3))
                        if b > maxscore: 
                            maxscore = b
            return maxscore

    def recoverLocalAlignment (self):
            """
            
            Pega na matriz construída na última função e implementa o processo de recuperação de alinhamento ótimo local.
            
            """
            res = ["", ""]
            maxscore = 0
            maxrow = 0
            maxcol = 0
            for i in range(1,len(self.S)):
                for j in range(1, len(self.S[i])):
                    if self.S[i][j] > maxscore:
                        maxscore = self.S[i][j]
                        maxrow = i
                        maxcol = j
            i = maxrow
            j = maxcol
            while i>0 or j>0:
                if self.T[i][j]==1:
                    res[0] = self.seq1[i-1] + res[0]
                    res[1] = self.seq2[j-1] + res[1]
                    i -= 1
                    j -= 1
                elif self.T[i][j] == 3:
                    res[0] = "-" + res[0];
                    res[1] = self.seq2[j-1] + res[1]; 
                    j -= 1
                elif self.T[i][j] == 2:
                    res[0] = self.seq1[i-1] + res[0];
                    res[1] = "-" + res[1]; 
                    i -= 1
                else: break
            return MyAlign(res, self.seq1.tipo)

    def recoverAlignLocal_with_ties (self):
            """
            
            Semelhante ao anterior mas tem em conta empates.
            
            """
            maxval = self.S[0][0]
            maxtups = []
            for i in range(0,len(self.S)):
                for j in range(0, len(self.S[i])):
                    if self.S[i][j] > maxval:
                        maxval = self.S[i][j]
                        maxtups = [(i,j)]
                    elif self.S[i][j] == maxval:
                        maxtups.append((i,j))
            alins = []
            for (i,j) in maxtups:
                alins.append(["", "", i,j])        
            res = []
            while alins:
                al = alins.pop(0)
                i = al[2]
                j = al[3]   
                if (i==0 and j==0) or (0 in self.T[i][j]):
                    res.append(al[:2])
                else:
                    for t in self.T[i][j]:
                        p = []
                        if t==1:
                            p.append(self.seq1[i-1] + al[0])
                            p.append(self.seq2[j-1] + al[1])
                            p.append(i-1)
                            p.append(j-1)
                        elif t == 3:
                            p.append("-" + al[0])
                            p.append(self.seq2[j-1] + al[1])
                            p.append(i)
                            p.append(j-1)
                        else:
                            p.append(self.seq1[i-1] + al[0])
                            p.append("-" + al[1])
                            p.append(i-1)
                            p.append(j)
                        alins.append(p)
            return res

    def max3t (v1, v2, v3):
        """
        
        Devolve o maior valor.
        
        """
        if v1 > v2:
            if v1 > v3: return 1
            else: return 3
        else:
            if v2 > v3: return 2
            else: return 3

    def max3t_with_ties(v1, v2, v3):
        """
        
        Devolve o maior valor mas tem em atenção empates.
        
        """
        if v1 > v2:
            if v1 > v3: 
                return [1]
            elif v1 == v3:
                return [1,3]
            else:
                return [3]
        elif v1 == v2:
            if v1 > v3: 
                return [1,2]
            elif v1 == v3:
                return [1,2,3]
            else:
                return [3]
        else:
            if v2 > v3: return [2]
            elif v2 == v3:
                return [2,3]
            else: return [3]
