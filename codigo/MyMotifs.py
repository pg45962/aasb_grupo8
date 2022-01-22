
class MyMotifs:
    """

    Implementa a estrutura de dados e metodos precisos para criar a PWM para determinar diferentes representações deterministicas dos motifs e para determinar a probabilidade de ocorrência do motif ao longo da sequência.
    
    """
    def __init__(self, seqs):
        self.tam = len(seqs[0])
        self.seqs = seqs 
        self.alfabeto = seqs[0].alfabeto()
        self.calculaContagens()
        self.criaPWM()
        
    
    def calculaContagens(self):
        """
        
        Cria a matriz de contagens a partir das instâncias (seqs).
        
        """
        self.counts = createMatZeros(len(self.alfabeto), self.tam)
        for s in self.seqs:
            for i in range(self.tam):
                lin = self.alfabeto.index(s[i])
                self.counts[lin][i] += 1
    
    
    def criaPWM(self):
        """
        
        Cria a matriz pwm a partir das contagens.
        
        """
        if self.counts == None: self.calculaContagens()
        self.pwm = createMatZeros(len(self.alfabeto), self.tam)
        for i in range(len(self.alfabeto)):
            for j in range(self.tam):
                self.pwm[i][j] = float(self.counts[i][j]) / len(self.seqs)
                
                
    def consenso(self):
        """

        Define o consenso que é dado pelos caracteres mais comuns em cada posição (coluna) do motif.
        
        """
        res = ""
        for j in range(self.tam):
            maxcol = self.counts[0][j]
            maxcoli = 0
            for i in range(1, len(self.alfabeto) ):
                if self.counts[i][j] > maxcol: 
                    maxcol = self.counts[i][j]
                    maxcoli = i
            res += self.alfabeto[maxcoli]        
        return res
    
    
    def probabSeq (self, seq):

        """Dá a probabilidade de uma dada sequência inst do mesmo tamanho das instâncias ser gerada pelo motif."""

        res = 1.0
        for i in range(self.tam):
            lin = self.alfabeto.index(seq[i])
            res *= self.pwm[lin][i]
        return res
    
    
    def mostProbableSeq(self, seq):

        """Dá a posição mais provável de ocorrência do motif numa sequência seq de tamanho maior do que o motif."""

        maximo = -1.0
        maxind = -1
        for k in range(len(seq)-self.tam):
            p = self.probabSeq(seq[k:k+ self.tam])
            if(p > maximo):
                maximo = p
                maxind = k
        return maxind
    
    
        
