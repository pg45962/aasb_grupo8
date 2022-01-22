class MyAlign:
    """
    
    Serve para guardar e manipular alinhamentos de sequências (duas ou mais sequências).
    
    """

    def __init__(self, lseqs, tipo="protein"):
        self.listseqs = lseqs
        self.tipo = tipo

    def __len__(self):  # number of columns
        return len(self.listseqs[0])

    def __getitem__(self, n):
        if type(n) is tuple and len(n) == 2:
            i, j = n
            return self.listseqs[i][j]
        elif type(n) is int:
            return self.listseqs[n]
        return None

    def __str__(self):
        res = ""
        for seq in self.listseqs:
            res += "\n" + seq
        return res

    def numSeqs(self):
        return len(self.listseqs)

    def column(self, indice):
        res = []
        for k in range(len(self.listseqs)):
            res.append(self.listseqs[k][indice])
        return res

    def consensus(self):
        """
        
        Usa o alinhamento representado com a lista de strings e retorna o caracter mais comum em cada coluna, ignorando espaçamentos.
        
        """
        
        cons = ""
        for i in range(len(self)):
            cont = {}
            for k in range(len(self.listseqs)):
                c = self.listseqs[k][i]
                if c in cont:
                    cont[c] = cont[c] + 1
                else:
                    cont[c] = 1
            maximum = 0
            cmax = None
            for ke in cont.keys():
                if ke != "-" and cont[ke] > maximum:
                    maximum = cont[ke]
                    cmax = ke
            cons = cons + cmax
        return cons