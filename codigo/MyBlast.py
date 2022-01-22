def readDatabase(filename):
    """

    Lê um conjunto de sequencias num ficheiro texto, uma por linha.
    
    """
    file= open (filename, "r+")
    db=[]
    for line in file.readlines():
        db.append(line)
    file.close()
    return db

def buildMap(query,w):
    """
    
    Dada a sequência (**query**) e o tamanho das palavras (**W**) mapea cada ocorrência de sub-sequências de tamanho W. O resultado é um dicionário com:

        - Chaves: sub-sequencias de tamanho W (que ocorrem na query);

        - Valores: lista de indices (posições onde a sub-sequência ocorre na query).
    
    """
    res={}
    i=0
    while i <= len(query)-w:
        q=query[i:i+w]
        if q not in res:
            res[q]=[i]
        else:
            res[q].append(i)
        i+=1
    return res

def getHits(seq,m,w):
    """
    
    Tem como entradas o Mapa da sequência query (criada pela função anterior), Sequência na qual procurar (e.g. da BD) e Valor do tamanho das palavras (W). O seu objetivo é descobrir todas as co-ocorrências de palavras de tamanho W entre a query e a sequência (só serão considerados matches perfeitos de tamanho W). O resultado é uma lista de hits; cada hit é um tuplo, com o 1º elemento sendo a posição da palavra na query e o 2º sendo a posição da palavra na sequência.
     
    """
    res=[]
    i=0
    while i<=len(seq)-w:
        s=seq[i:i+w]
        if s in m:
            for q in m[s]:
                res.append((q,i))
        i+=1
    return res

def extendHit (seq,hit,query,w):
    """
    Cada hit será estendido em ambas as direções. Um hit será estendido enquanto o nº de matches nessa direção for pelo menos metade do tamanho da extensão. A função terá como entradas: a sequência (da bd), a query, o hit e o valor de W. O resultado será um tuplo com 4 valores: início do hit estendido na query, início do hit estendido na sequência, tamanho do hit estendido, score do hit estendido (= nº de matches).
    
    Parâmetros:

    - stq: início do hit na query;
    - sts: início do hit na sequência;
    - matfw: matches à frente;
    - k: posições andadas em frente;
    - size: tamanho do hit estendido.

    """
    stq,sts =hit[0],hit[1]
    matfw =0
    k=0
    bestk =0
    while 2*matfw >=k and stq+w+k < len(query) and sts+w+k < len(seq):
        if query[stq+w+k] == seq[sts+w+k]:
            matfw+=1
            bestk = k+1
        k += 1
    size = w + bestk
    k = 0
    matbw = 0
    bestk = 0
    while 2*matbw >= k and stq > k and sts > k:
        if query[stq-k-1] == seq[sts-k-1]:
            matbw+=1
            bestk = k+1
        k+=1
    size += bestk
    return (stq-bestk, sts-bestk, size, w+matfw+matbw)

def hitBestScore (seq,query,m,w):
    """
    
    Devolve o hit com maior score de uma sequencia com a query. Em caso de empate devolve a de menor tamanho.
    
    """
    hits=getHits(seq,m,w)
    bestScore =-1.0
    best=()
    for h in hits:
        ext =extendHit(seq,h,query,w)
        score=ext[3]
        if score > bestScore or (score== bestScore and ext[2] < best[2]):
            bestScore = score
            best=ext
    return best