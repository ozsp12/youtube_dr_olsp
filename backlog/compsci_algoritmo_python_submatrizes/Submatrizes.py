import numpy as np
A=np.array(
    [[1,2,3,3],
     [3,4,5,2],
     [6,7,8,9],
     [9,10,11,2]]
)
a=3
b=2
def Submatrizes(A,a,b):
    '''Retorna todas as submatrizes contiguas axb de uma matriz A
        Parametros
        ----------
        A:Matriz base
        a:numero de linhas de cada submatriz
        b:numero de colunas de cada submatriz
        Retorna
        ----------
        Um tensor com cada submatriz e o numero teórico de submatrizes
        '''
    n=A.shape[0]
    m=A.shape[1]
    # Calculo teorico do numero de submatrizes
    N=(n+1-a)*(m+1-b)
    # Calculo do tensor das submatrizes
    B=np.zeros([N,a,b])
    k=0
    for i in range(n+1-a):
        for j in range(m+1-b):
            B[k]=A[i:a+i,j:b+j]
            k+=1

    return B, N
B,N=Submatrizes(A,a,b)
print('As submatrizes:\n',B)
print('Numero teórico de submatrizes:',N)