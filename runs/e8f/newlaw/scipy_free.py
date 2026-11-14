import numpy as np
def nnls(X, y, iters=20000):
    """non-negative least squares by projected gradient (no scipy here)"""
    n=X.shape[1]; c=np.zeros(n); G=X.T@X; g=X.T@y; L=np.linalg.eigvalsh(G).max(); z=c.copy(); t=1.0
    for _ in range(iters):
        cn=np.maximum(z-(G@z-g)/L,0); tn=(1+np.sqrt(1+4*t*t))/2; z=cn+(t-1)/tn*(cn-c); c,t=cn,tn
    return c
