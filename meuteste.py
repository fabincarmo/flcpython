
import flc
import numpy as np
import matplotlib.pyplot as pl
from scipy.io import mmread

def testshiftm():
    x = np.random.randint(10,99,(5,5)).astype(float)
    print x.astype(int)
    print "\n"
    y = flc.mshift(x, -2, "lr")
    print y.astype(int)
    return

def testbmap():
    t = np.random.randint(1,99,(100)).astype(float)
    flc.bmap(t, "Df.txt")
    D = mmread("Df.txt")
    print D.getnnz()
    print D.shape[0]*D.shape[1]
    pl.imshow(D.todense(), interpolation='none')
    pl.show()
    return

def testkmap2():
    t = np.linspace(0,1,20).tolist()
    t += np.linspace(0,1,40).tolist()
    t += np.linspace(0,1,60).tolist()
    t += np.linspace(0,1,80).tolist()
    t = np.asarray(t)
    flc.kmap2(t, 1e-2, 1e-2, "Dt.txt")
    D = mmread("Dt.txt")
    print D.getnnz()
    print D.shape[0]*D.shape[1]
    pl.imshow(D.todense(), interpolation='none')
    pl.show()
    return

if __name__=="__main__":
#    testshiftm()
#    testbmap()
#    testkmap2()