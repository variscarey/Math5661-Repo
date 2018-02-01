def dB(x,i,knot,degree):
    if degree==0:
        return np.zeros(x.shape())
    else:
        temp=degree/(knot[i+degree]-knot[i])*Brecur(x,i,knot,degree-1)
        temp-=degree/(knot[i+degree+1]-knot[i+1])*Brecur(x,i+1,knot,degree-1)
        return temp
    
def d2B(x,i,knot,degree):
    if degree<=1:
        return np.zeros(x.shape())
    else:
        temp=degree/(knot[i+degree]-knot[i])*dB(x,i,knot,degree-1)
        temp-=degree/(knot[i+degree+1]-knot[i+1])*dB(x,i+1,knot,degree-1)
        return temp
