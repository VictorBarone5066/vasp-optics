from matplotlib import pyplot as plt
import OpticsHeader as h

infile = "C://Users//baron//Desktop//opticsTest//OUTCAR_sicOptics"

#Make stuff easier to work with.  Each tensor function is named (chg/den)(Im/Re) depending on 
#what type of theory it was calculated from and if its the imaginary or real part.
#so chgIm is the imaginary part, calculated from the chg-chg theory
dieLis = h.GetDielecTensors(infile)
for i in range(0, len(dieLis)):
    if(dieLis[i][0].thry == "c-c"):
        if(dieLis[i][0].imag):
            chgIm = dieLis[i]
        if(dieLis[i][0].real):
            chgRe = dieLis[i]
    if(dieLis[i][0].thry == "d-d"):
        if(dieLis[i][0].imag):
            denIm = dieLis[i]
        if(dieLis[i][0].real):
            denRe = dieLis[i]

#Check the symmetries of the system.  If any of the output numbers are significant, will need
#to treat the system more specifically.  They'll probably all be small, though.  
h.CheckDielectricSyms(infile)

#Plot alpha(x) and refl(x)
x = [x_.energy for x_ in chgRe]
##Alpha
y = []
for r, i in zip(chgRe, chgIm):
    y.append(h.AbsCoeff(r, i)[0]*10**(-9)) #now in nm^-1       

x_, y_ = h.GetSplineInterp(x, y)
plt.xlim(0, 6)
plt.plot(x, y, 'ko')
plt.plot(x_, y_, 'c')
plt.show()
plt.clf()

##R
y = []
for r, i in zip(chgRe, chgIm):
    y.append(h.Reflec(r, i)[0])       

x_, y_ = h.GetSplineInterp(x, y)
plt.xlim(0, 6)
plt.plot(x, y, 'ko')
plt.plot(x_, y_, 'c')
plt.show()
plt.clf()


