#For Optical Calculations

class gridPoint:
    real = None
    imag = None
    thry = None
    
    energy = None
    x = None
    y = None
    z = None
    xy = None
    yz = None
    zx = None
    
    def __init__(self, cmplx, line, typ = "undf"):
        self.real = None
        self.imag = None
        self.thry = None
        
        self.energy = None
        self.x = None
        self.y = None
        self.z = None
        self.xy = None
        self.yz = None
        self.zx = None
        
        if(cmplx[0] == 'r' or cmplx[0] == 'R'):
            self.real = True
            self.imag = False
        elif(cmplx[0] == 'i' or cmplx[0] == 'I'):
            self.real = False
            self.imag = True
            
        self.thry = typ
        
        self.energy = float(line.split()[0])
        self.x = float(line.split()[1])
        self.y = float(line.split()[2])
        self.z = float(line.split()[3])
        self.xy = float(line.split()[4])
        self.yz = float(line.split()[5])
        self.zx = float(line.split()[6])
        
#------------------------------------------------------------------------------------------------- 
#Get Tensors from OUTCAR file
#-------------------------------------------------------------------------------------------------        
import re
def GetLinesContainingPhrase(filePath, key, lineNums = False):
    phrases = []
    lines = []
    count = 0
    with open (filePath, 'r') as infile: ##scan infile as read-only line-by-line
        for line in infile:
            count = count + 1
            if(re.search(key, line)):
                phrases.append(line)
                lines.append(count)
    infile.close()
    
    if(len(phrases) == 0):
        return False
    if(lineNums):
        return phrases, lines
    return phrases

#Returns all of the dielectric tensors that it can find from OUTCAR (presumably).  Will be stored
#like this:
#[tensor0,           tensorI = [gridPoint0,          N = number of tensors
# tensor1,                      gridPoint1,          M = NEDOS? It should be, but isnt...
# ...                           ...
# tensorN]                      gridPointM]
def GetDielecTensors(outcarLoc="OUTCAR", key=r"frequency dependent .* DIELECTRIC FUNCTION"):
    try:
        phrs, nums = GetLinesContainingPhrase(outcarLoc, key, lineNums = True)
    except(TypeError):
        print("Unable to find the key phrase in the file search location")
        return False
    nums.append(-1)
    phrs.append(-1)
    tens = []

    with open(outcarLoc, 'r') as infile:
        read = False
        for i, line in enumerate(infile):
            #Determins if we should begin reading in function info. Trust me, it works. Don't touch  
            if(not read and i >= nums[len(tens)] and nums[len(tens)] != -1):
                if(line.split()[0][0] != 'E' and line.split()[0][0:2] != "--"):
                    read = True
                    thisTens = []
                    ##Determine if real or imaginary part
                    if(re.search("REAL", phrs[len(tens)])):
                        complexType = "real"
                    elif(re.search("IMAGINARY", phrs[len(tens)])):
                        complexType = "imag"
                    else:
                        complexType = "undf :("
                    ##Determine if its density-density or charge-charge
                    if(re.search("density-density", phrs[len(tens)])):
                        thryType = "d-d"
                    elif(re.search("current-current", phrs[len(tens)])):
                        thryType = "c-c"
                    else:
                        thryType = "undf :("
            if(read and line.split() == []):
                read = False
                tens.append(thisTens)
            #Actually read in function info
            if(read and line.split()[0][0] != 'T'):
                thisTens.append(gridPoint(complexType, line, thryType))
                
        infile.close()
    return tens
            
#-------------------------------------------------------------------------------------------------    
#Calculations of interest:
#Reference: "Optical Properties of Solids" by Mark Fox, Oxford U Press
#-------------------------------------------------------------------------------------------------
hbar = 1.05457163e-34 #Js
sol = 299792458 #m/s
def sqrt(n):
    return (n)**(1.0/2.0)

def EvToJo(e):   #dont ask who Jo(e) is
    return e*1.60218e-19

#Returns the generally nonsymmetric parts of the normal refractive index tensor
def RefInd(re, im):
    if(re.energy != im.energy):
        print("ERR: The energy of the real part does not match the energy of the imaginary part!")
        return None
    n = [None]*6 #Will be nX, nY, nZ, nXY, nYZ, nZX
    n[0] = sqrt(sqrt(re.x**2 + im.x**2) + re.x)/sqrt(2)
    n[1] = sqrt(sqrt(re.y**2 + im.y**2) + re.y)/sqrt(2)
    n[2] = sqrt(sqrt(re.z**2 + im.z**2) + re.z)/sqrt(2)
    n[3] = sqrt(sqrt(re.xy**2 + im.xy**2) + re.xy)/sqrt(2)
    n[4] = sqrt(sqrt(re.yz**2 + im.yz**2) + re.yz)/sqrt(2)
    n[5] = sqrt(sqrt(re.zx**2 + im.zx**2) + re.zx)/sqrt(2)

    return n

#Returns the generally nonsymmetric parts of the extinction coefficient tensor
def ExtCoeff(re, im):
    if(re.energy != im.energy):
        print("ERR: The energy of the real part does not match the energy of the imaginary part!")
        return None
    k = [None]*6 #Will be nX, nY, nZ, nXY, nYZ, nZX
    k[0] = sqrt(sqrt(re.x**2 + im.x**2) - re.x)/sqrt(2)
    k[1] = sqrt(sqrt(re.y**2 + im.y**2) - re.y)/sqrt(2)
    k[2] = sqrt(sqrt(re.z**2 + im.z**2) - re.z)/sqrt(2)
    k[3] = sqrt(sqrt(re.xy**2 + im.xy**2) - re.xy)/sqrt(2)
    k[4] = sqrt(sqrt(re.yz**2 + im.yz**2) - re.yz)/sqrt(2)
    k[5] = sqrt(sqrt(re.zx**2 + im.zx**2) - re.zx)/sqrt(2)

    return k    

#Returns attenuation/absorption coefficient (tensor) in SI base units
def AbsCoeff(re, im):
    if(re.energy != im.energy):
        print("ERR: The energy of the real part does not match the energy of the imaginary part!")
        return None
    
    a = [None]*6 #Will be alphaX, alphaY, alphaZ, alphaXY, alphaYZ, alphaZX
    kap = ExtCoeff(re, im)
    for i in range(0, 6):
        a[i] = 2*EvToJo(re.energy)/(hbar*sol)*kap[i]
        
    return a
        
#Returns reflectivity (tensor) in arbitrary units
def Reflec(re, im):
    if(re.energy != im.energy):
        print("ERR: The energy of the real part does not match the energy of the imaginary part!")
        return None
    
    r = [None]*6 #Will be rX, rY, rZ, rXY, rYZ, rZX       
    n = RefInd(re, im)
    kap = ExtCoeff(re, im)
    for i in range(0, 6):
        r[i] = ((n[i] - 1)**2 + kap[i]**2)/((n[i] + 1)**2 + kap[i]**2)
        
    return r
  

#Get the real and imaginary parts of the functions
def CheckDielectricSyms(infileLoc, key=r"frequency dependent .* DIELECTRIC FUNCTION"):
    dieLis = GetDielecTensors(infileLoc, key)
        
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

    #For the Chg-Chg method, get the deviation in rel perms of the real part from the average
    print("Chg-Chg Method")
    xs = [xs_.x for xs_ in chgRe]
    ys = [ys_.x for ys_ in chgRe]
    zs = [zs_.x for zs_ in chgRe]
    avgs = [1.0/3.0*(a+b+c) for a, b, c in zip(xs, ys, zs)]
    print("Max Deviation from Average in X (real): " + 
          str(max([x - avg for x, avg in zip(xs, avgs)])))
    print("Max Deviation from Average in Y (real): " + 
          str(max([y - avg for y, avg in zip(ys, avgs)])))
    print("Max Deviation from Average in Z (real): " + 
          str(max([z - avg for z, avg in zip(zs, avgs)])))
    #And now the imaginary
    xs = [xs_.x for xs_ in chgIm]
    ys = [ys_.x for ys_ in chgIm]
    zs = [zs_.x for zs_ in chgIm]
    avgs = [1.0/3.0*(a+b+c) for a, b, c in zip(xs, ys, zs)]
    print("Max Deviation from Average in X (imag): " + 
          str(max([x - avg for x, avg in zip(xs, avgs)])))
    print("Max Deviation from Average in Y (imag): " + 
          str(max([y - avg for y, avg in zip(ys, avgs)])))
    print("Max Deviation from Average in Z (imag): " + 
          str(max([z - avg for z, avg in zip(zs, avgs)])))
    
    #For the Den-Den method, get the deviation in rel perms of the real part from the average
    print("\nDen-Den Method")
    xs = [xs_.x for xs_ in denRe]
    ys = [ys_.x for ys_ in denRe]
    zs = [zs_.x for zs_ in denRe]
    avgs = [1.0/3.0*(a+b+c) for a, b, c in zip(xs, ys, zs)]
    print("Max Deviation from Average in X (real): " + 
          str(max([abs(x - avg) for x, avg in zip(xs, avgs)])))
    print("Max Deviation from Average in Y (real): " + 
          str(max([abs(y - avg) for y, avg in zip(ys, avgs)])))
    print("Max Deviation from Average in Z (real): " + 
          str(max([abs(z - avg) for z, avg in zip(zs, avgs)])))
    #And now the imaginary
    xs = [xs_.x for xs_ in denIm]
    ys = [ys_.x for ys_ in denIm]
    zs = [zs_.x for zs_ in denIm]
    avgs = [1.0/3.0*(a+b+c) for a, b, c in zip(xs, ys, zs)]
    print("Max Deviation from Average in X (imag): " + 
          str(max([x - avg for x, avg in zip(xs, avgs)])))
    print("Max Deviation from Average in Y (imag): " + 
          str(max([y - avg for y, avg in zip(ys, avgs)])))
    print("Max Deviation from Average in Z (imag): " + 
          str(max([z - avg for z, avg in zip(zs, avgs)])))
    
    #Check the difference between values of rel perm between chg-chg and den-den methods
    print("\nReal")
    dx, cx = [dx_.x for dx_ in chgRe], [cx_.x for cx_ in chgRe]
    dy, cy = [dy_.y for dy_ in chgRe], [cy_.y for cy_ in chgRe]
    dz, cz = [dz_.z for dz_ in chgRe], [cz_.z for cz_ in chgRe]
    print("Max Deviation Between Den-Den and Chg-Chg Methods in X: " + 
          str(max(abs(cx_ - dx_) for cx_, dx_ in zip(cx, dx))))
    print("Max Deviation Between Den-Den and Chg-Chg Methods in Y: " + 
          str(max(abs(cy_ - dy_) for cy_, dy_ in zip(cy, dy))))
    print("Max Deviation Between Den-Den and Chg-Chg Methods in Z: " + 
          str(max(abs(cz_ - dz_) for cz_, dz_ in zip(cz, dz))))
    print("\nImaginary")
    dx, cx = [dx_.x for dx_ in chgIm], [cx_.x for cx_ in chgIm]
    dy, cy = [dy_.y for dy_ in chgIm], [cy_.y for cy_ in chgIm]
    dz, cz = [dz_.z for dz_ in chgIm], [cz_.z for cz_ in chgIm]
    print("Max Deviation Between Den-Den and Chg-Chg Methods in X: " + 
          str(max([abs(cx_ - dx_) for cx_, dx_ in zip(cx, dx)])))
    print("Max Deviation Between Den-Den and Chg-Chg Methods in Y: " + 
          str(max([abs(cy_ - dy_) for cy_, dy_ in zip(cy, dy)])))
    print("Max Deviation Between Den-Den and Chg-Chg Methods in Z: " + 
          str(max([abs(cz_ - dz_) for cz_, dz_ in zip(cz, dz)])))

    return
      
#-------------------------------------------------------------------------------------------------
#Plotting Utilities
#Specifically cubic spline interpolation
#-------------------------------------------------------------------------------------------------
from scipy.interpolate import CubicSpline
import numpy as np

#Returns interpolated x, interpolated f(x)
def GetSplineInterp(x, y, den = 0.01):
    cub = CubicSpline(x, y)
    x_ = np.arange(x[0], x[-1], den)
    return list(x_), list(cub(x_))
      


