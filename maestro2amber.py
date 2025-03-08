import sys
inputpdb = sys.argv[1]
outputpdb = sys.argv[2]
fr = open(inputpdb, "r")
fw = open(outputpdb, "w")
# only atom lines would be retained
# for GLU, if it has a H named HE2, change it's residue name as GLH
# for ASP, if it has a H named HD2, change it's residue name as ASH
# for LYS, if it doesn't have a H named HZ1, change it's residue name as LYN
# for CYS, if it doesn't have a H named HG, change it's residue name as CYX
# for HIS, three possiblities
# 1, have a H named HD1 and a H named HE2, change it's residue name as HIP
# 2, have a H named HD1, not a H HE2, HID
# 3, have a H named HE2, not a H HD1, HIE
# line[17:20] is residue name
# line[12:16] is atom name
# line[21:22] is chain id
# line[22:26] is residue number
class queue():
    def __init__(self):
        self.tail = -1
        self.Q = []
        self.isGlh = -1
        self.isAsh = -1
        self.isHie = -1
        self.isHid = -1
        self.isHip = -1
        self.isCys = -1
        self.isCy = -1
    def addQ(self, line):
        self.Q.append(line)
        self.tail += 1
        if line[17:20] == "GLU" and line[12:16].strip() == "HE2":
            #print (line)
            self.isGlh = 1
        elif line[17:20] == "ASP" and line[12:16].strip() == "HD2":
            #print (line)
            self.isAsh = 1
        elif line[17:20] == "HIS" and line[12:16].strip() == "HD1":
            #print (line)
            self.isHid = 1
        elif line[17:20] == "HIS" and line[12:16].strip() == "HE2":
            #print (line)
            self.isHie = 1
        elif line[17:20] == "CYS" and line[12:16].strip() == "HG":
            #print (line)
            self.isCys = 1
        else:
            pass
    def leaQ(self):
        if (self.tail!=-1):
            line = self.Q[0]
            self.Q = self.Q[1:self.tail+1]
            self.tail-=1
            if line[17:20] == 'CYS':
                self.isCy = 1
            if self.isGlh == 1:
                print("changed: " + line[0:17] + "GLH" + line[20:])
                return line[0:17] + "GLH" + line[20:]
            elif self.isAsh == 1:
                print("changed: " + line[0:17] + "ASH" + line[20:])
                return line[0:17] + "ASH" + line[20:]
            elif self.isCys == -1 and self.isCy == 1:
                print("changed: " + line[0:17] + "CYX" + line[20:])
                return line[0:17] + "CYX" + line[20:]
            elif self.isHid == 1 and self.isHie == 1:
                print("changed: " + line[0:17] + "HIP" + line[20:])
                return line[0:17] + "HIP" + line[20:]
            elif self.isHie == 1:
                print("changed: " + line[0:17] + "HIE" + line[20:])
                return line[0:17] + "HIE" + line[20:]
            elif self.isHid == 1:
                print("changed: " + line[0:17] + "HID" + line[20:])
                return line[0:17] + "HID" + line[20:]
            else:
                print("old one: " + line)
                return line
        else:
            self.isGlh = -1
            self.isAsh = -1
            self.isHie = -1
            self.isHid = -1
            self.isHip = -1
            self.isCys = -1
            self.isCy = -1
            return "nothing"
pdbQ = queue()
lrN = " "
lcI = " "
lrU = 999999
while True:
    line = fr.readline().strip()
    if not line:
        while True:
            lin1 = pdbQ.leaQ()
            if lin1 == "nothing":
                break
            #print (lin1)
            fw.write(lin1 + '\n')
        break
    if line[0:3] == "TER":
        while True:
            lin1 = pdbQ.leaQ()
            if lin1 == "nothing":
                break
            #print (lin1)
            fw.write(lin1 + '\n')
        fw.write(line + '\n')
        continue
    if line[0:4]!='ATOM' and line[0:6]!='HETATM':
        continue
    resName = line[17:20]
    atomName = line[12:16].strip()
    chainId = line[21:22]
    resNum = int(line[22:26])
    if ((resName != lrN) or (chainId != lcI) or (resNum != lrU)):
        # new residue, pop old residue from the queue
        while True:
            lin1 = pdbQ.leaQ()
            if lin1 == "nothing":
                break
            #print (lin1)
            fw.write(lin1 + '\n')
        pdbQ.addQ(line)
    else:
        pdbQ.addQ(line)
    lrN = resName
    lcI = chainId
    lrU = resNum
