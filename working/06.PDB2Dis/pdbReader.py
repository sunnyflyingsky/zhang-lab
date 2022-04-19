from settings import PDB_Entry,RingData,MAXNUM,RAD
import math
import re



class PDB(object):
    '''
    用于加载PDB结构，这是一个自定义的PDB结构解析器，如果你有更好的选择，比如调包，可以尝试更可靠的方式
    '''
    def __init__(self,filePath:str,pdbType:str='RNA'):
        self.ATOM_info = []  #简单序列排列的原子列表
        self.sequenceList = {}  #链
        self.Conformers = []  ##chain_id->atomnum->ATOM_PDB
        self.Conformers_Atom = []  ##chain_id->resnum->atomname->ATOM_PDB
        self.residueList = []  ##chain_id->resnum->resname
        self.EMPTY = PDB_Entry()  #储存空的原子结构
        itnum = 0
        flag = 0
        if pdbType=='RNA':
            with open(filePath) as read_object:
                for line in read_object:
                    if line.startswith('ENDMDL'):
                        break
                    if flag and line.startswith('ATOM'):
                        flag = 0
                    if flag: break
                    if line.startswith('ATOM') or line.startswith('HETATM'):
                        itnum+=1
                        num = int(self.getField(line,2))
                        name = self.getField(line,3)
                        residue = self.getField(line,4)
                        chain_id = self.getField(line,5)
                        seq = int(self.getField(line,6))
                        X = float(self.getField(line,7))
                        Y = float(self.getField(line,8))
                        Z = float(self.getField(line,9))
                        #occupancy = getField(line.strip(),10)
                        b_factor = float(self.getField(line.strip(),11))
                        pdb_temp = PDB_Entry(atomNum=num,resNum=seq,atomName=name,resName=residue,\
                            chainName=chain_id,X=X,Y=Y,Z=Z,B_Factor=b_factor,Coord=[X,Y,Z])
                        self.ATOM_info.append(pdb_temp)
                        if chain_id in self.sequenceList.keys():
                            self.sequenceList[chain_id]+=residue
                            self.Conformers[-1][itnum] = PDB_Entry(atomNum=itnum,resNum=seq,atomName=name,resName=residue,\
                                chainName=chain_id,X=X,Y=Y,Z=Z,B_Factor=b_factor,Coord=[X,Y,Z])
                            if seq not in self.Conformers_Atom[-1].keys():
                                self.Conformers_Atom[-1][seq] = {}
                            self.Conformers_Atom[-1][seq][name] = PDB_Entry(atomNum=itnum,resNum=seq,atomName=name,resName=residue,\
                                chainName=chain_id,X=X,Y=Y,Z=Z,B_Factor=b_factor,Coord=[X,Y,Z])
                            self.residueList[-1][seq] = residue
                        else:
                            self.sequenceList[chain_id]=residue
                            self.Conformers.append({})
                            self.Conformers_Atom.append({})
                            self.Conformers[-1][itnum] = PDB_Entry(atomNum=itnum,resNum=seq,atomName=name,resName=residue,\
                                chainName=chain_id,X=X,Y=Y,Z=Z,B_Factor=b_factor,Coord=[X,Y,Z])
                            if seq not in self.Conformers_Atom[-1].keys():
                                self.Conformers_Atom[-1][seq] = {}
                            self.Conformers_Atom[-1][seq][name] = PDB_Entry(atomNum=itnum,resNum=seq,atomName=name,resName=residue,\
                                chainName=chain_id,X=X,Y=Y,Z=Z,B_Factor=b_factor,Coord=[X,Y,Z])
                            self.residueList.append({seq:residue})
                    elif line.startswith('TER'):
                        flag = 1
                    elif line.startswith(''):
                        pass
        elif pdbType=='ligand':
            with open(filePath) as read_object:
                for line in read_object:
                    if line.startswith('ENDMDL'):
                        break
                    if flag and line.startswith('ATOM'):
                        flag = 0
                    if flag and line.startswith('HETATM'):
                        itnum+=1
                        num = int(self.getField(line,2))
                        name = self.getField(line,3)
                        residue = self.getField(line,4)
                        chain_id = self.getField(line,5)
                        seq = int(self.getField(line,6))
                        X = float(self.getField(line,7))
                        Y = float(self.getField(line,8))
                        Z = float(self.getField(line,9))
                        #occupancy = getField(line.strip(),10)
                        b_factor = float(self.getField(line.strip(),11))
                        pdb_temp = PDB_Entry(atomNum=num,resNum=seq,atomName=name,resName=residue,\
                            chainName=chain_id,X=X,Y=Y,Z=Z,B_Factor=b_factor,Coord=[X,Y,Z])
                        self.ATOM_info.append(pdb_temp)
                        if chain_id in self.sequenceList.keys():
                            self.sequenceList[chain_id].append(residue)
                            self.Conformers[-1][itnum] = PDB_Entry(atomNum=itnum,resNum=seq,atomName=name,resName=residue,\
                                chainName=chain_id,X=X,Y=Y,Z=Z,B_Factor=b_factor,Coord=[X,Y,Z])
                            if seq not in self.Conformers_Atom[-1].keys():
                                self.Conformers_Atom[-1][seq] = {}
                            self.Conformers_Atom[-1][seq][name] = PDB_Entry(atomNum=itnum,resNum=seq,atomName=name,resName=residue,\
                                chainName=chain_id,X=X,Y=Y,Z=Z,B_Factor=b_factor,Coord=[X,Y,Z])
                            self.residueList[-1][seq] = residue
                        else:
                            self.sequenceList[chain_id]=[residue]
                            self.Conformers.append({})
                            self.Conformers_Atom.append({})
                            self.Conformers[-1][itnum] = PDB_Entry(atomNum=itnum,resNum=seq,atomName=name,resName=residue,\
                                chainName=chain_id,X=X,Y=Y,Z=Z,B_Factor=b_factor,Coord=[X,Y,Z])
                            if seq not in self.Conformers_Atom[-1].keys():
                                self.Conformers_Atom[-1][seq] = {}
                            self.Conformers_Atom[-1][seq][name] = PDB_Entry(atomNum=itnum,resNum=seq,atomName=name,resName=residue,\
                                chainName=chain_id,X=X,Y=Y,Z=Z,B_Factor=b_factor,Coord=[X,Y,Z])
                            self.residueList.append({seq:residue})
                    elif line.startswith('TER'):
                        flag = 1
                    
                        
    #####
    # 载入PDB参数
    #####

    def getField(self,text:str,index:int):
        '''
        从PDB的原子field信息中获取对应的具体信息的方法
        '''
        if index==1:return text[0:6].strip()  ##"ATOM"
        elif index==2:return text[6:11].strip()  ##Atom number
        elif index==3:return text[11:16].strip()  ##Atom name
        elif index==4:return text[17:21].strip()  ##Residue name
        elif index==5:return text[21:22].strip()  ##Chain ID
        elif index==6:return text[22:26].strip()  ##Residue seq
        elif index==7:return text[30:38].strip()  ##X
        elif index==8:return text[38:46].strip()  ##Y
        elif index==9:return text[46:54].strip()  ##Z
        elif index==10:return text[54:60].strip()  ##Occupancy
        elif index==11:return text[60:66].strip()  ##B-factor
        else:
            return None
    
    def getEntry(self,conformerID:int,rNum:int,aName:str=''):
        '''
        获取对应的原子
        '''
        if aName=='':
            return self.Conformers[conformerID][rNum]
        else:
            return self.Conformers_Atom[conformerID][rNum][aName]
    
    #####
    # 计算氨基酸的键转角,
    #####
    def getBondAngle(self, A: list, B: list, C: list):
        '''
        得到键键之间的夹角
        获得的是角度值
        '''
        v,v1,v2 = B.copy(),A.copy(),C.copy()
        self.Vec3Sub(v1,v) #向量BA
        self.Vec3Sub(v2,v) #向量BC
        self.Vec3Norm(v1)
        self.Vec3Norm(v2)
        c = self.Vec3Scalar(v1,v2)  #内积
        self.Vec3Cross(v1,v2)  #外积
        s = self.Vec3Abs(v1)

        return math.atan2(s,c)*RAD
    
    def getBondAngle_2(self, a: PDB_Entry, b: PDB_Entry, c:PDB_Entry):
        '''
        函数重载，获取原子之间的键角
        '''
        return self.getBondAngle(a.Coord,b.Coord,c.Coord)
    
    def getDihedralAngle(self,a: PDB_Entry,b: PDB_Entry,c:PDB_Entry,d:PDB_Entry):
        '''
        求平面夹角
        返回的是角度值
        '''
        cb,n1,n2 = [0,0,0],[0,0,0],[0,0,0]
        co = 0
        if a.atomName == '' or  b.atomName == '' or c.atomName == '' or  d.atomName == '':
            return MAXNUM
        #bc向量
        cb[0],cb[1],cb[2] = c.X-b.X, c.Y-b.Y,c.Z-b.Z
        ##平面abc的法向量
        n1[0] = (b.Y - a.Y) * cb[2] + (a.Z - b.Z) * cb[1]
        n1[1] = (b.Z - a.Z) * cb[0] + (a.X - b.X) * cb[2]
        n1[2] = (b.X - a.X) * cb[1] + (a.Y - b.Y) * cb[0]
        ##平面bcd的法向量
        n2[0] = cb[1] * (d.Z - c.Z) + cb[2] * (c.Y - d.Y)
        n2[1] = cb[2] * (d.X - c.X) + cb[0] * (c.Z - d.Z)
        n2[2] = cb[0] * (d.Y - c.Y) + cb[1] * (c.X - d.X)
        TEMP,TEMP1,TEMP2,TEMP3,TEMP4,TEMP5 = \
            n1[0],n1[1],n1[2],n2[0],n2[1],n2[2]
        #计算两个平面的夹角，即abc平面与bcd平面的夹角的cos值
        #cos(a-b-c-d) = |n1·n2|/||*||
        co = float((n1[0]*n2[0] + n1[1]*n2[1] + n1[2]*n2[2])/\
            math.sqrt((TEMP**2+TEMP1**2+TEMP2**2)*(TEMP3**2+TEMP4**2+TEMP5**2)))
        return RAD * (math.acos(co) * \
            self.sgn(float((n1[1]*n2[2]-n1[2]*n2[1])*cb[0] + \
            (n1[2]*n2[0]-n1[0]*n2[2])*cb[1] + \
            (n1[0]*n2[1]-n1[1]*n2[0])*cb[2])))
        ##这个正负取得是什么值，两个法向量与BC垂直，则

    def getDist(self,A:list,B:list):
        '''
        计算距离
        '''
        v = [0,0,0]
        v = B.copy()
        self.Vec3Sub(v,A)
        return self.Vec3Abs(v)

    #####
    # 某一些运算符
    #####
    def sgn(self,x):
        '''
        正负指示器
        '''
        return ((x >= 0) - (x < 0))
    
    def Vec3Zero(self,v:list):
        '''
        向量归零
        '''
        for i in range(len(v)):
            v[i] = float(0.0)
    
    def Vec3Abs(self, v: list):
        '''
        取模
        '''
        return math.sqrt(v[0]**2+v[1]**2+v[2]**2)
    
    def Vec3DiffAbs(self,v1:list, v2:list):
        '''
        向量终点之间的距离
        '''
        return math.sqrt((v1[0]-v2[0])**2 + \
            (v1[1]-v2[1])**2 + \
            (v1[2]-v2[2])**2 )
 
    def Vec3Norm(self, v: list):
        '''
        归一化
        '''
        a = self.Vec3Abs(v)
        if a != 0:
            for i in range(len(v)):
                v[i] /= a

    def Vec3Scale(self,v:list,s:float):
        '''
        向量伸缩
        '''
        for i in range(len(v)):
            v[i] *= s

    def Vec3Add(self,v1:list,v2:list):
        '''
        向量加法
        '''
        for i in range(len(v1)):
            v1[i]+=v2[i]

    def Vec3Sub(self,v1:list,v2:list):
        '''
        向量减法
        v1 - v2
        '''
        for i in range(len(v1)):
            v1[i]-=v2[i]

    def Vec3Scalar(self,v1:list,v2:list):
        '''
        向量点乘（内积）
        '''
        return v1[0]*v2[0] + v1[1]*v2[1] + v1[2]*v2[2]

    def Vec3Cross(self,v1:list,v2:list):
        '''
        向量叉乘（外积）
        '''
        vRes = [0,0,0]
        vRes[0] = v1[1]*v2[2] - v1[2]*v2[1]
        vRes[1] = v1[2]*v2[0] - v1[0]*v2[2]
        vRes[2] = v1[0]*v2[1] - v1[1]*v2[0]
        v1 = vRes.copy()

    def Mat3VecMult(self,v:list,m: list):
        '''
        向量和矩阵(3x3)的外积
        '''
        vRes = [0,0,0]
        vRes[0] = v[0]*m[0][0] + v[1]*m[1][0] + v[2]*m[2][0]
        vRes[1] = v[0]*m[0][1] + v[1]*m[1][1] + v[2]*m[2][1]
        vRes[2] = v[0]*m[0][2] + v[1]*m[1][2] + v[2]*m[2][2]
        v = vRes.copy()

    def Vec3ScaleAdd(self,v1:list,s:float,v2:list):
        '''
        向量伸缩加法
        '''
        for i in range(len(v1)):
            v1[i] += s*v2[i]

    def isNum(self,a):
        try:
            t = float(a)
            return True
        except:
            return False










