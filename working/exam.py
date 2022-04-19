inputInfo = '/Share2/home/zhangqf7/yuhan/rotation_student/zhangjiasheng/Docking/RNA/data_set/xxxx.pdb '
info = inputInfo.strip().split('/')
inputPath = '/'.join(info[:-1])
name = info[-1][:-4]

print(inputPath,name)