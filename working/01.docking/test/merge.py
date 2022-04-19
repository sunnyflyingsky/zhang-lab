import os
import sys
import getopt

Usage = 'Usage: ' + sys.argv[0]
Usage += '''
    <Requires>
    -i input pdb file path
    -r ref pdb path
    -o output path

    [Options]
    None
'''
Usage += "EX: " + sys.argv[0] + ''' -i /Share2/home/zhangqf7/yuhan/rotation_student/zhangjiasheng/Docking/RNA/data_set/xxxx.pdb 
    -r /Share2/home/zhangqf7/yuhan/rotation_student/zhangjiasheng/Docking/RNA/data_set/xxxx.pdb 
    -o /Share2/home/zhangqf7/yuhan/rotation_student/zhangjiasheng/Docking/RNA/temp_set/xxxx.pdb'''
if len(sys.argv)<4 or not sys.argv[1].startswith('-'):sys.exit(Usage)

if __name__=='__main__':  
    print('run!')
    oplist,alist = getopt.getopt(sys.argv[1:],'hi:r:o:')
    for opt in oplist:
        if opt[0] == '-h':sys.exit(Usage)
        elif opt[0] == '-i':inpath = opt[1]
        elif opt[0] == '-r':refpath = opt[1]
        elif opt[0] == '-o':outpath = opt[1]
        else:
            sys.exit(Usage)

    #inpath = 'docking/test/4_results_final_noH.pdb' 
    #refpath = 'docking/test/4_results_tmp.pdb'
    #utpath = 'docking/test/4_results_final_merged.pdb'
    try:
        read_object = open(refpath)
    except:
        sys.exit('ref file path is not exit!')
    chain = {}
    for line in read_object:
        if line.startswith('ATOM'):
            seq = line[21:22]
            X = float(line[30:38].strip())
            Y = float(line[38:46].strip())
            Z = float(line[46:54].strip())
            chain[(X,Y,Z)] = seq
    read_object.close()
    try:
        write_object = open(outpath,'w')
    except:
        sys.exit('output folder path is not exit!')
    try:
        read_object = open(inpath)
    except:
        sys.exit('input file path is not exit!')
    for line in read_object:
        if line.startswith('TER'):
            line_merged = line[:21]+seq+line[22:]
            write_object.write(line_merged)
        elif line.startswith('ATOM'):
            X = float(line[30:38].strip())
            Y = float(line[38:46].strip())
            Z = float(line[46:54].strip())
            seq = chain[(X,Y,Z)]
            line_merged = line[:21]+seq+line[22:]
            write_object.write(line_merged)
        else:
            write_object.write(line)
    read_object.close()
    write_object.close()
    print('end!')

