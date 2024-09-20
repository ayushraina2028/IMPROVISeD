import argparse
import numpy as np
import scipy.spatial as spatial

from Bio.PDB import *

'''
python makeCrosslinks.py -p 'prtA.pdb,prtB.pdb,prtC.pdb,prtD.pdb' -b 'prtA,prtB,prtC,prtD'

python makeCrosslinks.py -p 'prtA.pdb,prtB.pdb,prtC.pdb,prtD.pdb' -b 'prtA,prtB,prtC,prtD' -d 25

python makeCrosslinks.py -p 'prtA.pdb,prtB.pdb,prtC.pdb,prtD.pdb' -b 'prtA,prtB,prtC,prtD' -d 20

python makeCrosslinks.py -p 'prtA.pdb,prtB.pdb,prtC.pdb,prtD.pdb' -b 'prtA,prtB,prtC,prtD' -d 15

python makeCrosslinks.py -p 'prtA.pdb,prtB.pdb,prtC.pdb,prtD.pdb' -b 'prtA,prtB,prtC,prtD' -d 10
'''


def giveCoords(residue, resname, resinum, atomname='CA'):
    '''

    '''
    res = []
    for atoms in residue:
        if atoms.get_name() == atomname:
             res = [resinum, resname, atoms.get_name(), atoms.get_serial_number(), atoms.get_coord()]

    return res

def getCalphaForResidues(pdbfile, obj1, residuelist=[]):
    '''

    '''

    parse = PDBParser()
    st = parse.get_structure(obj1, pdbfile)
    
    res = []
    for model in st:
        for chain in model:
            for residue in chain:
                if residue.has_id("CA"):
                    #print(residue.segid)
                    resi_id = residue.id
                    res_num = resi_id[1]
                    if residuelist:
                        if residue.get_resname() in residuelist:
                            res.append(giveCoords(residue, residue.get_resname(), res_num))
                    else:
                        res.append(giveCoords(residue, residue.get_resname(), res_num))

    return res

def giveNpArrayNDict(list_of_coords):
    '''

    '''

    Y = np.empty((0,3), float)
    Y_dict = {} 

    count = 0
    for items in list_of_coords:
        Y = np.insert(Y, count, items[4], axis=0)
        Y_dict[count] = (items[0], items[1], items[2], items[3])
        count = count+1

    return Y, Y_dict
                        

def giveCrosslinksPair(ref_coords, ref_dict, ref_name, targ_coords, targ_dict, targ_name, dist_radius, opfilename=''):
    '''

    '''
    
    point_tree = spatial.cKDTree(targ_coords)

    if not opfilename:
        opfilename = ref_name+'_crosslink_'+targ_name+'.csv'

    print(opfilename)
    count = 0
    with open(opfilename,'w') as op:
         op.write('res1,prot1,res2,prot2')
         for ref_item in ref_coords:
             pts_near = point_tree.query_ball_point(ref_item, dist_radius)
             ref_item = ref_dict.get(count)
             for i in pts_near:
                 targ_item = targ_dict.get(i)
                 #print('\n%d,%s,%d,%s'%(ref_item[0],ref_name,targ_item[0],targ_name))
                 op.write('\n%d,%s,%d,%s'%(ref_item[0],ref_name,targ_item[0],targ_name))

             count = count +1

def driverCode(list_of_pdbs, list_of_pdbObjs, residuelist=['LYS'], dist_cutoff=30):
    '''

    '''
    
    resis = "_".join([i.upper() for i in residuelist])

    for i in range(0,len(list_of_pdbs)):
        ref_resi_ca       = getCalphaForResidues(list_of_pdbs[i], list_of_pdbObjs[i], residuelist)
        ref_Y, ref_Y_dict = giveNpArrayNDict(ref_resi_ca)

        for j in range(i+1, len(list_of_pdbs)):
            print(list_of_pdbs[i]+' vs '+list_of_pdbs[j])

            targ_resi_ca        = getCalphaForResidues(list_of_pdbs[j], list_of_pdbObjs[j], residuelist)
            targ_Y, targ_Y_dict = giveNpArrayNDict(targ_resi_ca)

            oppfilename = '%s_crosslink%d_%s_%s_Ca.csv'%(list_of_pdbObjs[i],dist_cutoff,list_of_pdbObjs[j],resis)
            giveCrosslinksPair(ref_Y, ref_Y_dict, list_of_pdbObjs[i], targ_Y, targ_Y_dict, list_of_pdbObjs[j], dist_cutoff, oppfilename)

def parse_arguments():
    parser = argparse.ArgumentParser()
    
    parser.add_argument('-p', '--pdb',
                        default=None,
                        help="Enter pdb names separated by comma")

    parser.add_argument('-b', '--obj',
                        default=None,
                        help="Enter object names separated by comma")

    parser.add_argument('-r', '--resi',
                        default='LYS',
                        help="Enter 3 letter residue name separated by comma")

    parser.add_argument('-d', '--dist',
                        default=30,
                        help="Enter crosslinks distance")

    return parser.parse_args()


if __name__ == '__main__':

    args = parse_arguments()

    pdbs = args.pdb.split(',')
    objs = args.obj.split(',')
    resi = args.resi.split(',')
    dist = float(args.dist)

    if len(pdbs) != len(objs):
       print('No of pdbs and objs must be same')
    else:
       driverCode(pdbs, objs, resi, dist)
