import json
import scipy
import pandas as pd
import numpy as np
from random import SystemRandom
from pymol import cmd


'''
createEqnCross(['prtA','prtB','prtC','prtD'],'prtA_crosslink30_prtB_LYS_Ca.csv,prtA_crosslink30_prtC_LYS_Ca.csv,prtA_crosslink30_prtD_LYS_Ca.csv,prtB_crosslink30_prtC_LYS_Ca.csv,prtB_crosslink30_prtD_LYS_Ca.csv,prtC_crosslink30_prtD_LYS_Ca.csv')
'''

def convertToDf(crosslinkfilelist):
    '''

    '''

    df = pd.DataFrame(columns=['res1','prot1','res2','prot2'])
    crosslinklist = crosslinkfilelist.split(',')
    for files in crosslinklist:
        tmp = pd.read_csv(files)
        df = df.append(tmp, ignore_index=True)

    return df

def giveDictsFromCrosslinks(df_crosslink):
    '''
    the function assumes that 2nd and 4th col of crosslinkfile
    is the name of the dictionary
    '''

    unique_prots = list(set(df_crosslink['prot1'].to_list()).union(set(df_crosslink['prot2'].to_list())))

    dict_resi = {}
    count = 0
    for i in unique_prots:
        
        # create dictionary
        row_prot1 = df_crosslink.loc[df_crosslink['prot1']==i]
        row_prot2 = df_crosslink.loc[df_crosslink['prot2']==i]

        resis = sorted(set(row_prot1['res1'].to_list()).union(set(row_prot2['res2'].to_list())))

        dict_resi[i] = resis

    return dict_resi

def getEqDict(dict_resi):
    '''
    assumes that the pdbnames are dict_names.pdb
    '''

    # load 
    keys = dict_resi.keys()

    # eq_dict
    eq_dict = {}

    for item in keys:
        #cmd.load('%s.pdb'%(item), item)

        np_arr = np.empty((0,3), float)
        count = 0
        for resi in dict_resi.get(item):
            query  = '%s///%s/CA'%(item,resi)
            np_arr = np.insert( np_arr, count, cmd.get_coords(query)[0], axis=0)
            count  = count+1

        eq_dict[item] = np_arr

    return eq_dict

def createIndexForObjs(objectslist):
    '''

    '''
    count = 0
    dict_objs_start = {}
    for obj in objectslist:
        if count == 0:
            dict_objs_start[obj] = 0
        else:
            dict_objs_start[obj] = count
        count = count + cmd.count_atoms(obj)
    return dict_objs_start

def getIndxForAtoms(objname, dict_objs_start, residue_no, atom_name = ['CA']):
    '''

    '''
    
    query = '(%s and resi %d)'%(objname, residue_no)
    myspace = {'atom_nm'  :  [], 'atom_no': []}
    cmd.iterate(query, 'atom_nm.append(name), atom_no.append(index)', space=myspace, quiet=1)

    ret_atm_nm  = []
    ret_atm_pos = []
    for i in range(0,len(myspace.get('atom_nm'))):
        if myspace.get('atom_nm')[i] in atom_name:            
           ret_atm_nm.append(myspace.get('atom_nm')[i])

           offset_obj = dict_objs_start[objname]
           ret_atm_pos.append(offset_obj+myspace.get('atom_no')[i])

    #print(ret_atm_nm)
    #print(ret_atm_pos)

    return ret_atm_nm, ret_atm_pos

def getIndxForAtom(objname, dict_objs_start, residue_no, atom_name = 'CA'):
    '''

    '''
    
    query = '(%s and resi %d)'%(objname, residue_no)
    myspace = {'atom_nm'  :  [], 'atom_no': []}
    cmd.iterate(query, 'atom_nm.append(name), atom_no.append(index)', space=myspace, quiet=1)

    ret_atm_nm  = []
    ret_atm_pos = []
    for i in range(0,len(myspace.get('atom_nm'))):
        if myspace.get('atom_nm')[i] == atom_name:            
           ret_atm_nm.append(myspace.get('atom_nm')[i])

           offset_obj = dict_objs_start[objname]
           ret_atm_pos.append(offset_obj+myspace.get('atom_no')[i])

    #print(ret_atm_nm)
    #print(ret_atm_pos)

    return ret_atm_pos[0]



def crosslinksIndexed(df_crosslinks, dict_objs_start):
    '''

    '''
    df_op = pd.DataFrame(columns = ['prot1','res1','indx1','prot2','res2','indx2'])

    for indx, row in df_crosslinks.iterrows():
        indx_res1 = getIndxForAtom(row['prot1'], dict_objs_start, row['res1'])
        indx_res2 = getIndxForAtom(row['prot2'], dict_objs_start, row['res2'])

        new_row = pd.DataFrame([{'prot1':row['prot1'], 'res1':row['res1'] ,'indx1': indx_res1, 'prot2':row['prot2'], 'res2':row['res2'], 'indx2': indx_res2}])
        df_op = pd.concat([df_op, new_row], ignore_index=True)

    return df_op

def equalityIndexed(dict_cross, dict_eqs, dict_objs_start):
    '''

    '''

    dict_eq_indexed = {}
    for keys in dict_cross.keys():
        resi_nos         = dict_cross.get(keys)
        resi_nos_indexed = [ getIndxForAtom(keys, dict_objs_start, i) for i in resi_nos ]

        np_arr = dict_eqs.get(keys)
        dist = scipy.spatial.distance.cdist(np_arr, np_arr)

        tmp_list_eq_indexed = []
        for i in range(0,dist.shape[0]):
            for j in range(i+1, dist.shape[1]):
                  tmp = [ keys, resi_nos[i], resi_nos_indexed[i], keys, resi_nos[j], resi_nos_indexed[j], dist[i,j]]
                  tmp_list_eq_indexed.append(tmp)

        dict_eq_indexed[keys] = tmp_list_eq_indexed

    return dict_eq_indexed


def df_equalityIndexed(dict_cross, dict_eqs, dict_objs_start):
    '''

    '''

    dict_eq_indexed = {}
    df = pd.DataFrame(columns=['prot1','globalindx1','tmpindx1','prot2','globalindx2','tmpindx2', 'dist'])

    for keys in dict_cross.keys():
        resi_nos         = dict_cross.get(keys)
        resi_nos_indexed = [ getIndxForAtom(keys, dict_objs_start, i) for i in resi_nos ]

        np_arr = dict_eqs.get(keys)
        dist = scipy.spatial.distance.cdist(np_arr, np_arr)

        print('-------------')
        print(keys)
        print(np_arr)
        print(dist)

        for i in range(0,dist.shape[0]):
            for j in range(i+1, dist.shape[1]):
                  newrow = pd.DataFrame([{'prot1':keys, 'globalindx1':resi_nos_indexed[i], 'tmpindx1':resi_nos_indexed[i]-dict_objs_start.get(keys),
                                          'prot2':keys, 'globalindx2':resi_nos_indexed[j], 'tmpindx2':resi_nos_indexed[j]-dict_objs_start.get(keys),
                                          'dist':dist[i,j]}])

                  df = pd.concat([df, newrow], ignore_index=True)


    return df


def df_equalityIndexed_onlyresi(dict_cross, dict_eqs, dict_objs_start):
    '''

    '''
    
    dict_eq_indexed = {}
    #df = pd.DataFrame(columns=['prot1','globalindx1','tmpindx1','prot2','globalindx2','tmpindx2', 'dist'])
    df = pd.DataFrame(columns=['prot1','resi1','prot2','resi2', 'dist'])


    for keys in dict_cross.keys():
        resi_nos         = dict_cross.get(keys)
        resi_nos_indexed = [ getIndxForAtom(keys, dict_objs_start, i) for i in resi_nos ]

        np_arr = dict_eqs.get(keys)
        dist = scipy.spatial.distance.cdist(np_arr, np_arr)
    
        print('-------------')
        print(keys)
        print(np_arr)
        print(dist)

        for i in range(0,dist.shape[0]):
            for j in range(i+1, dist.shape[1]):
                  #newrow = pd.DataFrame([{'prot1':keys, 'globalindx1':resi_nos_indexed[i], 'tmpindx1':resi_nos_indexed[i]-dict_objs_start.get(keys),
                  #                        'prot2':keys, 'globalindx2':resi_nos_indexed[j], 'tmpindx2':resi_nos_indexed[j]-dict_objs_start.get(keys),
                  #                        'dist':dist[i,j]}])
                  newrow = pd.DataFrame([{'prot1':keys, 'resi1':resi_nos[i],
                                          'prot2':keys, 'resi2':resi_nos[j],
                                          'dist':dist[i,j]}])

                  df = pd.concat([df, newrow], ignore_index=True)


    return df


def mapPointsToIndex(dict_cross, dict_objs_start):
    '''

    '''
    dict_points_indexed = {}
    for key in dict_cross.keys():
        tmp = []
        for points in dict_cross.get(key):
            tmp.append(getIndxForAtom(key, dict_objs_start, points))
        dict_points_indexed[key] = tmp

    return dict_points_indexed

def createXnIndexProts(dict_pts_indexed, dict_cross, dict_objs_start):
    '''

    '''

    df = pd.DataFrame(columns = ['prot','globalindx','tmpindx','x','y','z'])
    for key in dict_pts_indexed.keys():
        global_indx   = dict_pts_indexed.get(key)
        tmp_indx_resi = dict_cross.get(key)
        for i in range(0,len(global_indx)):
            glo_ndx = global_indx[i]
            #tmp_ndx = tmp_indx_resi[i]
            tmp_ndx = glo_ndx - dict_objs_start.get(key)
            query   = '%s///%d/CA'%(key,tmp_indx_resi[i])
            nparr   = cmd.get_coords(query)[0]
            new_row = pd.DataFrame([{'prot':key, 'globalindx':glo_ndx, 'tmpindx':tmp_ndx,  'x':nparr[0], 'y':nparr[1], 'z':nparr[2]}])
            df      = pd.concat([df, new_row], ignore_index=True)

    return df

def createXnIndexCross(df_x_n_indx_prots, dict_objs_start):
    '''

    '''
    df = pd.DataFrame(columns=['prot','globalindx','tmpindx','x','y','z'])

    for indx, row in df_x_n_indx_prots.iterrows():
        key1 = row['prot1']
        key2 = row['prot2']

        gl_indx1 = row['indx1']
        gl_indx2 = row['indx2']

        res1 = row['res1']
        res2 = row['res2']
        tmp1_indx = gl_indx1 - dict_objs_start.get(key1)
        tmp2_indx = gl_indx2 - dict_objs_start.get(key2)

        query1 = '%s///%d/CA'%(key1,res1)
        query2 = '%s///%d/CA'%(key2,res2)

        coord1 = cmd.get_coords(query1)[0]
        coord2 = cmd.get_coords(query2)[0]

        r1 = pd.DataFrame([{'prot':key1, 'globalindx':gl_indx1, 'tmpindx':tmp1_indx, 'x':coord1[0], 'y':coord1[1], 'z':coord1[2]}])
        r2 = pd.DataFrame([{'prot':key2, 'globalindx':gl_indx2, 'tmpindx':tmp2_indx, 'x':coord2[0], 'y':coord2[1], 'z':coord2[2]}])

        df = pd.concat([df, r1], ignore_index=True)
        df = pd.concat([df, r2], ignore_index=True)

    return df.drop_duplicates()


def creatCrosslinkDist(df_cross_indxed, dict_objs_start):
    '''

    '''
    df = pd.DataFrame(columns=['prot1','globalindx1','tmpindx1','prot2','globalindx2','tmpindx2', 'dist'])
    for indx, row in df_cross_indxed.iterrows():
        key1 = row['prot1']
        key2 = row['prot2']
     
        gl_indx1 = row['indx1']
        gl_indx2 = row['indx2']

        res1 = row['res1']
        res2 = row['res2']
        tmpindx1 = gl_indx1 - dict_objs_start.get(key1)
        tmpindx2 = gl_indx2 - dict_objs_start.get(key2)

        query1 = '%s///%d/CA'%(key1,res1)
        query2 = '%s///%d/CA'%(key2,res2)
 
        coord1 = cmd.get_coords(query1)[0]
        coord2 = cmd.get_coords(query2)[0]

        temp = coord1-coord2
        dist = np.sqrt(np.sum(np.square(temp)))

        new_row = pd.DataFrame([{'prot1':key1, 'globalindx1':gl_indx1, 'tmpindx1':tmpindx1,
                                 'prot2':key2, 'globalindx2':gl_indx2, 'tmpindx2':tmpindx2,
                                 'dist':dist}])
        df = pd.concat([df, new_row], ignore_index=True)

    return df

def df_indexAllObjs(objs, dict_objs_start):
    '''

    '''
    df = pd.DataFrame(columns=['prot', 'globalindx','tmpindx','x','y','z'])
    for obj in objs:
        myspace = { 'prot': [], 'globalindx': [], 'tmpindx': [], 'x_': [], 'y_': [], 'z_': [] }

        cmd.iterate_state(1, obj, 'prot.append(model), globalindx.append(index), tmpindx.append(resv), x_.append(x), y_.append(y), z_.append(z)', space=myspace)

        offset = dict_objs_start.get(obj)
        for i in range(0, len(myspace.get('prot'))):
            new_row = pd.DataFrame([{'prot': myspace.get('prot')[i],
                                     'globalindx': myspace.get('globalindx')[i] + offset,
                                     'tmpindx': myspace.get('globalindx')[i],
                                     'x': myspace.get('x_')[i], 'y': myspace.get('y_')[i], 'z':myspace.get('z_')[i]}])
            df = pd.concat([df, new_row], ignore_index=True)
    return df

def createEqnCross(objects, crosslinkfiles, op_suffix):
    '''

    '''
    
    # -- get all the crosslinks as dataframe file 
    df_cross   = convertToDf(crosslinkfiles)

    # -- get the unions of the points each pdbwise from crosslinks
    dict_cross = giveDictsFromCrosslinks(df_cross)

    # -- find the xyz coords corr to the points from crosslinks pdbwise
    dict_eqs   = getEqDict(dict_cross)

    # -- assign a global indexing for all the pdbs combined together
    dict_objs_start = createIndexForObjs(objects)


    # -- index the croslinks points based on global indexing
    df_cross_indxed = crosslinksIndexed(df_cross, dict_objs_start)

    # -- index the points from each pdb which shares a 
    dict_points_indexed = mapPointsToIndex(dict_cross, dict_objs_start)

    ####################################### Data for registration ##############################################
    # -- local PDB coordinates with local and global indexing
    #df_x_n_indx_prots = createXnIndexProts(dict_points_indexed, dict_cross, dict_objs_start)                       # use for calpha only
    #filename = 'calpha_prots_%s.csv'%(op_suffix)
    #df_x_n_indx_prots.to_csv(filename)

    # -- local points from crosslinks with local and global indexing (dist field only for taylor made examples)
    #df_x_n_indx_cross = createXnIndexCross(df_cross_indxed, dict_objs_start)
    #filename = 'crosslink_atoms_%s.csv'%(op_suffix)
    #df_x_n_indx_cross.to_csv(filename)

    # -- local points from each pdb with global index 
    #df_x_n_index_allobjs = df_indexAllObjs(objects, dict_objs_start)                                               # use for entire protin
    #filename = 'allatoms_prots_%s.csv'%(op_suffix)
    #df_x_n_index_allobjs.to_csv(filename)

    ####################################### Data for localization ##############################################
    # -- distances for the crosslinks with global and local indexing (distance mentioned for tailor made eg for verfication only)
    #df_cross_links    = creatCrosslinkDist(df_cross_indxed, dict_objs_start)
    #filename = 'cross_links_pairs_%s.csv'%(op_suffix)
    #df_cross_links.to_csv(filename)

    # -- distance between atoms of the same protein with global and local indexing
    #df_eq_bounds      = df_equalityIndexed(dict_cross, dict_eqs, dict_objs_start)
    df_eq_bounds     = df_equalityIndexed_onlyresi(dict_cross, dict_eqs, dict_objs_start)

    filename = 'eq_dists_prots_%s.csv'%(op_suffix)
    df_eq_bounds.to_csv(filename, index=False)

cmd.extend("createEqnCross",createEqnCross)
