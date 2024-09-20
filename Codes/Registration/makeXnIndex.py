import pandas as pd
import numpy as np

from pymol import cmd

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

def getCoords(obj, dict_obj_start, selection='all'):
    '''

    '''


    #myspace = { 'prot': [], 'globalindx': [], 'tmpindx': [], 'resi': [], 'atom': [], 'x_': [], 'y_': [], 'z_': [] }
    #cmd.iterate_state(1, obj, 'prot.append(model), globalindx.append(index), tmpindx.append(resv), resi.append(resv), atom.append(name), x_.append(x), y_.append(y), z_.append(z)', space=myspace)

    myspace = { 'prot': [], 'globalindx': [], 'tmpindx': [], 'resino': [], 'atom_nm': [], 'x_': [], 'y_': [], 'z_': [] }
    cmd.iterate_state(1, obj, 'prot.append(model), globalindx.append(index), tmpindx.append(resv), resino.append(resv), atom_nm.append(name), x_.append(x), y_.append(y), z_.append(z)', space=myspace)

    df = pd.DataFrame(columns=['prot', 'globalindx','tmpindx','resi','atom','x','y','z'])
    offset = dict_obj_start.get(obj)
    for i in range(0, len(myspace.get('prot'))):
            new_row = pd.DataFrame([{'prot': myspace.get('prot')[i],
                                     'globalindx': myspace.get('globalindx')[i] + offset,
                                     'tmpindx': myspace.get('globalindx')[i],
                                     'resi': myspace.get('resino')[i], 'atom': myspace.get('atom_nm')[i],
                                     'x': myspace.get('x_')[i], 'y': myspace.get('y_')[i], 'z':myspace.get('z_')[i]}])
            df = pd.concat([df, new_row], ignore_index=True)
    return df

def getGlobalIndices(df_obj, crs_lnk_resis, atomname='CA'):
    '''

    '''
    ret_list = []
    print(crs_lnk_resis)
    print(df_obj)
    for resi in crs_lnk_resis:
        df = df_obj[ (df_obj['resi']==resi) & (df_obj['atom']==atomname)]
        ret_list.append(df['globalindx'].tolist()[0])

    return ret_list

def convertToLists(crosslink_resis):
    '''
    list of txt files
    '''

    ret_list = []
    for file in crosslink_resis:
        fp_  = open(file,'r')
        data = fp_.read()
        tmp  = data.split('\n')
        tmp_ = list(filter(None, tmp))
        res = [eval(i) for i in tmp_]
        ret_list.append(res)
        fp_.close()

    return ret_list

def convertToNp(crosslink_xyz):
    '''

    '''

    return np.loadtxt(crosslink_xyz, delimiter=',')

def makeXnIndex(objects, crosslink_xyz, crosslink_indices, opp_prefix):
    '''
              objects: list of pdb objects
        crosslink_xyz: xyz coordinates
    crosslink_indices: ".txt" files of indices in crosslink, 
                       total number of cumulative lines in all .txt is
                       same as that total no. of lines in crosslink_xyz
    '''
    
    # -- assign a global indexing for all the pdbs combined together --#
    dict_objs_start = createIndexForObjs(objects)

    # -- convert to list of lists -- #
    crs_lnk_resis = convertToLists(crosslink_indices)
    crs_lnk_xyz   = convertToNp(crosslink_xyz)
    print('-----------')
    print(crs_lnk_xyz)

    # --  write XnIndex -- #
    op_crs_indx = '%scrosslink_indx.txt'%(opp_prefix)
    if os.path.exists(op_crs_indx): os.remove(op_crs_indx)
    fp_crs_indx = open(op_crs_indx,'a')

    op_crs_xyz  = '%scrosslink_xyz.txt'%(opp_prefix)
    if os.path.exists(op_crs_xyz): os.remove(op_crs_xyz)
    fp_crs_xyz  = open(op_crs_xyz, 'a')

    count = 0
    for i in range(0,len(objects)):
        obj = objects[i]
         
        df_obj   = getCoords(obj, dict_objs_start)
        indx_obj = df_obj['globalindx'].to_list()
        xyz_obj  = np.array(df_obj[['x','y','z']])

        # -- write XnIndex of object-i into a file -- #
        opfile_indx = '%s_%s_indx.txt'%(opp_prefix ,obj)
        opfile_xyz  = '%s_%s_xyz.txt'%(opp_prefix, obj)

        with open(opfile_indx,'w') as fp:
            for indx in indx_obj:
                fp.write('%d\n'%(indx))

        header = 'x y z'
        np.savetxt(opfile_xyz, xyz_obj, fmt="%f")

        # -- write indiex of crosslinks for object-i -- #
        crs_indx_start = count
        crs_indx_end   = count + len(crs_lnk_resis[i])

        #opfile_crs_xyz = 'crosslink_xyz_%s'%(obj)
        np.savetxt(fp_crs_xyz, crs_lnk_xyz[crs_indx_start:crs_indx_end,:], fmt='%f')

        #opfile_crs_indx = 'crosslink_indx_%s'%(obj)
        #with open(opfile_crs_indx,'w') as fp:
        indx_crosslink_obj = getGlobalIndices(df_obj, crs_lnk_resis[i])
        for indx in indx_crosslink_obj:
            fp_crs_indx.write('%d\n'%(indx))

        count = count + len(crs_lnk_resis[i])

    fp_crs_xyz.close()
    fp_crs_indx.close()

cmd.extend('makeXnIndex',makeXnIndex)
