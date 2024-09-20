import numpy as np

from pymol import cmd

def combineIndices(index_files):
    '''

    '''
    list_files = index_files.split(',')

    indx_arr   = np.empty(0)
    indx_size  = np.empty(len(list_files),dtype=int)
    #for file in list_files:
    for i in range(0,len(list_files)):
        file = list_files[i]
        indx_file     = np.loadtxt(file)
        indx_arr      = np.append(indx_arr, indx_file, axis=0)
        indx_size[i]  = indx_file.shape[0]

    return indx_size, indx_arr #.sort() 

def showRegPDBs(objects, reg_xyz_file, index_files):
    '''
    objects: pdb object
    '''
    objects_ = objects.split(',')
    sizes, indices  = combineIndices(index_files)

    reg_xyz = np.loadtxt(reg_xyz_file, delimiter=',')
    
    count = 0
    indices_ = index_files.split(',')
    #for obj in objects:
    for i in range(0, len(objects_)):
        obj = objects_[i]        
        indx_begin = count        
        indx_end   = count + sizes[i]
        xyz = reg_xyz[indx_begin:indx_end,:]
        print(xyz.shape)
        cmd.load_coords(xyz, obj)

        count = count + sizes[i]

cmd.extend('showRegPDBs',showRegPDBs)
