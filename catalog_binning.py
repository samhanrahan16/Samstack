import numpy as np


class Field_catalog:
    '''Takes a data table and parameters file and will have methods to split the catalog into mass bins'''

    def __init__(self,tbl,params,zkey,mkey,rkey,dkey):

        self.catalog = tbl
        self.zkey = zkey
        self.mkey = mkey
        self.rkey = rkey
        self.dkey = dkey
        self.m_nodes = params['bins']['m_nodes']


    def make_bin_names(self):

        bin_names = []
        bin_names_to_midpoints = {}
        nodes = self.m_nodes
        for index in range(len(nodes)-1):
            bin_name = 'masses:' + str(nodes[index]) + '-' + str(nodes[index+1])
            bin_names.append(bin_name)
            diff = (nodes[index+1] - nodes[index])/2
            mid = nodes[index] + diff
            bin_names_to_midpoints[bin_name] = mid

        return bin_names,bin_names_to_midpoints

    def split_catalog_by_mass(self):
        '''Splits the catalog into mass bins, returns a dictionary of ra and decs for each bin'''

        binned_ra_dec = {}
        self.bin_names,self.bin_names_to_midpoints = self.make_bin_names()
        ras = self.catalog[self.rkey]
        decs =  self.catalog[self.dkey]
        masses = self.catalog[self.mkey]

        for i in range(len(self.m_nodes)-1):
            bin_name = self.bin_names[i]
            bin_ras = []
            bin_decs = []
            for j in range(len(masses)):
                if masses[j] >= self.m_nodes[i] and masses[j] < self.m_nodes[i+1]:
                    ra = ras[j]
                    dec = decs[j]
                    bin_ras.append(ra)
                    bin_decs.append(dec)

            binned_ra_dec[bin_name] = (bin_ras,bin_decs)

        return binned_ra_dec
                    
                    
