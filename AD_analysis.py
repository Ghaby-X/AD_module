import pandas as pd
from rdkit import Chem
from rdkit.Chem import AllChem

import utils

class AD:

    #Initialization of this class takes in three parameters
    # the base parameter is the base data that you want to compare other compounds to
    #base parameter either takes in an array of smiles or a pandas dataframe consisting of only the fingerprints you want to compare
    
    def __init__(self, base, threshold = 0.8, base_type = 'smiles'):
        self.base = base
        self.threshold = threshold
        self.base_type = base_type

        if(self.base_type == 'smiles'):
            #if your self.base is not an array, dont create class
            if(not isinstance(self.base, list)):
                raise TypeError("'base' should be an array for a base_type of 'smiles'")
            
            #if an array, let's compute descriptors and save as ExplicitBitVect
            fps = utils.compute_morganfps(self.base)
            self.base = fps
        
        if(self.base_type == "fingerprint"):
            #check if self.base is a pandas dataframe
            if(not isinstance(self.base, pd.DataFrame)):
                raise TypeError("'base' should be a pandas Dataframe for type fingerprint")


    #function to redefine the various parameters
    def setParams(self, base, threshold, base_type):
        if (base):
            self.base = base
        if(threshold):
            self.threshold = threshold
        if (base_type):
            self.base_type = base_type

    def base_similarity(self):
        #return the similarites of all compounds in your base
        for i in range(len(self.base)):
            for j in range(i, len(self.base)):
                pass



    def compare_similarity(test):
        #return the similarity of your test to base
        testbitVect = utils.compute_morganfps(test)

        pass

    def visualize(test):
        #return the a visual applicability domain of the test with respect to base
        pass