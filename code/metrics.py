# importing pandas module
import h5py
import pandas as pd
import numpy as np
def main():


    # reading csv file
    data = pd.DataFrame(np.array(h5py.File("C:\Users\rohan\PycharmProjects\SciFaitEmory\data\")['variable_1']))




if(__name__=="__main__"):
    main()