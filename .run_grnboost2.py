#!/usr/bin/env python
# -*- coding=utf8 -*-
"""
#Author:QianJun
#Created Time:Thu 27 Jan 2022 09:00:53 PM CST
#File Name:run_grnboost2.py
#Description:to run grnboost2
"""

import pandas as pd
from arboreto.utils import load_tf_names
from arboreto.algo import grnboost2
from distributed import LocalCluster, Client
ex_path = "/home/qj/Rproject/scCRISPR/int/1.1_exprMatrix_filtered_t.txt"
tf_path = "/home/qj/Rproject/scCRISPR/int/1.1_inputTFs.txt"
 
if __name__ == '__main__':
    # create custom LocalCluster and Client instances
    local_cluster = LocalCluster(n_workers=2,threads_per_worker=1)
    custom_client = Client(local_cluster)
 
    # load the data
    ex_matrix = pd.read_csv(ex_path, sep='\t')
    tf_names = load_tf_names(tf_path)

    # run GRN inference multiple times
    network = grnboost2(expression_data=ex_matrix,tf_names=tf_names,client_or_address=custom_client,seed=313)
    network.to_csv('/home/qj/Rproject/scCRISPR/int/output_grnboost.tsv', sep='\t', index=False, header=False)
    # close the Client and LocalCluster after use
    custom_client.close()
    local_cluster.close()

