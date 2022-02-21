run_grnboost2 <- function(output, cores, seed){
  header <- "#!/usr/bin/env python\n# -*- coding=utf8 -*-\n\"\"\"\n#Author:QianJun\n#Created Time:Thu 27 Jan 2022 09:00:53 PM CST\n#File Name:run_grnboost2.py\n#Description:to run grnboost2\n\"\"\"\n
import pandas as pd\nfrom arboreto.utils import load_tf_names\nfrom arboreto.algo import grnboost2\nfrom distributed import LocalCluster, Client\n"
  ex_path = paste("ex_path = \"", output, "/int/1.1_exprMatrix_filtered_t.txt\"\n", sep ="")
  tf_path = paste("tf_path = \"", output, "/int/1.1_inputTFs.txt\"\n \n", sep = "")

  mid1 <- "if __name__ == '__main__':\n    # create custom LocalCluster and Client instances\n"
  core <- paste("    local_cluster = LocalCluster(n_workers=",cores,",threads_per_worker=1)\n    custom_client = Client(local_cluster)\n \n",sep = "")

  mid2 <- paste("    # load the data\n    ex_matrix = pd.read_csv(ex_path, sep='\\t')\n    tf_names = load_tf_names(tf_path)\n
    # run GRN inference multiple times\n    network = grnboost2(expression_data=ex_matrix,tf_names=tf_names,client_or_address=custom_client,seed=",seed,")\n",sep="")

  savefile <- paste("    network.to_csv", "('", output, "/int/output_grnboost.tsv', sep='\\t', index=False, header=False)\n", sep = "")

  tail <- "    # close the Client and LocalCluster after use\n    custom_client.close()\n    local_cluster.close()\n"

  txt <- paste(header,ex_path,tf_path,mid1,core,mid2,savefile,tail,sep = "")
  writeLines(txt, ".run_grnboost2.py")
}

