#!/bin/python

import glob
import os
import re
import numpy as np


dirName_arr = ["/proj/steinlab/HTSF/220311_UNC41-A00434_0442_BHMNNTDSX2", "/proj/steinlab/HTSF/221014_UNC41-A00434_0550_BHNG7WDSX3"]


with open("jobList_alevin_fry_WTK5to6.txt", "w") as fw:   # order of project, R1, R2, out_dir, newYamalFileName
	for dirName in dirName_arr:
		os.chdir(dirName)
		fastq_arr = glob.glob("*.fastq.gz")

		projectName_arr = []
		for fastqFileName in fastq_arr:
			projectName_arr.extend([ re.split('_R[1-2]_', fastqFileName)[0] ])

		projectName_arr = np.unique(projectName_arr)

		for projectName in projectName_arr:
			line = [projectName]
			R1 = dirName + "/" + glob.glob(projectName + "_R1_*.fastq.gz")[0]
			R2 = dirName + "/" + glob.glob(projectName + "_R2_*.fastq.gz")[0]
			line.extend([R1, R2])

			fw.write('\t'.join([x for x in line]) + "\n")










