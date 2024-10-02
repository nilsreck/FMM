import sys
import os

for dim in snakemake.params.dimension_list:
    file_name = str(snakemake.params.resource_path)+"/dimension_"+str(dim)
    file = open(file_name, "w")
    file.close()

for i in range(len(snakemake.params.dimension_list)-1):
    file_name = str(snakemake.params.resource_path)+"/comparison_"+str(snakemake.params.dimension_list[i])+"-"+str(snakemake.params.dimension_list[i+1])
    file = open(file_name, "w")
    file.close()