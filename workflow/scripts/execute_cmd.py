

import os
import sys
import subprocess

#with open("/home/hemit101/testing/test2.txt", "a") as file:

env = os.environ["CONDA_PREFIX"]
cmd = snakemake.params.cmd

    # file.write(str(cmd))
    # file.write("\n")
    # file.write(str(env))
    # file.write("\n")

cmd = cmd.replace("CONDA_PREFIX", env)

    # file.write(str(cmd))
    # file.close()



with open("/home/hemit101/testing/test2.txt", "a") as file:
    process = subprocess.run(cmd, shell=True)

    file.write(str(process.returncode))
    file.write("\n")
    file.write(str(process.stdout))
    file.write("\n")
    file.write(str(process.stderr))
    file.write("\n")