import os
import sys

def GetCallCount(vcf_name):
    if not os.path.isfile(vcf_name): return 0
    with open(vcf_name) as file:
        content = file.readlines()
    vcf_header = [x for x in content if x[0] == "#"]
    return len(content) - len(vcf_header)
