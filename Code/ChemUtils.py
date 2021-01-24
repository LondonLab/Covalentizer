import os
import subprocess
def convert(in_file, out_file):
	input_ter = in_file.split('.')[-1]
        output_ter = out_file.split('.')[-1]
        subprocess.call(["/work/londonlab/software/openbabel-2.3.2/bin/obabel", "-i", input_ter, in_file, "-o", output_ter, "-O", out_file])
