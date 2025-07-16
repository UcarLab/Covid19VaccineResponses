# To match 10x style, we transpose and rename the kite output

import os.path as path
from scipy import io as sp_io
import argparse

parser = argparse.ArgumentParser()
parser.add_argument('input_path', type=str, help='Input path')
parser.add_argument('output_path', type=str, help='Output path')
args = parser.parse_args()

kite_out=sp_io.mmread(args.input_path)
print(kite_out.shape)
kite_out_T = kite_out.transpose()
print(kite_out_T.shape)
sp_io.mmwrite(args.output_path, kite_out_T)
print(f'Read from {args.input_path}, saved at {args.output_path}.')