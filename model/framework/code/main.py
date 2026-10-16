# imports
import sys

from ersilia_pack_utils.core import read_smiles, write_out

from sampler import AnalogSampler, N_OUTPUTS

# parse arguments
input_file = sys.argv[1]
output_file = sys.argv[2]

# read SMILES from .csv file, assuming one column with header
_, smiles_list = read_smiles(input_file)

# load the model once (its checkpoint is fetched from eosvc into model/checkpoints/),
# then sample scaffold-hopping analogues for every input
sampler = AnalogSampler()
outputs = []
for smi in smiles_list:
    hops = sampler.hops(smi)
    outputs.append(hops + [""] * (N_OUTPUTS - len(hops)))  # empty cells where fewer than 100 were found

# check input and output have the same length
assert len(smiles_list) == len(outputs)

header = [f"smi_{str(i).zfill(2)}" for i in range(N_OUTPUTS)]

# write output in a .csv file
write_out(outputs, header, output_file, str)
