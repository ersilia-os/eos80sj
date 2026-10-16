"""Scaffold-hopping analogue sampler built on the ANNalog Python API.

ANNalog is a transformer sequence-to-sequence model that turns one SMILES into analogues.
This module keeps only the analogues whose Bemis-Murcko scaffold differs from the input's
(the "hops"), draws more until 100 are found or the draw budget is spent, and returns a
uniform random sample of them. Nothing is ranked or truncated by ANNalog's score: it is a
length-biased log-probability, so keeping the best-scoring analogues would skew the output
towards small, conservative edits.
"""
import os
import random
import sys

from annalog.model_handler import SMILESModelHandler
from annalog.SMILES_generator import SMILESGenerator
from rdkit import Chem, RDLogger
from rdkit.Chem.Scaffolds import MurckoScaffold

RDLogger.DisableLog("rdApp.*")

ROOT = os.path.dirname(os.path.abspath(__file__))
# Lev_extended.pt and the vocabulary are fetched from eosvc into model/checkpoints/
CHECKPOINTS = os.path.join(ROOT, "..", "..", "checkpoints")

N_OUTPUTS = 100      # columns in run_columns.csv
FIRST_ROUND = 120    # draws in the first round
TOP_UP = 20          # draws in every later round
MAX_DRAWS = 500      # give up on an input after this many draws
TEMPERATURE = 1.2    # the authors' value

# Fresh entropy for every seed. generate_smiles seeds python, numpy and torch globally on each
# call, so reusing a seed would replay the same samples, and a fixed one would make the model
# deterministic. SystemRandom is also unaffected by that global seeding.
_rng = random.SystemRandom()


def flat_scaffold(mol):
    """Murcko scaffold SMILES without stereochemistry ('' for an acyclic molecule)."""
    return Chem.MolToSmiles(MurckoScaffold.GetScaffoldForMol(mol), isomericSmiles=False)


class AnalogSampler:
    def __init__(self):
        vocab = os.path.join(CHECKPOINTS, "stereo_experiment_vocab_ttf.pkl")
        model = os.path.join(CHECKPOINTS, "Lev_extended.pt")
        handler = SMILESModelHandler(vocab, vocab, model, device="cpu")
        self.generator = SMILESGenerator(handler)

    def hops(self, smiles):
        """Up to N_OUTPUTS analogues of `smiles` with a different scaffold, in random order."""
        mol = Chem.MolFromSmiles(smiles) if isinstance(smiles, str) else None
        if mol is None:
            return []
        source = Chem.MolToSmiles(mol, canonical=False)  # the form the annalog CLI feeds the model
        scaffold = flat_scaffold(mol)

        pool, draws, n = {}, 0, FIRST_ROUND
        while len(pool) < N_OUTPUTS and draws < MAX_DRAWS:
            n = min(n, MAX_DRAWS - draws)
            try:
                # the method must be 'sampling'; the generate_smiles docstring says 'sample' and raises
                generated = self.generator.generate_smiles(
                    source,
                    generation_number=n,
                    generation_method="sampling",
                    temperature=TEMPERATURE,
                    prefix=0,
                    filter_invalid=False,
                    seed=_rng.randrange(2**31),
                )
            except Exception as exc:  # one failing molecule must not sink the batch
                print(f"ANNalog failed on {smiles!r}: {exc}", file=sys.stderr)
                break
            for smi, _score in generated:
                m = Chem.MolFromSmiles(smi)
                # a differing scaffold also removes the input itself; an empty molecule would
                # otherwise count as a hop
                if m is None or not m.GetNumAtoms() or flat_scaffold(m) == scaffold:
                    continue
                pool.setdefault(Chem.MolToSmiles(m), None)  # canonical isomeric SMILES
            draws += n
            n = TOP_UP

        hits = list(pool)
        return _rng.sample(hits, min(N_OUTPUTS, len(hits)))
