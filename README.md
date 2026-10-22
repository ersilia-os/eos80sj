# ANNalog Scaffold Hopping

ANNalog is a transformer sequence-to-sequence model trained on pairs of molecules tested in the same ChEMBL33 assay, which it treats as analogues. The Ersilia model proposes scaffold hops of an input compound, meaning analogues whose Bemis-Murcko scaffold differs from its own. It samples at the authors temperature of 1.2 and resamples to reach 100 hops where possible. Acyclic inputs have no scaffold, so any ring-containing output counts as a hop. Training excluded molecules above 600 Da, with over seven rotatable bonds, or with elements beyond C, N, O, S, F, Cl and Br.

This model was incorporated on 2026-10-14.Last packaged on 2026-10-16.

## Information
### Identifiers
- **Ersilia Identifier:** `eos80sj`
- **Slug:** `annalog`

### Domain
- **Task:** `Sampling`
- **Subtask:** `Generation`
- **Biomedical Area:** `Any`
- **Target Organism:** `Any`
- **Tags:** `Compound generation`, `ChEMBL`, `Chemical language model`

### Input
- **Input:** `Compound`
- **Input Dimension:** `1`

### Output
- **Output Dimension:** `100`
- **Output Consistency:** `Variable`
- **Interpretation:** Up to 100 sampled analogues whose Bemis-Murcko scaffold differs from the inputs; unused slots stay empty.

Below are the **Output Columns** of the model:
| Name | Type | Direction | Description |
|------|------|-----------|-------------|
| smi_00 | string |  | Scaffold-hopping analogue index 0 sampled with ANNalog in random order |
| smi_01 | string |  | Scaffold-hopping analogue index 1 sampled with ANNalog in random order |
| smi_02 | string |  | Scaffold-hopping analogue index 2 sampled with ANNalog in random order |
| smi_03 | string |  | Scaffold-hopping analogue index 3 sampled with ANNalog in random order |
| smi_04 | string |  | Scaffold-hopping analogue index 4 sampled with ANNalog in random order |
| smi_05 | string |  | Scaffold-hopping analogue index 5 sampled with ANNalog in random order |
| smi_06 | string |  | Scaffold-hopping analogue index 6 sampled with ANNalog in random order |
| smi_07 | string |  | Scaffold-hopping analogue index 7 sampled with ANNalog in random order |
| smi_08 | string |  | Scaffold-hopping analogue index 8 sampled with ANNalog in random order |
| smi_09 | string |  | Scaffold-hopping analogue index 9 sampled with ANNalog in random order |

_10 of 100 columns are shown_
### Source and Deployment
- **Source:** `Local`
- **Source Type:** `External`
- **DockerHub**: [https://hub.docker.com/r/ersiliaos/eos80sj](https://hub.docker.com/r/ersiliaos/eos80sj)
- **Docker Architecture:** `AMD64`, `ARM64`
- **S3 Storage**: [https://ersilia-models-zipped.s3.eu-central-1.amazonaws.com/eos80sj.zip](https://ersilia-models-zipped.s3.eu-central-1.amazonaws.com/eos80sj.zip)

### Resource Consumption
- **Model Size (Mb):** `16`
- **Environment Size (Mb):** `1418`
- **Image Size (Mb):** `1361.52`

**Computational Performance (seconds):**
- 10 inputs: `52.13`
- 100 inputs: `-1`
- 10000 inputs: `-1`

### References
- **Source Code**: [https://github.com/DVNecromancer/ANNalog](https://github.com/DVNecromancer/ANNalog)
- **Publication**: [https://doi.org/10.1186/s13321-026-01186-6](https://doi.org/10.1186/s13321-026-01186-6)
- **Publication Type:** `Peer reviewed`
- **Publication Year:** `2026`
- **Ersilia Contributor:** [arnaucoma24](https://github.com/arnaucoma24)

### License
This package is licensed under a [GPL-3.0](https://github.com/ersilia-os/ersilia/blob/master/LICENSE) license. The model contained within this package is licensed under a [MIT](LICENSE) license.

**Notice**: Ersilia grants access to models _as is_, directly from the original authors, please refer to the original code repository and/or publication if you use the model in your research.


## Use
To use this model locally, you need to have the [Ersilia CLI](https://github.com/ersilia-os/ersilia) installed.
The model can be **fetched** using the following command:
```bash
# fetch model from the Ersilia Model Hub
ersilia fetch eos80sj
```
Then, you can **serve**, **run** and **close** the model as follows:
```bash
# serve the model
ersilia serve eos80sj
# generate an example file
ersilia example -n 3 -f my_input.csv
# run the model
ersilia run -i my_input.csv -o my_output.csv
# close the model
ersilia close
```

## About Ersilia
The [Ersilia Open Source Initiative](https://ersilia.io) is a tech non-profit organization fueling sustainable research in the Global South.
Please [cite](https://github.com/ersilia-os/ersilia/blob/master/CITATION.cff) the Ersilia Model Hub if you've found this model to be useful. Always [let us know](https://github.com/ersilia-os/ersilia/issues) if you experience any issues while trying to run it.
If you want to contribute to our mission, consider [donating](https://www.ersilia.io/donate) to Ersilia!
