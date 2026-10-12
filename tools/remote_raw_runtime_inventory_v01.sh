set -euo pipefail
python3 - <<'PY'
from pathlib import Path
root=Path('/vast/projects/Immune_Aging/Itokawa_et_al_2022_RNAseq/trimmed')
for p in root.glob('*HSC*cutadapt.log'):
 print('CUTADAPT_LOG',p)
 with p.open(errors='replace') as f:
  for _,line in zip(range(4),f):print(line.rstrip())
root=Path('/home/users/allstaff/stacpoole.q')
for folder in ['Scripts','scripts','20210719_Scripts','20210719_scripts','bin']:
 p=root/folder
 if p.is_dir():
  for f in p.rglob('*presplit*'):print('HELPER_CANDIDATE',f)
PY
module load STAR/2.7.11b bowtie2/2.5.3 picard-tools/3.1.1 openjdk/21.0.2 MACS3/3.0.1 cutadapt
for name in STAR bowtie2 cutadapt macs3 java; do
 command -v "$name"
done
STAR --version
bowtie2 --version | head -2
cutadapt --version
macs3 --version
java -version
module show cutadapt MACS3/3.0.1 openjdk/21.0.2 picard-tools/3.1.1 2>&1
