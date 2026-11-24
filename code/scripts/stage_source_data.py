"""Copy preserved study outputs into the paths expected by the analysis scripts."""
import argparse
import hashlib
import json
import shutil
from pathlib import Path

parser=argparse.ArgumentParser(description=__doc__)
parser.add_argument('--source',type=Path,required=True,help='Path to the source data folder.')
parser.add_argument('--code-root',type=Path,default=Path.cwd())
args=parser.parse_args()
source=args.source.resolve()
root=args.code_root.resolve()
records=json.loads((source/'manifest.json').read_text())
jobs=[]
for rec in records:
    if not rec.get('repository_path'):
        continue
    src=(source/rec['path']).resolve()
    dst=(root/rec['repository_path']).resolve()
    if not src.is_relative_to(source) or not dst.is_relative_to(root):
        raise ValueError('A manifest path lies outside its root directory.')
    digest=hashlib.sha256(src.read_bytes()).hexdigest()
    if digest!=rec['sha256']:
        raise ValueError(f'Source checksum mismatch: {rec["path"]}')
    if dst.exists() and hashlib.sha256(dst.read_bytes()).hexdigest()!=digest:
        raise FileExistsError(f'Existing output differs: {dst}. Use a clean working copy.')
    jobs.append((src,dst))
for src,dst in jobs:
    dst.parent.mkdir(parents=True,exist_ok=True)
    if not dst.exists():
        shutil.copy2(src,dst)
print(f'Staged or verified {len(jobs)} source files.')
