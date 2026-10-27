"""Create or verify the tracked public-file SHA-256 manifest."""
from pathlib import Path
import argparse,hashlib,subprocess
ROOT=Path(__file__).resolve().parents[1];MANIFEST=ROOT/'checksums_public.sha256'
def digest(path):return hashlib.sha256(path.read_bytes()).hexdigest()
def main():
    ap=argparse.ArgumentParser();ap.add_argument('--write',action='store_true');args=ap.parse_args()
    if args.write:
        paths=subprocess.check_output(['git','-c','core.fsmonitor=false','ls-files','-z'],cwd=ROOT).decode().split('\0')
        lines=[f'{digest(ROOT/p)}  {p}' for p in sorted(paths) if p and p!=MANIFEST.name]
        MANIFEST.write_text('\n'.join(lines)+'\n',encoding='utf-8',newline='\n')
        print(f'Wrote {len(lines)} tracked-file checksums')
    else:
        lines=MANIFEST.read_text().splitlines()
        for line in lines:
            expected,path=line.split('  ',1);assert digest(ROOT/path)==expected,path
        print(f'PASS: {len(lines)} public-file checksums')
if __name__=='__main__':main()
