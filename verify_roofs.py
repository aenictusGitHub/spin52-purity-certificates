"""Verify both exact convex roofs and their four-variable simplification."""
from pathlib import Path
import subprocess
import sys
import time

def main():
    if not __debug__:
        raise RuntimeError('Do not use Python -O: verification relies on assertions.')
    root=Path(__file__).resolve().parent/'convex_roofs'
    start=time.monotonic()
    for name in ['verify_support.py','verify_q1.py','verify_q2_algebraic.py','verify_q2_reduced.py']:
        print('\nVERIFY convex_roofs/'+name,flush=True)
        subprocess.run([sys.executable,str(root/name)],cwd=root,check=True)
    print(f'\nAll exact convex-roof checks passed in {time.monotonic()-start:.1f} seconds.',flush=True)

if __name__=='__main__':
    main()
