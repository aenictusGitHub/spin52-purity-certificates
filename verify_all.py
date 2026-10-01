"""Run the complete exact proof checks; optional independent numerical checks."""
from pathlib import Path
import argparse,subprocess,sys,time

def main():
    if not __debug__:
        raise RuntimeError('Do not use Python -O: verification relies on assertions.')
    p=argparse.ArgumentParser(description=__doc__)
    p.add_argument('--quick',action='store_true',help='Skip only regeneration of the order-two cached coefficient matrices.')
    p.add_argument('--numerical-crosschecks',action='store_true',help='Also run direct complex-state and three-copy numerical cross-checks.')
    args=p.parse_args();root=Path(__file__).resolve().parent
    jobs=[['verify_optimizers.py'],['prove_purity2.py']+(['--quick'] if args.quick else []),
          ['prove_purity3.py'],['prove_purity4_sixth.py'],['verify_roofs.py']]
    if args.numerical_crosschecks:jobs += [['verify_purity2_independently.py'],['verify_purity3_independently.py']]
    start=time.monotonic()
    print('Mode:', 'QUICK (cached order-two matrices)' if args.quick else 'FULL EXACT COEFFICIENT REGENERATION',flush=True)
    for job in jobs:
        print('\nVERIFY',*job,flush=True)
        subprocess.run([sys.executable,str(root/job[0]),*job[1:]],cwd=root,check=True)
    print(f'\nAll requested certificate checks passed in {time.monotonic()-start:.1f} seconds.',flush=True)
if __name__=='__main__':main()
