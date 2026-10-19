"""Regenerate the manuscript figures. See README.md for the two calculations."""
from pathlib import Path
import runpy,sys,shutil

ROOT=Path(__file__).resolve().parent

def main():
    legacy=ROOT/'supplementary/legacy'
    runpy.run_path(str(legacy/'reproduce.py'),run_name='__main__')
    for name in ['phase_portrait','melnikov','bifurcation_lyapunov']:
        for ext in ['pdf','png']:
            shutil.copyfile(legacy/'figures'/f'{name}.{ext}',ROOT/'figures'/f'{name}.{ext}')
    updated=ROOT/'supplementary/updated'
    sys.path.insert(0,str(updated))
    import pendulo_figuras as calculations
    calculations.figure3(directory=updated)
    calculations.figure5(directory=updated)
    calculations.convergence_figure(directory=updated)
    for source,target in [('figura_3_seis_regimes','poincare'),
                          ('figura_5_cruzamentos_e_fractal','fractal_dimension')]:
        for ext in ['pdf','png']:
            shutil.copyfile(updated/f'{source}.{ext}',ROOT/'figures'/f'{target}.{ext}')

if __name__=='__main__':main()
