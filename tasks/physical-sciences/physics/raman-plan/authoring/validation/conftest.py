from pathlib import Path
import sys,os
ROOT=Path(__file__).resolve().parents[2]
sys.path[:0]=[str(ROOT/'solution'),str(ROOT/'tests'),str(ROOT/'environment/data')]
os.environ['RAMAN_DATA']=str(ROOT/'environment/data')
os.environ['OPENBLAS_NUM_THREADS']='1'
