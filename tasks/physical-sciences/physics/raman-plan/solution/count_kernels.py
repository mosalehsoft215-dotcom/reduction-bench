from pathlib import Path
import sys,os
sys.path.insert(0,os.environ.get("RAMAN_DATA","/app/data"))
from count_engine import comparisons as fourier_kernels
