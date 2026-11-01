"""Explicit imports of existing physical simulator and joint estimator helpers."""
import os
for name in ('OPENBLAS_NUM_THREADS','OMP_NUM_THREADS','MKL_NUM_THREADS'):
    os.environ.setdefault(name,'1')
import sys
from pathlib import Path
HERE=Path(__file__).resolve().parent
PARENT=HERE.parent
for path in (HERE.parents[1],PARENT/'imm_dpc',PARENT/'imm_ci',PARENT/'raw_sensor_maritime',PARENT/'joint_gp',HERE):
    sys.path.insert(0,str(path))
# Existing joint_gp.run_joint imports its historical reference as run_full.
# Bind that name explicitly so the new runner's identical basename cannot
# shadow the raw-sensor reference or recursively import itself.
import importlib.util
_spec=importlib.util.spec_from_file_location('run_full',PARENT/'raw_sensor_maritime'/'run_full.py')
_reference=importlib.util.module_from_spec(_spec)
sys.modules['run_full']=_reference
_spec.loader.exec_module(_reference)
