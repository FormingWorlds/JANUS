import os
from pathlib import Path

if not os.environ.get("RAD_DIR"):
    raise Exception(
        "Socrates environment variables not set! Have you installed Socrates and sourced set_rad_env?"
    )
if not os.environ.get("FWL_DATA"):
    raise Exception(
        "The FWL_DATA environment variable where spectral and evolution tracks data will be downloaded needs to be set up!"
    )
