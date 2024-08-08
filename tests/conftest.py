import os

if not os.environ.get('FWL_DATA'):
    raise Exception(
        'The FWL_DATA environment variable where spectral and evolution tracks data will be downloaded needs to be set up!'
    )
