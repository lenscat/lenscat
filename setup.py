import setuptools
from pathlib import Path

import re
VERSIONFILE="lenscat/_version.py"
verstrline = open(VERSIONFILE, "rt").read()
VSRE = r"^__version__ = ['\"]([^'\"]*)['\"]"
mo = re.search(VSRE, verstrline, re.M)
if mo:
    verstr = mo.group(1)
else:
    raise RuntimeError("Unable to find version string in %s." % (VERSIONFILE,))


setuptools.setup(
    name="lenscat",
    version=verstr,
    description="A public and community-maintained catalog of known strong gravitational lenses",
    long_description=Path("README.md").read_text(encoding="utf-8"),
    long_description_content_type="text/markdown",
    packages=[
        "lenscat",
    ],
    package_data={
        "lenscat": [
            "data/catalog.csv",
        ]
    },
    install_requires=[
        "astropy >= 6.0.0",
        "ligo.skymap",
    ],
    classifiers=[
        "License :: OSI Approved :: MIT License",
    ],
    python_requires='>=3.9',
)