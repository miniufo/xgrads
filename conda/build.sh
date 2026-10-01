#!/usr/bin/env bash
set -euo pipefail
echo "Building xgrads with pip"
${PYTHON} -m pip install . --no-deps --no-build-isolation -vv
