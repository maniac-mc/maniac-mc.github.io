#!/bin/bash
set -euo pipefail
# shellcheck source=../../lib.sh
source "$(dirname "${BASH_SOURCE[0]}")/../../lib.sh"
compile_and_run_all_tests