#!/usr/bin/env bash
set -Eeuo pipefail

readonly EXPECTED_GLIBC="2.28"
readonly PYTHON="python3.9"

if [[ "$(uname -m)" != "x86_64" ]]; then
  echo "error: builder architecture is $(uname -m); expected x86_64" >&2
  exit 1
fi

actual_glibc="$(getconf GNU_LIBC_VERSION | awk '{print $2}')"
if [[ "$actual_glibc" != "$EXPECTED_GLIBC" ]]; then
  echo "error: builder has glibc $actual_glibc; expected exactly $EXPECTED_GLIBC" >&2
  exit 1
fi

echo "Builder: $(cat /etc/redhat-release)"
echo "glibc:   $actual_glibc"
echo "gcc:     $(gcc -dumpfullversion -dumpversion)"
echo "Python:  $($PYTHON --version 2>&1)"
nim --version

# The workspace is mounted from the host. Keep every generated file there so
# upload-artifact can read it after the container exits.
make clean PYTHON="$PYTHON"
make build NIM=nim PYTHON="$PYTHON"

for binary in bin/slimmc bin/slimmc-summary; do
  test -x "$binary"
  file "$binary"
  ldd "$binary"
done

$PYTHON scripts/check_linux_binary.py \
  --glibc-max "$EXPECTED_GLIBC" \
  bin/slimmc bin/slimmc-summary

./bin/slimmc --version
$PYTHON scripts/check_build_provenance.py bin/slimmc --require-git
./bin/slimmc --check homo/tests/regression/frp_all_channels_seeded/FRP_ALLCHAN01.model
./bin/slimmc --check copo/tests/validation/phase_a/models/C01_init.model

$PYTHON scripts/package_binary_release.py \
  --platform linux-x86_64-glibc_2.28

archive="dist/slimmc-$(tr -d '\r\n' < VERSION)-linux-x86_64-glibc_2.28.zip"
test -s "$archive"
$PYTHON -m zipfile --test "$archive"
echo "glibc 2.28 release artifact: $archive"
