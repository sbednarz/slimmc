#!/usr/bin/env bash
set -euo pipefail
mkdir -p bin

FLAGS="-d:release --mm:arc --passC:-march=native"

echo ">>> slimmc"
nim c $FLAGS -o:bin/slimmc ./slimmc/slimmc.nim

echo ">>> slimmc-turbo"
nim c $FLAGS -d:danger -d:prg=slimmc-turbo \
  -d:extra="WARNING: turbo build, runtime checks off" \
  -o:bin/slimmc-turbo ./slimmc/slimmc.nim
