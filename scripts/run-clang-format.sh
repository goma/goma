#!/usr/bin/env bash
CLANG_FORMAT="${CLANG_FORMAT:=clang-format}"
cd "$(dirname "$BASH_SOURCE")/.."
find src/ include/ tests/ -iname "*.h" -print0 | xargs -t -n 1 -P 8 -0 $CLANG_FORMAT -i
find src/ include/ tests/ -iname "*.cpp" -print0 | xargs -t -n 1 -P 8 -0 $CLANG_FORMAT -i
find src/ include/ tests/ -iname "*.c" -print0 | xargs -t -n 1 -P 8 -0 $CLANG_FORMAT -i
