#!/bin/bash
# Compile each header of the H2 code and the GPU backends alone (-fsyntax-only)
# with the flags of the double library: does it include what it uses?
# Run it after header changes, in both modes (cpu needs build/, gpu build_gpu/,
# each already configured so its SRC_DOUBLE flags.make exists).
# Usage: selfcheck_headers.sh gpu|cpu [header ...]
mode=$1; shift
REPO=$(cd "$(dirname "${BASH_SOURCE[0]}")/.." && pwd)
export S=${TMPDIR:-/tmp}/$USER/bpack_selfcheck_$mode
mkdir -p $S
rm -f $S/*.result
if [ $mode = gpu ]; then B=$REPO/build_gpu; else B=$REPO/build; fi
F=$B/SRC_DOUBLE/CMakeFiles/dbutterflypack.dir/flags.make
export CHECK_FLAGS="$(sed -n 's/^CXX_DEFINES = //p' $F) $(sed -n 's/^CXX_INCLUDES = //p' $F) $(sed -n 's/^CXX_FLAGS = //p' $F | sed 's/-O3//')"
export REPO
cd $REPO
if [ $# -eq 0 ]; then
  set -- $(git ls-files 'h2_parallel/*.hpp' 'GPU_BACKEND/*.hpp' | grep -v '/tests/')
fi
one() {
  h=$1
  name=$(echo $h | tr '/' '_')
  printf '#include "%s"\n' "$REPO/$h" > $S/$name.cpp
  if CC $CHECK_FLAGS -fsyntax-only -Wno-unused -x c++ $S/$name.cpp > $S/$name.log 2>&1; then
    echo "ok   $h" > $S/$name.result
  else
    echo "FAIL $h  ($(grep -c 'error' $S/$name.log) errors; first: $(grep -m1 'error' $S/$name.log | sed 's#.*/##' | cut -c1-160))" > $S/$name.result
  fi
}
export -f one
printf '%s\n' "$@" | xargs -P 16 -I{} bash -c 'one {}'
cat $S/*.result | sort -k2
echo "failures: $(grep -c '^FAIL' $S/*.result | awk -F: '{s+=$2} END{print s}')"
