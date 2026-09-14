#!/bin/sh
# Checks that `coral-c/` is the CORAL reference C implementation at the revision
# we measured, unmodified.
#
# `coral-c/` is vendored rather than fetched so that the artifact builds without
# network access; this script is what makes the provenance checkable.  Only
# `coral-c/` is vendored: the repository also carries a 368 MB Sage
# specification, which is not what we measured.  `README.md` and `median.py` are
# ours; `UPSTREAM-README.md` is upstream's top-level README, checked below.
set -e
cd "$(dirname "$0")"
REV=61acaeb832243da630f316de340949d9d091983d   # the revision we measured
D=$(mktemp -d)
git clone -q https://github.com/CORAL-nike/CORAL "$D/CORAL"
git -C "$D/CORAL" checkout -q "$REV"
ok=1
diff -r -x build -x 'build*' coral-c "$D/CORAL/coral-c" || ok=0
diff UPSTREAM-README.md "$D/CORAL/README.md" || ok=0
if [ "$ok" = 1 ]; then
    echo "coral-c/ matches CORAL-nike/CORAL at $REV ($(find coral-c -type f -not -path '*/build/*' | wc -l | tr -d ' ') files)"
else
    echo "coral-c/ DIFFERS from upstream, see above" >&2
fi
rm -rf "$D"
[ "$ok" = 1 ]
