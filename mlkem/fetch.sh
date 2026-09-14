#!/bin/sh
# Checks that `ref/` is the pq-crystals reference C implementation of ML-KEM-768
# at the revision we measured, unmodified.
#
# The `standard` branch is the FIPS 203 version.  The default branch of that
# repository is still round-3 Kyber, which is a different scheme, so the branch
# matters.
#
# `ref/` is vendored rather than fetched so that the artifact builds without
# network access; this script is what makes the provenance checkable.  `bench.c`
# and `README.md` are ours and are not checked here.
set -e
cd "$(dirname "$0")"
REV=d5b791c0c601b543233daccbae2845c6197a9e77   # `standard`, the revision we measured
D=$(mktemp -d)
git clone -q -b standard https://github.com/pq-crystals/kyber "$D/kyber"
git -C "$D/kyber" checkout -q "$REV"
ok=1
diff -r ref "$D/kyber/ref" || ok=0
for f in LICENSE AUTHORS SHA256SUMS; do diff "$f" "$D/kyber/$f" || ok=0; done
if [ "$ok" = 1 ]; then
    echo "ref/ matches pq-crystals/kyber at $REV ($(find ref -type f | wc -l | tr -d ' ') files)"
else
    echo "ref/ DIFFERS from upstream, see above" >&2
fi
rm -rf "$D"
[ "$ok" = 1 ]
