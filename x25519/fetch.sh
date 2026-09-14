#!/bin/sh
# Checks that `ref10/` is the upstream reference implementation of X25519,
# unmodified.
#
# `ref10/` is D. J. Bernstein's `crypto_scalarmult/curve25519/ref10` from
# SUPERCOP, in the public domain (see `ref10/README`).  It is vendored rather
# than fetched because it is 1 300 lines and SUPERCOP is distributed as dated
# tarballs rather than as a git repository, so there is no revision to pin.  This
# script fetches the same files from a mirror of SUPERCOP and diffs them, which is
# the check that vendoring would otherwise skip.
#
# The five headers in `glue/` are ours, not upstream, and are therefore not
# checked here; each says in a comment what it does and why it is needed.
set -e
cd "$(dirname "$0")"
U=https://raw.githubusercontent.com/jedisct1/supercop/master/crypto_scalarmult/curve25519/ref10
D=$(mktemp -d)
for f in $(ls ref10); do
    curl -sSf "$U/$f" -o "$D/$f"
done
if diff -r ref10 "$D"; then
    echo "ref10/ matches upstream ($(ls ref10 | wc -l | tr -d ' ') files)"
else
    echo "ref10/ DIFFERS from upstream, see above" >&2
    exit 1
fi
rm -rf "$D"
