#!/bin/bash
# pack.sh BUNDLE OUTDIR: pack the relocatable library bundle BUNDLE as OUTDIR/<bundle name>.tar.gz with its
# .sha256 beside it. The archive holds the bundle directory itself, so it unpacks to <bundle name>/.
set -euo pipefail
bundle=$(cd "$1" && pwd); out=$(mkdir -p "$2" && cd "$2" && pwd)
name=$(basename "$bundle")
tar -C "$(dirname "$bundle")" -czf "$out/$name.tar.gz" "$name"
if command -v sha256sum > /dev/null; then sum=(sha256sum); else sum=(shasum -a 256); fi
(cd "$out" && "${sum[@]}" "$name.tar.gz" > "$name.tar.gz.sha256")
echo "packed $out/$name.tar.gz ($(wc -c < "$out/$name.tar.gz" | tr -d ' ') bytes)"
