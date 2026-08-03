#!/usr/bin/env bash
# Download STAR 2.7.6a (Bioconda build) into tests/tools/star-2.7.6a/bin/STAR
set -euo pipefail
ROOT="$(cd "$(dirname "$0")/.." && pwd)"
DIR="$ROOT/tests/tools/star-2.7.6a"
PKG="$DIR/star-2.7.6a-0.tar.bz2"
mkdir -p "$DIR"

if [[ -x "$DIR/bin/STAR" ]]; then
  echo "STAR 2.7.6a already present at $DIR/bin/STAR"
  exit 0
fi

os="$(uname -s)"
arch="$(uname -m)"
url=""
case "$os-$arch" in
  Darwin-arm64|Darwin-x86_64)
    url="https://conda.anaconda.org/bioconda/osx-64/star-2.7.6a-0.tar.bz2"
    ;;
  Linux-x86_64)
    url="https://conda.anaconda.org/bioconda/linux-64/star-2.7.6a-0.tar.bz2"
    ;;
  *)
    echo "No bundled STAR 2.7.6a conda package for $os $arch." >&2
    echo "Install STAR 2.7.6a yourself and set STAR_276A_DIR to a directory containing a 'STAR' executable (see Makefile)." >&2
    exit 1
    ;;
esac

echo "Downloading STAR 2.7.6a from $url"
curl -fsSL -o "$PKG" "$url"
tar -xjf "$PKG" -C "$DIR" bin/STAR
rm -f "$PKG"
chmod +x "$DIR/bin/STAR"
echo "Installed $DIR/bin/STAR"
"$DIR/bin/STAR" --version
