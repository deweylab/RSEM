#!/usr/bin/env bash
# Assert RSEM honored a --<aligner>-path flag instead of falling back to PATH.
# Gold diffs cannot detect this: pointing the flag at a directory resolving the
# same binary PATH would have found produces identical output.
#
# Usage: assert_path_flag_honored.sh <probe-log> <probe-dir> <tool> <flag>
set -euo pipefail

log="${1:?probe log path required}"
probe_dir="${2:?probe dir required}"
tool="${3:?tool name required}"
flag="${4:?flag name required}"

if [[ ! -s "$log" ]]; then
  echo "FAIL: $probe_dir/$tool was never invoked -- $flag was ignored." >&2
  exit 1
fi

# Without this the check above proves nothing: the probe could have been found
# by an ordinary PATH lookup, with the flag contributing nothing.
resolved="$(command -v "$tool" 2>/dev/null || true)"
if [[ -n "$resolved" && "$resolved" == "$probe_dir/"* ]]; then
  echo "FAIL: $probe_dir is on PATH ($resolved); keep it off PATH." >&2
  exit 1
fi

echo "  probe: $flag honored -- $(wc -l <"$log" | tr -d ' ') invocation(s) of $probe_dir/$tool"
