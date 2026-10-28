#!/usr/bin/env bash
# Usage: append_to_description.sh <issue-or-pr-number> <section-file>
# Append a dated section to an issue or PR description (edits send no notification).
# Posts it as a comment instead near GitHub's 65536-character limit.
set -euo pipefail

endpoint="repos/$GITHUB_REPOSITORY/issues/$1"
section=$(mktemp)
printf '_Added by the run on %s:_\n\n' "$(date +%Y-%m-%d)" > "$section"
cat "$2" >> "$section"

# Byte counts are never below character counts, so these checks are conservative.
if [ "$(wc -c < "$section")" -gt 60000 ]; then
  # iconv -c drops a character cut in half; it exits 1 when it does.
  head -c 60000 "$section" | iconv -c -f UTF-8 -t UTF-8 > "$section.cut" || true
  printf '\n\n_Cut to fit the size limit; the full findings are in this run'"'"'s analysis artifact._\n' \
    >> "$section.cut"
  mv "$section.cut" "$section"
fi

body=$(mktemp)
gh api "$endpoint" --jq '.body // ""' > "$body"
printf '\n---\n\n' >> "$body"
cat "$section" >> "$body"

if [ "$(wc -c < "$body")" -gt 60000 ]; then
  gh api "$endpoint/comments" -F body=@"$section" > /dev/null
else
  gh api -X PATCH "$endpoint" -F body=@"$body" > /dev/null
fi
