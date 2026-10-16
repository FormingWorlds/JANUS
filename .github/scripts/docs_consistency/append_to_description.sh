#!/usr/bin/env bash
# Usage: append_to_description.sh <issue-or-pr-number> <section-file>
# Append a dated section to an issue or PR description (edits send no notification).
# Posts a comment instead near GitHub's 65536-character limit.
set -euo pipefail

endpoint="repos/$GITHUB_REPOSITORY/issues/$1"
body=$(mktemp)
gh api "$endpoint" --jq '.body // ""' > "$body"
printf '\n---\n\n_Added by the run on %s:_\n\n' "$(date +%Y-%m-%d)" >> "$body"
cat "$2" >> "$body"

# Bytes are never fewer than characters, so this check is conservative.
if [ "$(wc -c < "$body")" -gt 60000 ]; then
  gh api "$endpoint/comments" -F body=@"$2" > /dev/null
else
  gh api -X PATCH "$endpoint" -F body=@"$body" > /dev/null
fi
