#!/usr/bin/env bash
# Validate line limits on the repository's agent-instruction files.
#
# The cap on .github/copilot-instructions.md exists so the file stays
# readable as an entry point; the detailed rule files live under
# .github/.claude/rules/ and carry a separate, more generous cap.

set -e

AGENTS_MAX=750
RULES_MAX=1500

EXIT_CODE=0

if [ -f ".github/copilot-instructions.md" ]; then
    AGENTS_LINES=$(wc -l < .github/copilot-instructions.md | tr -d ' ')
    if [ "$AGENTS_LINES" -gt "$AGENTS_MAX" ]; then
        echo "ERROR: .github/copilot-instructions.md exceeds $AGENTS_MAX lines (current: $AGENTS_LINES)"
        EXIT_CODE=1
    else
        echo "OK: .github/copilot-instructions.md has $AGENTS_LINES lines (max: $AGENTS_MAX)"
    fi
else
    echo "WARNING: .github/copilot-instructions.md not found"
fi

if [ -d ".github/.claude/rules" ]; then
    for rule_file in .github/.claude/rules/*.md; do
        [ -f "$rule_file" ] || continue
        RULE_LINES=$(wc -l < "$rule_file" | tr -d ' ')
        if [ "$RULE_LINES" -gt "$RULES_MAX" ]; then
            echo "ERROR: $rule_file exceeds $RULES_MAX lines (current: $RULE_LINES)"
            EXIT_CODE=1
        else
            echo "OK: $rule_file has $RULE_LINES lines (max: $RULES_MAX)"
        fi
    done
fi

exit $EXIT_CODE
