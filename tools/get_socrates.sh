#!/bin/bash
# Download and compile SOCRATES radiative transfer code.
# Optionally pass a target directory as the first argument (default: ./socrates).
#
# Usage:
#   bash get_socrates.sh
#   bash get_socrates.sh /path/to/install/socrates

set -e

# ── Target directory ────────────────────────────────────────────────────────
socpath="${1:-socrates}"
socpath="$(realpath "$socpath")"

# ── SSH check ────────────────────────────────────────────────────────────────
# ssh -T git@github.com exits with code 1 when auth succeeds (no shell granted)
# and with code 255 when auth fails or no SSH key is present.
use_ssh=false
if ssh -T git@github.com 2>/dev/null; then
    : # exit 0 would be unusual, but treat as success
elif [ $? -eq 1 ]; then
    use_ssh=true
fi

if [ "$use_ssh" = true ]; then
    uri="git@github.com:FormingWorlds/SOCRATES.git"
else
    uri="https://github.com/FormingWorlds/SOCRATES.git"
fi

# ── Clone ────────────────────────────────────────────────────────────────────
if [ -d "$socpath" ]; then
    echo "Directory $socpath already exists — removing it"
    rm -rf "$socpath"
fi

echo "Cloning SOCRATES"
echo "  $uri -> $socpath"
git clone "$uri" "$socpath"

# ── Build ────────────────────────────────────────────────────────────────────
cd "$socpath"

logfile="$socpath/build.log"
echo ""
echo "Build log: $logfile"

# Configure — output to log only
echo ""
echo "Configuring SOCRATES..."
./configure >> "$logfile" 2>&1

# Build — full output goes to log, only actual errors printed to terminal.
# Disable set -e around the pipe so grep non-matches don't trigger exit.
echo "Building SOCRATES (this may take a few minutes)..."
set +e
./build_code 2>&1 | tee -a "$logfile" | grep -iE "error:|fatal error:|undefined reference|cannot find|ld returned" >&2
build_exit=${PIPESTATUS[0]}
set -e

if [ "$build_exit" -ne 0 ]; then
    echo ""
    echo "Build failed (exit code $build_exit). See full log:"
    echo "  $logfile"
    exit "$build_exit"
fi

# ── Done ─────────────────────────────────────────────────────────────────────
echo ""
echo "SOCRATES has been built in: $socpath"
echo "Full build log: $logfile"
echo ""
echo "Add the following to your ~/.bashrc (or ~/.zshrc):"
echo ""
echo "    export RAD_DIR=$socpath"
echo ""