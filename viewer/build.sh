#!/bin/sh
# Build the real-time viewer: compile the engine to WebAssembly and inline it.
#
# The result is one self-contained HTML file with no external requests except
# the web fonts. Open it directly — no server, no toolchain, no install.
set -e
cd "$(dirname "$0")/.."
cargo rustc --release --lib --target wasm32-unknown-unknown --crate-type cdylib
WASM=target/wasm32-unknown-unknown/release/phys.wasm
echo "engine: $(wc -c < "$WASM") bytes"
# Substitution is done in python rather than sed: a 290 kB base64 string is far
# past the shell's argument limit.
python3 - "$WASM" <<'PY'
import base64, pathlib, sys
root = pathlib.Path(__file__).resolve().parent if "__file__" in dir() else pathlib.Path(".")
wasm = pathlib.Path(sys.argv[1]).read_bytes()
b64 = base64.b64encode(wasm).decode()
head = pathlib.Path("viewer/template.head.html").read_text()
tail = pathlib.Path("viewer/template.tail.html").read_text()
out = pathlib.Path("viewer/index.html")
out.write_text(head + tail.replace("__WASM_B64__", b64))
print(f"viewer: {out.stat().st_size} bytes -> {out}")
PY
