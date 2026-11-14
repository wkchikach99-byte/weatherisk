"""Write the upload copy of the v24 manuscript: revision macros inlined, comments removed."""
import re
import sys

src, dst = sys.argv[1], sys.argv[2]
s = open(src).read()

s = re.sub(r"\n% --- THORSTEN-18 TRACKING VERSION ---.*?\\newcommand\{\\revv\}\[1\]\{\{#1\}\}\n",
           "\n", s, flags=re.S)
assert "\\newcommand{\\rev" not in s and "RevisionRed" not in s

macro = re.compile(r"\\rev[a-z]?\{")
while True:
    m = macro.search(s)
    if not m:
        break
    i, depth = m.end(), 1
    while depth:
        c = s[i]
        if c == "\\":
            i += 2
            continue
        depth += {"{": 1, "}": -1}.get(c, 0)
        i += 1
    s = s[:m.start()] + s[m.end():i - 1] + s[i:]

lines = []
for line in s.split("\n"):
    if line.lstrip().startswith("%"):
        continue
    m = re.search(r"(?<!\\)%", line)
    if m and not line.rstrip().endswith("%"):
        line = line[:m.start()].rstrip()
    lines.append(line)
s = re.sub(r"\n{3,}", "\n\n", "\n".join(lines))
open(dst, "w").write(s)
