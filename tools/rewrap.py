"""Rewrap a markdown draft's prose paragraphs to the repository's width (140), leaving indented display lines, tables and headings."""
import sys, textwrap
p = sys.argv[1]
W = int(sys.argv[2]) if len(sys.argv) > 2 else 140
out, para = [], []

def flush():
    if para:
        out.extend(textwrap.wrap(' '.join(s.strip() for s in para), width=W, break_long_words=False, break_on_hyphens=False))
        para.clear()

for line in open(p).read().split('\n'):
    if not line.strip() or line.startswith('    ') or line.startswith('|') or line.startswith('#'):
        flush(); out.append(line)
    else:
        para.append(line)
flush()
open(p, 'w').write('\n'.join(out))
