"""Extract a subagent's final report from its JSONL output into a file, printing only its size and first line."""
import json, sys
src, dst, marker = sys.argv[1], sys.argv[2], sys.argv[3]
best = ''
def walk(x):
    global best
    if isinstance(x, str):
        if marker in x and len(x) > len(best):
            best = x
    elif isinstance(x, dict):
        for v in x.values():
            walk(v)
    elif isinstance(x, list):
        for v in x:
            walk(v)
for line in open(src, encoding='utf-8', errors='replace'):
    try:
        walk(json.loads(line))
    except Exception:
        pass
i = best.find(marker)
best = best[i:] if i >= 0 else best
open(dst, 'w').write(best.rstrip() + '\n')
print(len(best), best.splitlines()[0] if best else '(none)')
