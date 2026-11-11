"""The numbers a draft's added lines quote that its results md does not print — a check that a draft quotes the md verbatim (the md writes
a·10ⁿ as a e±nn). Numbers also in the removed lines are skipped; citations (arXiv ids, DOIs, years) show up and are judged by eye.

    venv/bin/python tools/numcheck.py review/drafts/D58.md results/chaos/width_mismatch.md [git revision, default HEAD]"""
import collections
import re
import subprocess
import sys

SUP = str.maketrans('⁰¹²³⁴⁵⁶⁷⁸⁹⁻⁺', '0123456789-+')


def norm(s):
    s = re.sub(r'(\d+(?:\.\d+)?)·10([⁻⁺]?[⁰¹²³⁴⁵⁶⁷⁸⁹]+)', lambda m: m.group(1) + 'e' + m.group(2).translate(SUP), s).replace('−', '-')
    out = []
    for m in re.finditer(r'(?<![\w.])(\d+(?:\.\d+)?)(?:e([-+]?)(\d+))?', s):
        b, sg, ex = m.groups()
        out.append(b if ex is None else f"{b}e{'-' if sg == '-' else '+'}{int(ex):02d}")
    return out


draft, md = sys.argv[1], sys.argv[2]
rev = sys.argv[3] if len(sys.argv) > 3 else 'HEAD'
diff = subprocess.run(['git', 'diff', '-U0', rev, '--', draft], capture_output=True, text=True, check=True).stdout.split('\n')
added = ' '.join(l[1:] for l in diff if l.startswith('+') and not l.startswith('+++'))
removed = collections.Counter(norm(' '.join(l[1:] for l in diff if l.startswith('-') and not l.startswith('---'))))
known = set(norm(open(md, encoding='utf-8').read()))
miss = [t for t in dict.fromkeys(norm(added)) if t not in removed and t not in known]
print('not in the md:', miss if miss else 'none')
