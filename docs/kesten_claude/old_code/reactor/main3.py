import re
from pathlib import Path
import tomllib


def _to_float_list(s):
    nums = []
    for tok in s.strip().split():
        tok = tok.replace('D','E').replace('d','e')
        nums.append(float(tok))
    return nums

def _parse_interval(lineA, lineB, lineC):
    # Header: T_low T_high 7  exponents(8)  [maybe trailing number]
    header_vals = _to_float_list(lineA)
    if len(header_vals) < 2:
        raise ValueError("Interval header too short: %r" % lineA)
    T_low, T_high = header_vals[0], header_vals[1]

    # Try to find the "7" marker and read next 8 numbers as exponents
    exponents = None
    if 7.0 in header_vals:
        i7 = header_vals.index(7.0)
        exponents = header_vals[i7+1:i7+9]
    if not exponents or len(exponents) < 8:
        # Fallback: take last 8 numbers in header as exponents
        exponents = header_vals[-8:]

    # Coefficients: NASA-7 (a1..a7): 5 on next line, then 2 at start of the following line
    valsB = _to_float_list(lineB)
    valsC = _to_float_list(lineC)
    a = valsB[:5] + valsC[:2]
    extras = valsC[2:] if len(valsC) > 2 else []

    # Coerce exponent list to ints when they look integral (e.g., -2.0 -> -2)
    exp_clean = []
    for x in exponents:
        xi = int(round(x))
        exp_clean.append(xi if abs(x - xi) < 1e-12 else x)

    return {
        'T_low': T_low,
        'T_high': T_high,
        'exponents': exp_clean,
        'a': a,
        'extras': extras,
    }

def _parse_elements_and_tail(line2):
    # Grab element/count pairs flexibly (H 2.00, O 1.00, E -1.00, etc.)
    elem_pairs = re.findall(r'([A-Za-z][a-z]?)\s*([+-]?\d+(?:\.\d*)?)', line2)
    elements = {}
    for sym, val in elem_pairs:
        try:
            elements[sym] = float(val)
        except ValueError:
            pass

    # MW and Hf298: take the last two floats on the line (robust across variants)
    floats_l2 = _to_float_list(re.sub(r'[A-Za-z/]', ' ', line2))
    mw, hf298 = (floats_l2[-2], floats_l2[-1]) if len(floats_l2) >= 2 else (None, None)
    return elements, mw, hf298

def parse_cea_block(block_lines):
    """
    Parse one species block (header + line2 + N interval triplets).
    Returns dict with 'name','comment','elements','mw','hf298','ranges'.
    """
    if len(block_lines) < 5:
        raise ValueError("Block too short")

    # Line 1: name + the rest is comment
    l1 = block_lines[0].rstrip('\n')
    m = re.match(r'\s*([^\s]+)\s*(.*)$', l1)
    name = m.group(1)
    comment = m.group(2).rstrip()

    # Line 2: elements, MW, Hf298 (layout varies slightly)
    l2 = block_lines[1].rstrip('\n')
    elements, mw, hf298 = _parse_elements_and_tail(l2)

    # Remaining lines should be groups of 3 per interval
    ranges = []
    i = 2
    while i + 2 < len(block_lines):
        A = block_lines[i].rstrip('\n')
        B = block_lines[i+1].rstrip('\n')
        C = block_lines[i+2].rstrip('\n')
        ranges.append(_parse_interval(A, B, C))
        i += 3

    return {
        'name': name,
        'comment': comment,
        'elements': elements,
        'mw': mw,
        'hf298': hf298,
        'ranges': ranges
    }

def _is_new_record_header(line):
    # Heuristic: species header starts in column 0 with a non-space char (name),
    # whereas all interval lines are indented. Works for standard thermo.inp.
    return bool(line) and not line.startswith(' ')

def parse_cea_text(text):
    """
    Parse an entire thermo text (multiple species).
    Returns dict keyed by species name; each value is the species dict.
    """
    lines = [ln.rstrip('\n') for ln in text.splitlines() if ln.strip() != '']
    blocks = []
    cur = []
    for ln in lines:
        if _is_new_record_header(ln):
            if cur:
                blocks.append(cur)
                cur = []
        cur.append(ln)
    if cur:
        blocks.append(cur)

    out = {}
    for blk in blocks:
        sp = parse_cea_block(blk)
        out[sp['name']] = sp
    return out

def parse_cea_file(path):
    text = Path(path).read_text(encoding='utf-8', errors='ignore')
    return parse_cea_text(text)

# ---------- quick demo ----------
if __name__ == "__main__":
    sample = """H2O               CODATA,1989. JRNBS v92,1987,p35. TRC tuv-25,10/88.            
 2 l 8/89 H   2.00O   1.00    0.00    0.00    0.00 0     18.01528    -241826.000
    200.000  1000.000 7 -2.0 -1.0  0.0  1.0  2.0  3.0  4.0  0.0         9904.092
 -3.94795999D+04  5.75572977D+02  9.31783351D-01  7.22271091D-03 -7.34255448D-06
  4.95504134D-09 -1.33693261D-12  0.00000000D+00 -3.30397425D+04  1.72420539D+01
   1000.000  6000.000 7 -2.0 -1.0  0.0  1.0  2.0  3.0  4.0  0.0         9904.092
  1.03497224D+06 -2.41269895D+03  4.64611114D+00  2.29199814D-03 -6.83683007D-07
  9.42646842D-11 -4.82238028D-15  0.00000000D+00 -1.38428625D+04 -7.97815119D+00
H2O+              Cons & Hf298: TPIS,v1,pt1,1989,p125.                          
 3 tpis89 H   2.00O   1.00E  -1.00    0.00    0.00 0     18.01473     981600.000
    200.000  1000.000 7 -2.0 -1.0  0.0  1.0  2.0  3.0  4.0  0.0         9934.466
 -4.21440824D+04  6.75359651D+02 -1.96925514D-02  1.07169694D-02 -1.28088946D-05
  9.21018233D-09 -2.65338672D-12  0.00000000D+00  1.13695145D+05  2.29367287D+01
   1000.000  6000.000 7 -2.0 -1.0  0.0  1.0  2.0  3.0  4.0  0.0         9934.466
  6.26589105D+05 -2.87353484D+03  7.72634464D+00 -9.06918399D-04  6.18773526D-07
 -1.20274764D-10  7.41410714D-15  0.00000000D+00  1.34269245D+05 -2.64305226D+01
   6000.000 20000.000 7 -2.0 -1.0  0.0  1.0  2.0  3.0  4.0  0.0         9934.466
 -1.75326556D+07 -2.63657336D+03  1.16561964D+01 -7.73121836D-04  5.63107918D-08
 -1.96513165D-12  2.68205807D-17  0.00000000D+00  1.16847509D+05 -5.76343070D+01
H2O2              Cons: TPIS,v1,pt1,p121,1978. Hf298: TRC w-31,6/88.            
 2 l 2/93 H   2.00O   2.00    0.00    0.00    0.00 0     34.01468    -136310.000
    200.000  1000.000 7 -2.0 -1.0  0.0  1.0  2.0  3.0  4.0  0.0        11158.835
 -9.27966807D+04  1.56476669D+03 -5.97655534D+00  3.27076867D-02 -3.93222494D-05
  2.50927583D-08 -6.46509773D-12  0.00000000D+00 -2.49918528D+04  5.87722783D+01
   1000.000  6000.000 7 -2.0 -1.0  0.0  1.0  2.0  3.0  4.0  0.0        11158.835
  1.48940442D+06 -5.17075504D+03  1.12819764D+01 -8.03840050D-05 -1.81951973D-08
  6.94886179D-12 -4.82870187D-16  0.00000000D+00  1.41303656D+04 -4.65080305D+01
"""
    db = parse_cea_text(sample)
    # Minimal sanity check prints
    for k, v in db.items():
        print(k, "MW=", v['mw'], "Hf298=", v['hf298'], "#ranges=", len(v['ranges']))
        for r in v['ranges']:
            print("  ", r['T_low'], "->", r['T_high'], "a1..a7 len:", len(r['a']))
