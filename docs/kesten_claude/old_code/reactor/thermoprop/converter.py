def parse_thermo_inp_file(path):
    reactants = {}
    products = {}

    current_section = 'products'
    current_species_name = None
    current_species_lines = []

    with open(path, 'r') as f:
        for line in f:
            stripped = line.strip()

            # Skip completely empty lines
            if not stripped:
                continue

            # Skip the thermo header block
            if stripped.lower().startswith('thermo'):
                continue

            # Detect section boundaries
            if stripped.upper().startswith('END PRODUCTS'):
                if current_species_name and current_species_lines:
                    reactants[current_species_name] = ''.join(current_species_lines)
                current_section = 'reactants'
                current_species_name = None
                current_species_lines = []
                continue

            if stripped.upper().startswith('END REACTANTS'):
                if current_species_name and current_species_lines:
                    products[current_species_name] = ''.join(current_species_lines)
                break  # End of file
                        
            # Detect the start of a species record
            if line.lstrip() and line.lstrip()[0].isalpha():
                # Save the previous species before starting a new one
                if current_species_name and current_species_lines:
                    if current_section == 'reactants':
                        reactants[current_species_name] = ''.join(current_species_lines)
                    elif current_section == 'products':
                        products[current_species_name] = ''.join(current_species_lines)
                
                current_species_name = stripped.split()[0]
                current_species_lines = [line]
                # If section not set yet, assume reactants
                if current_section is None:
                    current_section = 'reactants'
            else:
                current_species_lines.append(line)

    print(f"Reactants found: {len(reactants)}")
    print(f"Products found: {len(products)}")
    return reactants, products


def parse_record_lines(record_str):
    """
    Parse a full NASA-thermo record into structured fields.

    Returns a dict with:
      - name: species name (first 24 chars of first line)
      - comment: comment/data source info (chars 25–80 of first line)
      - lines: list of all lines in the record
    """
    if not record_str:
        return {"name": "", "comment": "", "lines": []}

    lines = record_str.splitlines()
    first_line = lines[0] if lines else ""
    name    = first_line[:24].strip()
    comment = first_line[24:80].rstrip()

    return {
        "name": name,
        "comment": comment,
        "lines": lines
    }

def build_species_dict(reactants, products):
    """
    Combine reactants and products into a structured dictionary.
    """
    species_data = {}
    for section, d in (("reactants", reactants), ("products", products)):
        for _, record in d.items():
            parsed = parse_record_lines(record)
            species_data[parsed["name"]] = {
                "section": section,
                "comment": parsed["comment"],
                "lines": parsed["lines"]
            }
    return species_data

def _safe_float(s):
    s = (s or "").strip().replace("D", "E").replace("d", "E")
    if not s:
        return None
    try:
        return float(s)
    except ValueError:
        return None

def _pad80(s):
    return (s or "").rstrip("\n").ljust(80)

def parse_record_second_line(lines):
    """
    Parse the 2nd line of a NASA-thermo record using 1-based columns:
      col 2                   -> number of T intervals (int)
      cols 4–9                -> optional identification code (str)
      cols 11–50              -> chemical formulas (str)
      col 52                  -> gas/condensed flag (str, 1 char)
      cols 53–65              -> molecular weight (float)
      cols 66–80              -> heat of formation at 298.15 K (float)
    """
    if not lines or len(lines) < 2:
        return {
            "n_T_intervals": None,
            "id_code": "",
            "formula": "",
            "phase_flag": "",
            "molecular_weight": None,
            "hf_298": None,
            "second_line": ""
        }

    line2 = _pad80(lines[1])

    # 1-based -> Python slices
    # col 2
    n_T_intervals_raw = line2[1:2]     # single char
    n_T_intervals = int(n_T_intervals_raw) if n_T_intervals_raw.strip().isdigit() else None

    # cols 4–9
    id_code = line2[3:9].strip()

    # cols 11–50
    formula = line2[10:50].strip()

    # col 52
    phase_flag = line2[51:52].strip()

    # cols 53–65
    molecular_weight = _safe_float(line2[52:65])

    # cols 66–80
    hf_298 = _safe_float(line2[65:80])

    return {
        "n_T_intervals": n_T_intervals,
        "id_code": id_code,
        "formula": formula,
        "phase_flag": phase_flag,
        "molecular_weight": molecular_weight,
        "hf_298": hf_298,
        "second_line": lines[1]
    }


def enrich_with_second_line_fields(species_dict):
    for sp_name, rec in species_dict.items():
        meta = parse_record_second_line(rec.get("lines", []))
        rec.update(meta)
    return species_dict

def _safe_float(s):
    s = (s or "").strip().replace("D", "E").replace("d", "E")
    if not s:
        return None
    try:
        return float(s)
    except ValueError:
        return None

def _split_floats(s):
    s = (s or "").strip().replace("D", "E").replace("d", "E")
    if not s:
        return []
    out = []
    for tok in s.split():
        try:
            out.append(float(tok))
        except ValueError:
            pass
    return out

def _pad80(line):
    return (line or "").rstrip("\n").ljust(80)

def parse_record_line3(lines):
    """
    Line 3 (1-based columns):
      2–21  -> temperature range (two floats; we'll parse all floats we find)
      23    -> number of Cp coefficients (int)
      24–63 -> T exponents for Cp (floats list)
      66–80 -> heat of formation (float)
    """
    if len(lines) < 3:
        return {
            "t_range": [],
            "n_cp_coeff": None,
            "t_exponents": [],
            "hf_line3": None,
            "line3": lines[2] if len(lines) >= 3 else ""
        }

    line3 = _pad80(lines[2])

    # 1-based -> Python slices
    t_range_raw   = line3[1:21]     # cols 2–21
    n_cp_raw      = line3[22:23]    # col 23
    t_exps_raw    = line3[23:63]    # cols 24–63
    hf_line3_raw  = line3[65:80]    # cols 66–80

    t_range = _split_floats(t_range_raw)         # expect [T_low, T_high] typically
    n_cp_coeff = int(n_cp_raw) if n_cp_raw.strip().isdigit() else None
    t_exponents = _split_floats(t_exps_raw)
    hf_line3 = _safe_float(hf_line3_raw)

    return {
        "t_range": t_range,
        "n_cp_coeff": n_cp_coeff,
        "t_exponents": t_exponents,
        "hf_line3": hf_line3,
        "line3": lines[2]
    }

def parse_record_line4(lines):
    """
    Line 4: first five Cp coefficients.
    Fortran format: 5D16.8 over columns 1–80 (i.e., five 16-char fields).
    """
    if len(lines) < 4:
        return {"cp_coeffs_1_5": [], "line4": lines[3] if len(lines) >= 4 else ""}

    line4 = _pad80(lines[3])  # reuse your helper: pad to 80 cols, keep raw too
    fields = [line4[i:i+16] for i in range(0, 80, 16)]
    coeffs = []
    for f in fields:
        v = _safe_float(f)  # handles D/E exponents; returns None if empty/bad
        coeffs.append(v)

    return {"cp_coeffs_1_5": coeffs, "line4": lines[3]}


def parse_record_line5(lines):
    """
    Line 5:
      cols 1–48  -> 3 Cp coefficients (Fortran 3D16.8)  -> cp_coeffs_6_8
      cols 49–80 -> 2 integration constants (2D16.8)    -> b1, b2
    """
    if len(lines) < 5:
        return {
            "cp_coeffs_6_8": [],
            "b1": None,
            "b2": None,
            "line5": lines[4] if len(lines) >= 5 else ""
        }

    line5 = _pad80(lines[4])

    # First 3 * 16-char fields -> cols 1–48
    cp_fields = [line5[i:i+16] for i in range(0, 48, 16)]
    cp_coeffs_6_8 = [_safe_float(f) for f in cp_fields]

    # Next 2 * 16-char fields -> cols 49–80
    b_fields = [line5[i:i+16] for i in range(48, 80, 16)]
    b_vals = [_safe_float(f) for f in b_fields]
    b1 = b_vals[0] if len(b_vals) > 0 else None
    b2 = b_vals[1] if len(b_vals) > 1 else None

    return {
        "cp_coeffs_6_8": cp_coeffs_6_8,
        "b1": b1,
        "b2": b2,
        "line5": lines[4]
    }

def enrich_with_lines3to5(species_dict):
    for sp, rec in species_dict.items():
        lines = rec.get("lines", [])
        rec.update(parse_record_line3(lines))
        rec.update(parse_record_line4(lines))
        rec.update(parse_record_line5(lines))
    return species_dict


def parse_interval_block(lines, start_idx):
    """
    Parse one (line3, line4, line5) block starting at start_idx.
    Returns a dict with fields from line3/4/5 plus combined cp list.
    """
    # Reuse your existing parsers, but pointing at the correct offsets
    # Temporarily slice a view so parse_record_lineX can work unchanged
    block = lines[start_idx:start_idx+3]

    # --- line 3 ---
    # Using the same logic as parse_record_line3 but on 'block'
    if len(block) >= 1:
        l3 = _pad80(block[0])
        t_range      = _split_floats(l3[1:21])        # cols 2–21
        n_cp_coeff   = int(l3[22:23]) if l3[22:23].strip().isdigit() else None
        t_exponents  = _split_floats(l3[23:63])       # cols 24–63
        hf_line3     = _safe_float(l3[65:80])         # cols 66–80
        raw_line3    = block[0]
    else:
        t_range, n_cp_coeff, t_exponents, hf_line3, raw_line3 = [], None, [], None, ""

    # --- line 4 ---
    if len(block) >= 2:
        l4 = _pad80(block[1])
        fields4 = [l4[i:i+16] for i in range(0, 80, 16)]
        cp_coeffs_1_5 = [_safe_float(f) for f in fields4]
        raw_line4 = block[1]
    else:
        cp_coeffs_1_5, raw_line4 = [], ""

    # --- line 5 ---
    if len(block) >= 3:
        l5 = _pad80(block[2])
        cp_fields = [l5[i:i+16] for i in range(0, 48, 16)]   # 3D16.8
        cp_coeffs_6_8 = [_safe_float(f) for f in cp_fields]
        b_fields = [l5[i:i+16] for i in range(48, 80, 16)]   # 2D16.8
        b_vals = [_safe_float(f) for f in b_fields]
        b1 = b_vals[0] if len(b_vals) > 0 else None
        b2 = b_vals[1] if len(b_vals) > 1 else None
        raw_line5 = block[2]
    else:
        cp_coeffs_6_8, b1, b2, raw_line5 = [], None, None, ""

    cp_coeffs_all = (cp_coeffs_1_5 or []) + (cp_coeffs_6_8 or [])

    return {
        "t_range": t_range,
        "n_cp_coeff": n_cp_coeff,
        "t_exponents": t_exponents,
        "hf_line3": hf_line3,
        "line3": raw_line3,
        "cp_coeffs_1_5": cp_coeffs_1_5,
        "line4": raw_line4,
        "cp_coeffs_6_8": cp_coeffs_6_8,
        "b1": b1,
        "b2": b2,
        "cp_coeffs_all": cp_coeffs_all,
        "line5": raw_line5,
    }

def enrich_with_intervals(species_dict):
    """
    For each species record, parse repeated (line3,4,5) blocks based on n_T_intervals.

    - Stores all blocks under rec["intervals"] as a list.
    - Also mirrors the *first* interval's fields back to top-level keys for
      backward compatibility (t_range, cp_coeffs_all, b1/b2, etc.).
    """
    for sp, rec in species_dict.items():
        lines = rec.get("lines", [])
        # prefer explicit value from previously parsed line 2; else infer
        n_int = rec.get("n_T_intervals")
        if n_int is None:
            # infer how many full triplets exist after the first two lines
            n_int = max(0, (len(lines) - 2) // 3)

        intervals = []
        base = 2  # line3 starts at index 2
        for i in range(n_int):
            idx = base + 3*i
            intervals.append(parse_interval_block(lines, idx))

        rec["intervals"] = intervals

        # mirror the first interval to top-level for convenience
        if intervals:
            first = intervals[0]
            rec.update({
                "t_range": first["t_range"],
                "n_cp_coeff": first["n_cp_coeff"],
                "t_exponents": first["t_exponents"],
                "hf_line3": first["hf_line3"],
                "cp_coeffs_1_5": first["cp_coeffs_1_5"],
                "cp_coeffs_6_8": first["cp_coeffs_6_8"],
                "cp_coeffs_all": first["cp_coeffs_all"],
                "b1": first["b1"],
                "b2": first["b2"],
                "line3": first["line3"],
                "line4": first["line4"],
                "line5": first["line5"],
            })
        else:
            rec.setdefault("intervals", [])
    return species_dict

import re

def _phase_flag_to_int(ch):
    # You asked for an int flag; treat 'G' as gas=1, everything else condensed=0
    return 1 if (ch or "").strip().upper() == "G" else 0

_formula_token = re.compile(r"([A-Z][a-z]?)([-+]?\d*\.?\d*(?:[Ee][+-]?\d+)?)")

def parse_formula_to_dict(formula_str):
    """
    Convert a formula string like 'AL' or 'C 1 H 4' or 'Al2O3' into {'AL':1.0}, {'C':1.0,'H':4.0}, {'Al':2,'O':3}
    Loose but effective: ignores dots, separators, and parentheses content (no isotope handling).
    """
    if not formula_str:
        return {}
    s = formula_str.replace("·", " ").replace(".", " ").replace(",", " ").replace("(", " ").replace(")", " ")
    parts = {}
    for sym, num in _formula_token.findall(s):
        count = float(num) if num not in ("", "+", "-") else 1.0
        parts[sym] = parts.get(sym, 0.0) + count
    return parts

def _range_key(tlow, thigh):
    # match your example style like '200-1000'
    def fmt(x):
        return str(int(x)) if abs(x - int(x)) < 1e-9 else f"{x:g}"
    return f"{fmt(tlow)}-{fmt(thigh)}"

def to_template_species_dict(species_dict):
    """
    Transform your parsed species_dict into the requested template structure.
    """
    out = {}

    for _, rec in species_dict.items():
        name = rec.get("name", "")
        comment = rec.get("comment", "")

        # line 2 fields
        n_intervals = rec.get("n_T_intervals")
        ident_code  = rec.get("id_code", "")
        formula_raw = rec.get("formula", "")
        mw          = rec.get("molecular_weight")
        hf298       = rec.get("hf_298")
        phase_ch    = rec.get("phase_flag", "")

        # build ranges
        ranges = {}
        for it in rec.get("intervals", []):
            tlo, thi = (it.get("t_range") + [None, None])[:2]
            if tlo is None or thi is None:
                continue
            key = _range_key(tlo, thi)
            coeff = list(it.get("cp_coeffs_all") or [])
            # append integration constants if present
            if it.get("b1") is not None:
                coeff.append(it["b1"])
            if it.get("b2") is not None:
                coeff.append(it["b2"])

            ranges[key] = {
                "t_low": tlo,
                "t_high": thi,
                "num_coeff": it.get("n_cp_coeff"),
                "hf": it.get("hf_line3"),
                "coeff": coeff,
                "T_exp": it.get("t_exponents") or []
            }

        # final object for this species
        sp_key = (name.split() or [name])[0]  # e.g., 'AL' from 'AL                Cons:'
        out[sp_key] = {
            "name": name,
            "comment": comment,
            "num_T_intervals": n_intervals,
            "ident_code": ident_code,
            "formula": parse_formula_to_dict(formula_raw),
            "phase_flag": _phase_flag_to_int(phase_ch),
            "mw": mw,
            "HF": hf298,
            "ranges": ranges,
        }

    return out



if __name__ == '__main__':
    path = 'raw/thermo.inp'
    
    reactants, products = parse_thermo_inp_file(path)
    
    species = build_species_dict(reactants, products)
    species = enrich_with_second_line_fields(species)
    species = enrich_with_intervals(species)
    
    data = to_template_species_dict(species)
    
    # sanity check one species
    k = next(iter(data))
    print(k, data[k]["name"])
    print("interval keys:", list(data[k]["ranges"].keys()))






    data = {'AL':   {'name': 'AL                Cons:',
                     'comment': 'JPCRD v20 n5 p775 1991. Hf:CODATA,1989,p24. T-fit',
                     'num_T_intervals': 3,
                     'ident_code': 'l 5/97',
                     'formula': {'AL': 1.0},
                     'phase_flag': 0,
                     'mw': 26.98154,
                     'HF': 330000.000,
                     'ranges': {'200-1000'  : {'t_low': 200.000,  't_high': 1000.000, 'num_coeff': 7, 'hf': 6918.671,  
                                               'coeff': [5.006608890e+03, 1.861304407e+01, 2.412531111e+00, 1.987604647e-04, -2.432362152e-07, 1.538281506e-10, -3.944375734e-14, 0.000000000e+00],
                                               'T_exp': [-2.0, -1.0, 0.0, 1.0, 2.0, 3.0, 4.0, 0.0]},
                                '1000-6000' : {'t_low': 1000.000, 't_high': 6000.000, 'num_coeff': 7, 'hf': 6918.671,
                                               'coeff': [828097610e+04, 1.140929691e+02, 2.359891025e+00, 7.574401020e-05, -1.483585474e-08, -1.060430572e-12, 5.086638598e-16, 0.000000000e+00 , 3.825003380e+04, 6.579504480e+00],
                                               'T_exp': [-2.0, -1.0, 0.0, 1.0, 2.0, 3.0, 4.0, 0.0]},
                                '6000-20000': {'t_low': 6000.000, 't_high': 20000.000, 'num_coeff': 7, 'hf': 6918.671,
                                               'coeff': [-5.038989040e+08, 3.801217920e+05, -1.082073372e+02, 1.549111860e-02, -1.069893863e-06, 3.591451830e-11, -4.695230800e-16, 0.000000000e+00, -2.900178789e+06, 9.489506830e+02],
                                               'T_exp': [-2.0, -1.0, 0.0, 1.0, 2.0, 3.0, 4.0, 0.0]}}
                     
                     
                    }
            
           }    
    
    
    
    
    
    
    
    
    
    
    