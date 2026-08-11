import re

class ThermoDBConverter():
   
    def __init__(self, path, autoload=True, include_b_in_coeff=True, skip_empty=False):
        self.path               = path              # [str] Location of thermo.inp file
        self.reactants          = {}                # [dict] Reactants Species Records
        self.products           = {}                # [dict] Products Species Records
        self.lines              = []                # [list] List of string lines
        self.include_b_in_coeff = include_b_in_coeff
        self.skip_empty         = skip_empty
        
        # Compile regex for parsing chemical formulas
        self._formula_token     = re.compile(r"([A-Z][a-z]?)([-+]?\d*\.?\d*(?:[Ee][+-]?\d+)?)")
        
        if autoload:
            self.load_datafile()                    # Load the thermo.inp file into a list of line strings
            self.parse_file()                       # Generate a list of of reactant records and product recoards
    
    # ----- Public Loaders ----------------------------------------------------------------------------------
    def load_datafile(self):
        
        with open(self.path, 'r') as f:
            self.lines = f.readlines()

    # ----- Parsing Methods ---------------------------------------------------------------------------------
    def parse_file(self):
        current_section       = 'products'   # start in products; first END flips to reactants
        current_species_name  = None
        current_species_lines = []
    
        # keep your header trimming
        self._remove_header_rows(num_rows=2)
    
        for line in self.lines:
            stripped = line.strip()
            if not stripped:
                continue
    
            # ----- Section boundaries -----
            if stripped.upper().startswith('END PRODUCTS'):
                # flush the last *product* before switching
                if current_species_name and current_species_lines:
                    self.products[current_species_name] = ''.join(current_species_lines)
                current_section = 'reactants'
                current_species_name  = None
                current_species_lines = []
                continue
    
            if stripped.upper().startswith('END REACTANTS'):
                # flush the last *reactant* and stop
                if current_species_name and current_species_lines:
                    self.reactants[current_species_name] = ''.join(current_species_lines)
                break
    
            # ----- Species record start? -----
            if line.lstrip() and line.lstrip()[0].isalpha():
                # save the previous species (to whichever section we're in)
                if current_species_name and current_species_lines:
                    target = self.reactants if current_section == 'reactants' else self.products
                    target[current_species_name] = ''.join(current_species_lines)
    
                # start a new species block
                current_species_name  = stripped.split()[0]
                current_species_lines = [line]
            else:
                # continuation line
                current_species_lines.append(line)
    
        # Flush any trailing species if file didn’t end with END REACTANTS
        # I think this can be deleted since there will not be any species after the 'END REACTANTS' line
        if current_species_name and current_species_lines:
            target = self.reactants if current_section == 'reactants' else self.products
            if current_species_name not in target:
                target[current_species_name] = ''.join(current_species_lines)
    
    # ----- Parsing Helper Methods --------------------------------------------------------------------------
    def _remove_header_rows(self, num_rows):
        if num_rows > 0 and len(self.lines) >= num_rows:            # Make sure the file isn't empty
            self.lines = self.lines[num_rows:]                      # Remove the first 'num_rows' of the list of lines

    
    # ----- Record Parsing Helper Methods -------------------------------------------------------------------
    def _pad80(self, s): 
        return (s or "").rstrip("\n").ljust(80)

    def _safe_float(self, s):
        s = (s or "").strip().replace("D", "E").replace("d", "E")
        if not s:
            return None
        try:
            return float(s)
        except ValueError:
            return None

    def _split_floats_field(self, s):
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
    
    def _range_key(self, tlow, thigh):
        def fmt(x):
            return str(int(x)) if x is not None and abs(x - int(x)) < 1e-9 else f"{x:g}"
        return f"{fmt(tlow)}-{fmt(thigh)}"
    
    def _parse_formula_dict(self, formula_str):
        if not formula_str:
            return {}
        s = (formula_str.replace("·", " ")
                          .replace(".", " ")
                          .replace(",", " ")
                          .replace("(", " ")
                          .replace(")", " "))
        parts = {}
        for sym, num in self._formula_token.findall(s):
            count = 1.0
            if num not in ("", "+", "-"):
                try:
                    count = float(num)
                except ValueError:
                    count = 1.0
            parts[sym] = parts.get(sym, 0.0) + count
        return parts
    
    def _phase_flag_to_int(self, ch):
        return 1 if (ch or "").strip().upper() == "G" else 0

    # ----- Record Line Parsing Methods ---------------------------------------------------------------------
    def _parse_line1(self, line):
        l1 = self._pad80(line)
        # NASA-thermo style: name in cols 1-24, comment in 25-80
        name    = l1[:24].strip()
        comment = l1[24:80].rstrip()
        return name, comment

    def _parse_line2(self, line):
        l2 = self._pad80(line)
        n_T_intervals_raw = l2[1:2]       # col 2
        n_T_intervals = int(n_T_intervals_raw) if n_T_intervals_raw.strip().isdigit() else None
        ident_code    = l2[3:9].strip()   # cols 4–9
        formula_raw   = l2[10:50].strip() # cols 11–50
        phase_flag    = l2[51:52].strip() # col 52
        mw            = self._safe_float(l2[52:65])  # cols 53–65
        hf298         = self._safe_float(l2[65:80])  # cols 66–80
        return {
            "num_T_intervals": n_T_intervals,
            "ident_code": ident_code,
            "formula_raw": formula_raw,
            "phase_flag": self._phase_flag_to_int(phase_flag),
            "mw": mw,
            "HF": hf298,
        }

    def _parse_interval_block(self, lines3to5):
        """
        lines3to5: [line3, line4, line5] for one temperature interval
        Returns a dict for that range.
        """
        # line 3
        l3 = self._pad80(lines3to5[0])
        t_range     = self._split_floats_field(l3[1:21])       # cols 2–21
        n_cp_coeff  = int(l3[22:23]) if l3[22:23].strip().isdigit() else None
        t_exponents = self._split_floats_field(l3[23:63])      # cols 24–63
        hf_line3    = self._safe_float(l3[65:80])              # cols 66–80
        t_low, t_high = (t_range + [None, None])[:2]

        # line 4: 5D16.8
        l4 = self._pad80(lines3to5[1])
        f4 = [l4[i:i+16] for i in range(0, 80, 16)]
        cp_1_5 = [self._safe_float(x) for x in f4]

        # line 5: 3D16.8 then 2D16.8
        l5 = self._pad80(lines3to5[2])
        f5_cp = [l5[i:i+16] for i in range(0, 48, 16)]
        f5_b  = [l5[i:i+16] for i in range(48, 80, 16)]
        cp_6_8 = [self._safe_float(x) for x in f5_cp]
        b_vals = [self._safe_float(x) for x in f5_b]
        b1 = b_vals[0] if len(b_vals) > 0 else None
        b2 = b_vals[1] if len(b_vals) > 1 else None

        cp_all = (cp_1_5 or []) + (cp_6_8 or [])
        if self.include_b_in_coeff:
            if b1 is not None: cp_all.append(b1)
            if b2 is not None: cp_all.append(b2)

        return {
            "key": self._range_key(t_low, t_high),
            "t_low": t_low,
            "t_high": t_high,
            "num_coeff": n_cp_coeff,
            "hf": hf_line3,
            "coeff": cp_all,
            "T_exp": t_exponents,
            "b1": b1,
            "b2": b2,
        }

    def _parse_species_record(self, record_lines):
        """
        record_lines: list[str] for ONE species.
        Returns {sp_key: {...}} or None if skip_empty and there are no intervals.
        """
        if len(record_lines) < 2:
            raise ValueError("Record too short to be a valid NASA-thermo species.")

        # line 1
        name, comment = self._parse_line1(record_lines[0])

        # line 2
        meta = self._parse_line2(record_lines[1])

        # Declared vs inferred intervals
        n_declared = meta["num_T_intervals"]
        n_inferred = max(0, (len(record_lines) - 2) // 3)

        if n_declared is None:
            n_intervals = n_inferred
        elif n_declared <= 0:
            # some files set 0 even when blocks exist — prefer the real blocks if present
            n_intervals = n_inferred if n_inferred > 0 else 0
        else:
            n_intervals = n_declared

        ranges = {}
        has_coeff = False

        if n_intervals > 0:
            base = 2
            for i in range(n_intervals):
                blk = record_lines[base + 3*i : base + 3*i + 3]
                if len(blk) < 3:
                    break
                interval = self._parse_interval_block(blk)
                ranges[interval["key"]] = {
                    "t_low": interval["t_low"],
                    "t_high": interval["t_high"],
                    "num_coeff": interval["num_coeff"],
                    "hf": interval["hf"],
                    "coeff": interval["coeff"],
                    "T_exp": interval["T_exp"],
                }
            has_coeff = len(ranges) > 0

        # Optionally skip truly empty species (no coefficient blocks at all)
        if not has_coeff and self.skip_empty:
            return None

        sp_key = (name.split() or [name])[0]
        return {
            sp_key: {
                "name": name,
                "comment": comment,
                "num_T_intervals": n_intervals,
                "ident_code": meta["ident_code"],
                "formula": self._parse_formula_dict(meta["formula_raw"]),
                "phase_flag": meta["phase_flag"],
                "mw": meta["mw"],
                "HF": meta["HF"],
                "ranges": ranges,
                "has_coeff": has_coeff,   # <- useful flag for downstream logic
            }
        }


    # ----- Public API: Reactants/Products -> Structured Dict -----------------------------------------------
    def parse_species_records(self, section='reactants'):
            """
            Parse all species in the chosen section ('reactants' or 'products') into
            the structured dict format.
            """
            if section not in ('reactants', 'products'):
                raise ValueError("section must be 'reactants' or 'products'")
    
            src = self.reactants if section == 'reactants' else self.products
            out = {}
            for name, rec_str in src.items():
                lines = rec_str.splitlines(True)  # keep newlines
                parsed = self._parse_species_record(lines)
                out.update(parsed)
            return out
        
    def build_database(self):
        """
        Parse BOTH reactants and products and combine into a single dict.
        If the same key appears in both, reactants will overwrite products (tweak if desired).
        """
        db = {}
        db.update(self.parse_species_records('products'))
        db.update(self.parse_species_records('reactants'))
        return db

if __name__ == '__main__':
    
    tc = ThermoDBConverter('raw/thermo.inp', include_b_in_coeff=True)
    reactant_db = tc.parse_species_records('reactants')
    product_db  = tc.parse_species_records('products')
    full_db     = tc.build_database()

    # quick sanity print
    print(f"reactants parsed: {len(reactant_db)}")
    print(f"products parsed:  {len(product_db)}")

    