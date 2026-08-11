import json
from importlib import resources
from pathlib import Path

try:
    import jsonschema
except Exception:
    jsonschema = None       # Validation Optional

class SpeciesDB:
    def __init__(self, autoload=True, data_pkg='thermoprop.data', validate=True):
        self.by_name            = {}                # Canonical -> Record (Allows Primary Name to get species record)
        self.alias              = {}                # Alias -> Canonical (Allows Aliases to get to canonical name)
        self.paths              = {}                # Canonical -> Source File Path (handy for debug/reload)
        self._schema            = None              # JSON schema for species data files
        self._validate_enabled  = bool(validate)    # Flag to determine if we validate the json files during loading
        if autoload:
            self.load_package_data(data_pkg)
 
    # ----- JSON Schema Loading -----------------------------------------------------------------------------
    def _load_schema(self):
        if self._schema is not None or not self._validate_enabled or jsonschema is None:
            return
        try:
            schema_pkg  = 'thermoprop.schema'
            schema_file = 'species.schema.json'
            with resources.as_file(resources.files(schema_pkg) / schema_file) as p:
                text = Path(p).read_text(encoding='utf-8')  # <-- always read text first
            try:
                self._schema = json.loads(text)
            except json.JSONDecodeError as e:
                lines = text.splitlines()
                bad = lines[e.lineno-1] if 0 < e.lineno <= len(lines) else ''
                pointer = ' ' * (e.colno-1) + '^'
                print(f'[speciesdatabase] schema JSON error at line {e.lineno}, col {e.colno}: {e.msg}\n{bad}\n{pointer}')
                self._schema = None
                self._validate_enabled = False
        except Exception as e:
            print(f'[speciesdatabase] schema load skipped: {e}')
            self._schema = None
            self._validate_enabled = False


    def _validate_record(self, rec, path):
        
        if not self._validate_enabled or jsonschema is None:
            return
        self._load_schema()
        
        if not self._schema:
            return
        
        try:
            jsonschema.validate(instance=rec, schema=self._schema)
        except jsonschema.exceptions.ValidationError as e:
            raise ValueError(f'schema validation failed for {path}:\n{e.message}') from e
        # extra domain checks: ordered, non-overlapping ranges
        ranges = sorted(rec['thermo_ranges'], key=lambda r: r['t_low'])
        for i, r in enumerate(ranges):
            if not (r['t_low'] < r['t_high']):
                raise ValueError(f'{path}: t_low >= t_high in range {i}')
            if i > 0:
                prev = ranges[i-1]
                if r['t_low'] < prev['t_high']:
                    raise ValueError(f'{path}: overlapping ranges {i-1} and {i}')

    # ----- Public Loaders ----------------------------------------------------------------------------------    
    def load_package_data(self, directory, package='thermoprop.data'):
        try:
            root = resources.files(package)
        except Exception as e:
            raise RuntimeError(f'Could not access package data at {package!r}: {e}') from e
        
        for entry in root.iterdir():
            if not entry.name.lower().endswith('.json'):
                continue
            with resources.as_file(entry) as fp:
                self._load_and_add(Path(fp))

    # ----- Core Helpers ------------------------------------------------------------------------------------
    def _load_and_add(self, path):
        # parse JSON with a helpful error if it’s malformed
        try:
            text = Path(path).read_text(encoding='utf-8')
            rec = json.loads(text)
        except json.JSONDecodeError as e:
            lines = text.splitlines()
            bad = lines[e.lineno-1] if 0 < e.lineno <= len(lines) else ''
            pointer = ' ' * (e.colno-1) + '^'
            raise ValueError(
                f'JSON parse error in {path} at line {e.lineno}, col {e.colno}: {e.msg}\n{bad}\n{pointer}'
            ) from e

        self._validate_record(rec, path)

        self._add_record(rec, str(path))
        
    def _add_record(self, rec, path=None):
        if 'species' not in rec:
            raise ValueError(f'species record missing "species": {path}')
        canonical = str(rec['species']).upper()
        self.by_name[canonical] = rec
        if path:
            self.paths[canonical] = path
        self.alias[canonical.lower()] = canonical
        for a in rec.get('aliases', []):
            self.alias[str(a).strip().lower()] = canonical

    # ----- API ---------------------------------------------------------------------------------------------
    def resolve(self, token):
        key = str(token).strip().lower()
        return self.alias.get(key, str(token).upper())

    def get(self, canonical):
        c = str(canonical).upper()
        if c not in self.by_name:
            raise ValueError(f'unknown species: {canonical!r}')
        return self.by_name[c]

    def list_species(self):
        """List all loaded species by canonical name and aliases."""
        if not self.by_name:
            print('(no species loaded)')
            return []
    
        # Build a simple table: name | aliases
        rows = []
        for name, record in sorted(self.by_name.items()):
            aliases = ', '.join(record.get('aliases', []))
            rows.append((name, aliases))
    
        # Show count
        print('\n')
        print(f'{len(rows)} species loaded:\n')
    
        # Figure out column widths
        name_w = max(len(r[0]) for r in rows)
        alias_w = max(len(r[1]) for r in rows) if any(r[1] for r in rows) else 0
    
        # Print header
        header = f'{"Species":<{name_w}}  {"Aliases":<{alias_w}}'
        print(header)
        print('-' * len(header))
    
        # Print rows
        for name, aliases in rows:
            print(f'{name:<{name_w}}  {aliases:<{alias_w}}')
    
        return [name for name, _ in rows]

    def reload(self, canonical=None):
        if canonical is None:
            items = list(self.paths.items())
        else:
            c = str(canonical).upper()
            items = [(c, self.paths[c])]
        for c,_ in items:
            if c in self.by_name:
                del self.by_name[c]
            for k,v in list(self.alias.items()):
                if v == c:
                    del self.alias[k]
        for c,p in items:
            self._load_and_add(Path(p))