
class SpeciesDB:
    def __init__(self):
        self.by_name    = {}            # Allows using the primary name to get information on the chemical
        self.alias      = {}            # chemical species can have multiple common names or abbriviations