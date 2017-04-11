"""Exposes classes and functions for interacting with the database
folders via a simple configuration file.
"""
from os import path
from matdb import msg
class DatabaseCollection(object):
    """Represents a collection of databases (all inheriting from
    :class:`Database`) that are all related be the atomic species that they
    model.

    .. note:: See the list of attributes below which correspond to the sections
    in the YAML database specification file.

    Args:
        name (str): name of the configuration that this database collection is
          operating for.
        poscar (str): name of the POSCAR file in `root` to extract atomic
          configuration information from.
        root (str): root directory in which all other database collections for
          the configurations in the same specification will be stored.
        parent (Controller): instance controlling multiple configurations.
        database (dict): with keys and values describing the kinds of
          sub-configuration databases to setup.

    Attributes:
        title (str): title for the alloy system that these databases work with.
        species (list): of `str` element names in the alloy system.
        potcars (dict): keys are lower-case element names; values are *suffixes*
          for the various pseudo-potentials that ship with VASP.
        root (str): root directory in which all other database directories will
          be stored. Defaults to the current directory.
        incar (dict): keys are valid settings in the INCAR file that should be
          used as defaults for *all* database calculations (with their
          corresponding values).
        kpoints (dict): keys are valid settings for the Mueller k-point PRECALC
          file that should be used as default for *all* database calculations.
        databases (dict): keys are database types (e.g. `liquid`, `phonon`,
          etc.); values are the corresponding class instances.
        order (list): of string database names (attribute `name` of each
          database class) in the order in which they should be processed.
        parent (Controller): instance controlling multiple configurations.
    """
    order = ["phonbase", "phoncalib", "phonons"]
    
    def __init__(self, name, poscar, root, parent, database):
        from quippy.atoms import Atoms
        self.name = name
        self.atoms = Atoms(path.join(root, poscar), format="POSCAR")
        self.root = path.join(root, name)

        if not path.isdir(self.root):
            from os import mkdir
            mkdir(self.root)

        parrefs = ["species", "potcars", "incar", "kpoints", "execution"]
        for ref in parrefs:
            setattr(self, ref, getattr(parent, ref))
        self.parent = parent
        
        from importlib import import_module
        self.databases = {}
        for dbspec in database:
            modname, clsname = dbspec["type"].split('.')
            fqdn = "matdb.database.{}".format(modname)
            module = import_module(fqdn)
            if not hasattr(module, clsname):
                #We haven't implemented this database type yet, just skip the
                #initialization for now.
                continue
            
            cls = getattr(module, clsname)

            #Make a copy of the original dictionary so that we don't mess up the
            #pointers; then add in the keyword arguments that are missing.
            cpspec = dbspec.copy()
            del cpspec["type"]
            cpspec["atoms"] = self.atoms
            cpspec["root"] = self.root
            cpspec["parent"] = self
            instance = cls(**cpspec)
            self.databases[instance.name] = instance

    def execute(self):
        """Submits job array files for any of the databases that are ready to
        execute, but which haven't been submitted yet.
        """
        ready = True
        for dbname in self.order:
            if dbname not in self.databases:
                continue
            if ready:
                ready = (self.databases[dbname].cleanup() or
                         self.databases[dbname].execute())
            else:
                imsg = "Database {} is not ready to execute yet. Done."
                msg.info(imsg.format(dbname))
                break        
            
    def cleanup(self):
        """Runs the cleanup methods of each database in the collection, in the
        correct order.
        """
        ready = True
        for dbname in self.order:
            if dbname not in self.databases:
                continue
            if ready:
                ready = self.databases[dbname].cleanup()
            else:
                msg.info("Database {} is not ready yet. Done.".format(dbname))
                break
            
    def setup(self):
        """Sets up the database collection by generating the POTCAR file and
        initializing any databases that haven't already been initialized.

        .. note:: The db setup functions are setup to only execute once, and then
           only if their dependencies have completed their calculations. This
           method can, therefore, be safely called repeatedly between different
           terminal sessions.
        """
        for dbname, db in self.databases.items():
            msg.info("Setting up database {}".format(dbname))
            db.setup()
        
class Controller(object):
    """Implements methods for tying a configuration dictionary (in
    YAML format) to instances of various databases.

    Args:
        config (str): path to the configuration YAML file that
          specifies all information for constructing the set of databases.

    Attributes:
        specs (dict): the parsed settings from the YAML configuration file.
        collections (dict): keys are configuration names listed in attribute
          `configs` of the YAML file. Values are the :class:`DatabaseCollection`
          instances.
    """
    def __init__(self, config):
        import yaml
        with open(config, 'r') as stream:
            self.specs = yaml.load(stream)

        self.root = path.abspath(path.expanduser(self.specs["root"]))
        self.title = self.specs["title"]
        self.collections = {}
        self.species = [s.lower() for s in self.specs["species"]]
        self.potcars = self.specs["potcars"]
        self.incar = self.specs.get("incar", {})
        self.kpoints = self.specs.get("kpoints", {})
        self.execution = self.specs["execution"]
        
        for cspec in self.specs["configs"]:
            name, poscar = cspec["name"], cspec["poscar"]
            instance = DatabaseCollection(name, poscar, self.root, self,
                                          self.specs.get("databases", []))
            self.collections[name] = instance

        #Extract the POTCAR that all the databases are going to use. TODO: pure
        #elements can't use this POTCAR, so we have to copy the single POTCAR
        #directly for those databases; update the create().
        self.POTCAR()

    def setup(self):
        """Sets up each of configuration's databases.
        """
        for name, coll in self.collections.items():
            coll.setup()

    def cleanup(self):
        """Runs cleanup on each of the configuration's databases.
        """
        for name, coll in self.collections.items():
            coll.cleanup()
        
    def POTCAR(self):
        """Creates the POTCAR file using the pseudopotential and version
        specified in the file.
        """
        target = path.join(self.root, "POTCAR")
        if not path.isfile(target):
            # Make sure that the POTCAR version and pseudopotential type match
            # up so that we don't get nasty surprises.
            potsrc = path.join(self.potcars["directory"],
                               "pot{}".format(self.potcars["pseudo"]))
            potcars = []
            
            for element in self.species:
                lel = element.lower()
                pp = element
                if lel in self.potcars:
                    pp += self.potcars[lel]

                ppot = path.join(potsrc, pp, "POTCAR")
                with open(ppot) as f:
                    first = f.readline()
                    if "version" in self.potcars:
                        assert self.potcars["version"] in first

                potcars.append(ppot)

            if len(potcars) == len(self.species):
                from matdb.utility import cat
                cat(potcars, target)
            else:
                msg.err("Couldn't create POTCAR for system.")
