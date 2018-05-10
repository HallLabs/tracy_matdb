# How to Contribute a Database to `matdb`.

#### Table of Contents

[Introduction](#introduction)

[Source Code](#source-code)
  * [Group Name](#group-name-and-default-name)
  * [Initialization](#initialization)
  * [Fitting Configs](#fitting-configs)
  * [Dictionary](#dictionary)
  * [Result Set](#result-set)
  * [Ready](#ready)
  * [Additional Methods](#additional-methods)

## Introduction

Inside `matdb` the Database types are organized in groups. Each group
has a unique name and is a subclass of a generic `Group` class. To
create a new group it is necessary to make a new module in the
`matdb/database/` package with the desired group name. For example if
you wanted to create a group to calculate and store the phonon
properties of a set of configurations then you would create the module
`matdb/database/phonon.py` which would contain [code](#source-code) for
the `Phonon` group.

## Source Code

Below is the source code that you will need to modify and extend in
order to create your group. The code is filled with place holders that
you will need to replace with code that matches you're use
case. They are indicated by `<<` and `>>` tokens.

Examples of completed groups that make good references for
developers are [manual group](../matdb/matdb/database/simple.py), [enumeration
group](../matdb/matdb/database/enumeration.py) and the [hessian
group](../matdb/matdb/database/hessian.py).

Details for what each replacement should be are given below. If your
group can generate configurations that would be used only as seed
configurations for other groups, i.e., no calculations need to be
performed on your group, then use the code relating to the
`extractable` key word that`s surrounded in <<Optional>> <<End
Optional>> flags.

```
class <<Group_Name>>(Group):
    """<<Your description here>>.
    Args:
        name (str): the database name.
	<<Additional arguments your database group will require>>.
	<<Optional>>
        extractable (bool): True if a calculation is to be performed.
	<<End Optional>>
	dbargs (dict): dictionary of arguments to be passed to the 
	  `Group` class.
	  
    .. note:: Additional attributes are also exposed by the super class
      :class:`Group`.

    Attributes:
        name (str): name of this database type relative to the over database
          collection. This is also the name of the folder in which all of its
          calculations will be performed.
	<<Additional attributes your database group will have>>.
    """

    <<All groups must have an __init__ function>>.
    def __init__(self, name=<<your groups default name>>, <<your group args>>, 
    	        <<Optional>>extractable=True,<<End Optional>> **dbargs):
        self.name = name
        self.seeded = <<True if your group uses a seed configuration, False otherwise>>
	<<Optional>>
	self.extractable = extractable
	<<End Optional>>
        dbargs["prefix"] = <<your group prefix>>
        dbargs["cls"] = <<Group_Name>>
        if <<Group_Name>> not in dbargs['root']:
            from os import mkdir
            new_root =path.join(dbargs['root'],<<Group_Name>>)
            if not path.isdir(new_root):
                mkdir(new_root)
            dbargs['root'] = new_root
        super(<<Group_Name>>, self).__init__(**dbargs)
        
        #Make sure that we override the global calculator default values with
        #those settings that we know are needed for good phonon calculations.
        calcargs = self.database.calculator.copy()
        if "calculator" in dbargs:
            if dbargs["calculator"] is not None and "name" in dbargs["calculator"]:
                calcargs.update(dbargs["calculator"])
                self._set_calc_defaults(calcargs)
                dbargs["calculator"] = calcargs            

	<<Whatever else your group needs to do during the initialization.>>
		

    @property
    def fitting_configs(self):
        """Returns a :class:`matdb.atoms.AtomsList` for all configs in this
        group. 
        """

        if len(self.sequence) == 0:
            <<Some method that returns the configurations from the group
	    that would get used by the fitters.>>
        else:
            result = []
            for g in self.sequence.values():
                result.extend(g.fitting_configs)
            return result
                
    def sub_dict(self):
        """Returns a dict needed to initialize the class.
        """
        args = <<dictionary of the parameters used in the initialization>>
        return args

    @property
    def rset(self):
        """
        Returns:
            list: of :class:`matdb.atoms.Atoms`
        """
        if len(self.sequence) == 0:
            #We are at the bottom of the stack; 
	    <<find the rset for this level of the database>>
            return [self.atoms]
        else:
            #Check where we are in the stack. If we are just below the database,
            #then we want to return <<your description of the rset here>>
	    #If we are not, then we must a parameter grid of sequences
            #to select from.
            result = []
            for g in self.sequence.values():
                result.extend(g.rset)
	    return AtomsList(result)

    def ready(self):
        """Returns True if all the calculations have been completed.
        """
        self._expand_sequence()
        if len(self.sequence) == 0:
	    <<code to determine if calculations were complete.>>
            if not result:
                msg.std("{} is not ready. Exiting.".format(self.root), 2)
            return result
        else:
            ready = False
            for p in self.sequence.values():
                if not p.ready():
                    msg.std("{} is not ready. Exiting.".format(p.root), 2)
                    break
            else:
                ready = True
            return ready

    def setup(self, rerun=False):
        """<<explanation of setup function.>>

        Args:
            rerun (bool): when True, recreate the folders even if they
              already exist. 
        """        
        super(<<Group_Name>>, self).setup(self._setup_configs, rerun)
            
    def _setup_configs(self, rerun):
        """<<explanation of function.>>

        Args:
            rerun (bool): when True, recreate the folders even if they
              already exist. 
        """        
        #We also don't want to setup again if we have the results already.
        if self.ready():
            return

        if not self.is_setup():
	    <<code to setup the group's configurations. This code must
	    make use of the self.crate(atoms_obj<<Optional>>,
	    extractable=self.extractable<<End Optional>>) method implemnted in
	    the maine `Group` class.>>

        # Last of all, create the job file to execute the job array.
        self.jobfile(rerun)

    <<Optional>>
    def can_extract(self):
        """Runs post-execution routines to clean-up the calculations. This super class
        implementation only checks that each of the sub-config directories has
        the necessary files needed for extraction of output.
        """
        if not self.extractable:
            return self.is_setup()
        else:
            return super(Manual, self).can_extract()
    <<End Optional>>
    <<Any additional methods, or properties your group will need.>>
```

### Group Name and default name

You need to choose a name for your group that follows python best
practices, i.e., an uppercase name that is `_` separated. This will
replace Group_Name in many places in the code. Additionally in the
`__init__` method has a name argument that you need to assign. This
will be the name of the group in the database file structure, by
default it should match the `Group_Name` except be lower case.

### Initialization

Each group must be initialized appropriately. In the default
`__init__` function you need to replace <<Group_Name>>, everywhere it
appears. You will also need to select a prefix for your group, the
prefix is usually 1-2 camel case letters that match the Group_Name,
for example the `Hessian` group has a prefix of H and the
`Enumeration` group has a prefix of E. Other than that you need to
determine the input arguments and set the attributes of the class to
meet the needs of your group (document both in the appropriate
places).

### Fitting Configs

Each group must have a `fitting_configs` property. This method selects
which configurations from the group will be used to train and test
potentials later by `matdb`. This can be all the configurations or a
subset depending on the nature of the configurations and the type of
calculations being performed.

### Dictionary

Each group must have a method `sub_dict` that converts the input
arguments to a dictionary as output. Thus ensuring that the `Group`
instance can be recreated.

### Result Set

Each group must have an `rset` property. This is a set of
configurations that will be used by the next group in the database for
it's configurations, for example the `Aflow` group's `rset` is
everything it finds, if the `Hessian` group follows the `Aflow` group
in the database then the `Hessian` group can use the `Aflow` group's
rset as it's seed configurations. The `rset` that exists for the
database you build needs to be the configurations that will get passed
on to the next set of calculations, i.e., the best configurations.

### Ready

Each group must have a `ready` method. This method determines if the
groups calculations have been completed so that the group can be
processed for use in fitting and other operations.

### Setup

Each group must have a `setup` method. The only changes you need to
make to the one provided are to update the documentation and to
replace <<Group_Name>> with your groups class name.

Each group must also have a `_setup_configs` method. This method will
perform the necessary steps to setup an `atoms` object for each
configuration in the group that calculations can be performed on. Once
the code has successfully created an `atoms` object it needs to call
the `self.create(atoms_obj)` function in order to have `matdb` create
the directory the calculations will be performed in and save the
`atoms` object to file.

### Additional Methods

Your group may need additional methods that will help it
function. 
