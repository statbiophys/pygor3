#      Author: Quentin Marcou & Carlos Olivares
#
#  This source code is distributed as part of the IGoR software.
#  IGoR (Inference and Generation of Repertoires) is a versatile software to
#  analyze and model immune receptors generation, selection, mutation and all
#  other processes.
#   Copyright (C) 2017  Quentin Marcou
#
#   This program is free software: you can redistribute it and/or modify
#   it under the terms of the GNU General Public License as published by
#   the Free Software Foundation, either version 3 of the License, or
#   (at your option) any later version.
#
#   This program is distributed in the hope that it will be useful,
#   but WITHOUT ANY WARRANTY; without even the implied warranty of
#   MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
#   GNU General Public License for more details.

#   You should have received a copy of the GNU General Public License
#   along with this program.  If not, see <https://www.gnu.org/licenses/>.

###   from . import IgorAligns
###   from . import aligns
###   from . import counters
###   from . import models
###   from . import utils
###   
###   __all__ = ["aligns","counters","models","utils"]


###   from . import IgorAlignment_data
###   from . import IgorModel
###   from . import IgorBestScenarios
###   from . import IgorIndexedSequencesDB
###   from . import IgorSqliteDBBestScenarios
###   from . import IgorSqliteDB


###   import IgorAlignment_data
###   import IgorBestScenarios
###   import IgorIndexedSequencesDB
###   import IgorModel
###   import IgorSqliteDBBestScenarios
###   import IgorSqliteDB
###   __all__ = ["IgorAlignment_data","IgorBestScenarios","IgorIndexedSequencesDB","IgorModel", "IgorSqliteDBBestScenarios", "IgorSqliteDB"]

# __version__ =

class RcParams(dict):
    #validate = { key : converter
    #            for key, (default, converter) in defaultParams.items()}
    def __init__(self, *args, **kwargs):
        self.update(*args, **kwargs)

    def __setitem__(self, key, val):
        #try:
        #    cval = self.validate[key](val)
        #except ValueError as ve:
        #    raise ValueError(f"Key {key} : {ve}") from None
        dict.__setitem__(self, key, val)
    
    def __getitem__(self, key):
        return dict.__getitem__(self, key)

defaultParams = {
    # default system paths
    'paths.igor_prefix' : None, #"/home/alfaceor/.local",
    'paths.igor_exec' : None, #"/home/alfaceor/.local/bin/igor",
    'paths.igor_data' : None, #"",
    'paths.igor_models': None, #"/home/alfaceor/.local/share/igor/models",
    'paths.igor_src' : None #""
}

rcParams = RcParams(defaultParams)

import pkg_resources

# __name__ in case you're within the package
# - otherwise it would be 'lidtk' in this example as it is the package name
sql_path = 'IgorDB.sql'  #
sql_filepath = pkg_resources.resource_filename(__name__, sql_path)

#from .IgorModel import *
from .IgorIO import *
from .IgorSqliteDB import *
# from .IgorSqliteDBBestScenarios import *
# from .IgorBestScenarios import *
from . import imgt

from .AIRR import *

#__all__ = ["IgorModel"]


