class RcParams(dict):
    # validate = { key : converter
    #            for key, (default, converter) in defaultParams.items()}
    def __init__(self, *args, **kwargs):
        self.update(*args, **kwargs)

    def __setitem__(self, key, val):
        # try:
        #    cval = self.validate[key](val)
        # except ValueError as ve:
        #    raise ValueError(f"Key {key} : {ve}") from None
        dict.__setitem__(self, key, val)

    def __getitem__(self, key):
        return dict.__getitem__(self, key)


defaultParams = {
    # default system paths
    'paths.igor_prefix': None,  # "/home/alfaceor/.local",
    'paths.igor_exec': None,  # "/home/alfaceor/.local/bin/igor",
    'paths.igor_data': None,  # "",
    'paths.igor_models': None,  # "/home/alfaceor/.local/share/igor/models",
    'paths.igor_src': None  # ""
}


# TODO: READ CONFIGURATION FROM FILE
import appdirs
import os

def create_config_files():
    try:
        import pathlib
        import pkg_resources
        import json

        dirs = appdirs.AppDirs('pygor3', 'alfaceor') #, version=__version__) # FIXME I'M USING SCM for versioning so?
        fln_config_json = 'config.json'

        default_config_path = os.path.join(dirs.user_data_dir, fln_config_json)

        pathlib.Path(dirs.user_data_dir).mkdir(parents=True, exist_ok=True)

        # TODO: if default_config_path exist just load it, else create it first
        if pathlib.Path(default_config_path).is_file():
            # asdf
            pass
        else:
            # pkg_resources.resource_filename('pygor3', fln_config_json)
            with open(pkg_resources.resource_filename('pygor3', fln_config_json), 'r') as f:
                json_object = json.load(f)

            # TODO: LOOOK FOR LOCATIONS IN SYSTEM
            # import copy
            # dictParams = copy.deepcopy(defaultParams)
            try:
                import subprocess
                p1 = subprocess.run(["which", "igor"], capture_output=True, text=True)
                igor_exec_path = p1.stdout.replace('\n', '')
                json_object['paths.igor_exec'] = igor_exec_path

                p2 = subprocess.run([igor_exec_path, "-getdatadir"], capture_output=True, text=True)
                json_object['paths.igor_data'] = p2.stdout.replace('\n', '')

            except Exception as e:
                print(e)
                pass

            with open(default_config_path, 'w') as f:
                json.dump(json_object, f, indent=4)

    except Exception as e:
        raise e

def load_config_files():
    try:
        dirs = appdirs.AppDirs('pygor3', 'alfaceor')
        default_config_path = os.path.join(dirs.user_data_dir, 'config.json')
        with open(default_config_path, 'r') as f:
            import json
            user_config_Params = json.load(f)
        return RcParams(user_config_Params)
    except Exception as e:
        raise e

rcParams = RcParams(defaultParams)

try:
    create_config_files()
    rcParams = load_config_files()
except Exception as e:
    print(e)
    pass


# TODO:
#  - CHECK IF AppDirs creates the directories
#  - LOAD data from files at the initial loading
#  - (optional) WRITE with AppDirs default locations.
# try:
#     os.makedirs(directory)
# except OSError as e:
#     sys.exit(_('Failed to create directory: {}').format(e))
#
# # FIXME CREATE A DIRECTORY
# # Create all necessary directories.
# for path in [args.log, args.token, args.config]:
#     directory = os.path.dirname(path)
#     if directory and not os.path.isdir(directory):
#         try:
#             os.makedirs(directory)
#         except OSError as e:
#             sys.exit(_('Failed to create directory: {}').format(e))



