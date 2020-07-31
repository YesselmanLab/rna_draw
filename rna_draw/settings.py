import os
import platform


def get_lib_path():
    file_path = os.path.realpath(__file__)
    spl = file_path.split("/")
    base_dir = "/".join(spl[:-2])
    return base_dir


def get_os():
    OS = None
    if platform.system() == "Linux":
        OS = "linux"
    elif platform.system() == "Darwin":
        OS = "osx"
    else:
        raise SystemError(platform.system() + " is not supported currently")
    return OS


class Paths:
    LIB_PATH = get_lib_path()
    RESOURCES_PATH = LIB_PATH + "/rna_draw/resources/"
    UNITTEST_PATH = LIB_PATH + "/tests/"
    EXAMPLE_PATH = LIB_PATH + "/examples/"
