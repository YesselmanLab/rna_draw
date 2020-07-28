import matplotlib.cm

class Data(object):
    def __init__(self, data_str=None, data_file=None, palette=None,
                 vmin=None, vmax=None, ignore_restype=None):
        self.vals = self.__parse_vals(data_str, data_file)
        self.min = min(self.vals)
        self.max = max(self.vals)
        self.ignore_restype = []
        self.palette = palette
        if vmin is not None:
            self.min = vmin
        if vmax is not None:
            self.max = vmax
        if ignore_restype is not None:
            self.ignore_restype = list(ignore_restype)
        if self.palette is None:
            self.palette = matplotlib.cm.get_cmap("Reds")

    def __parse_vals(self, data_str, data_file):
        if data_str is None and data_file is None:
            raise ValueError("must supply a str or file")

        if data_str:
            return [float(x) for x in data_str.split(";")]

        else:
            f = open(data_file)
            lines = f.readlines()
            f.close()

            return [float(x.rstrip()) for x in lines]



