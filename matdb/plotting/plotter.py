#Copyright (C) 2019  HALL LABS
#
#This program is free software: you can redistribute it and/or modify
#it under the terms of the GNU General Public License as published by
#the Free Software Foundation, either version 3 of the License, or
#(at your option) any later version.
#
#This program is distributed in the hope that it will be useful,
#but WITHOUT ANY WARRANTY; without even the implied warranty of
#MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
#GNU General Public License for more details.
#
#You should have received a copy of the GNU General Public License
#along with this program.  If not, see <https://www.gnu.org/licenses/>.
#
#If you have any questions contact: wmorgan@tracy.com
"""This module exposes a generic plotter that can generically plot quantities
against one another. Note that the actual plotting is carried out using the
:class:`~matdb.plotting.matd3.PointDetailImage`. This class wraps that one to
provide access to configurations from the global `matdb.yml` file, which can
include file directives to other folders.
"""
from matdb.plotting.matd3 import PointDetailImage as PDI

class PlotManager(object):
    """Represents a collection of plots that can be made generically for
    multiple objects.
    """
    def __init__(self, controller):
        self.controller = controller
        self.types = controller.specs.get("plotting")
        self.configs = {}
        for tdict in self.types:
            cdict = tdict.copy()
            name = cdict["name"]
            self.configs[name] = cdict

    def plot(self, objs, name):
        """Prepares a plot using the specified list of objects and a plotting
        configuration.

        Args:
            objs (list): list of trainer or database objects that support the methods
              described in the plotting configuration.
            name (str): name of a plotting configuration in :attr:`configs`.
        """
        config = self.configs[name].copy()
        method = config.pop("method")
        xc, yc = config.pop("x"), config.pop("y")
        plotargs = config.get("plotargs", {})
        x, y, label = [], [], []

        subplot_kw = {}
        subplot_kw.update({'x'+k: v for k, v in xc.items() if k != "data"})
        subplot_kw.update({'y'+k: v for k, v in yc.items() if k != "data"})
        if "subplot_kw" in plotargs:
            plotargs["subplot_kw"].update(subplot_kw)
        else:
            plotargs["subplot_kw"] = subplot_kw

        if "savefmt" not in plotargs:
            plotargs["savefmt"] = "pdf"

        for obj in objs:
            caller = getattr(obj, method)
            data = caller(**config)
            x.append(data[xc["data"]])
            y.append(data[yc["data"]])
            label.append(obj.fqn)

            plot = PDI(x, y, labels=label, folder=self.controller.plotdir,
                       index=name, name=obj.fqn, **plotargs)
