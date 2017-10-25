"""Custom plug-ins to create interactive plots with `matplotlib` and
`d3.js` for exploring fits and material properties.
"""
from mpld3 import utils, plugins
from os import path

class PointDetailImage(object):
    """Container for an image associated with a single point in the main plot.

    .. note:: To add labels to x- and y- axes, use `subplot_kw`
      arguments. Adding a title uses that same set. See
      :meth:`~matplotlib.figure.Figure.add_subplot` for details.

    .. note:: If `figsize` is not included in the `fig_kw` it will be added as
      (3, 2) inches.

    Args:
        x (numpy.ndarray): of x-values used for the 2D plot.
        y (numpy.ndarray): of y-values used for the 2D plot.
        plot (str): type of plot to perform; one of the plot types on
          :class:`matplotlib.axes._subplots.AxesSubplot`.
        subplot_kw (dict): keywords passed to the
          :meth:`~matplotlib.figure.Figure.add_subplot` call used to create each
          subplot. 
        gridspec_kw (dict): keywords passed to the
          :class:`~matplotlib.gridspec.GridSpec` constructor used to create the
          grid the subplots are placed on.
        fig_kw (dict): keyword arguments are passed to the
          :func:`~matplotlib.pyplot.figure` call.
        plot_kw (dict): keyword arguments passed to the particular plot method.
        base64 (bool): when True, use a base-64 encoded `src` for the output of
          this image; otherwise, the image is saved to `folder` and only its
          name is returned.
        name (str): name of the point to display in the interactive plotter.
        imgtype (str): name of the image type in the interactive plotter.
        index (int): integer index of this plot in the parent collection.
        folder (str): path to the folder where the saved image should be stored.
        save_kw (dict): keyword arguments passed to
          :func:`~matplotlib.pyplot.savefig`.

    Attributes:
        base64 (bool): when True, use a base-64 encoded `src` for the output of
          this image; otherwise, the image is saved to `folder` and only its
          name is returned.
        name (str): name of the point to display in the interactive plotter.
        imgtype (str): name of the image type in the interactive plotter.
        index (int): integer index of this plot in the parent collection.
        folder (str): path to the folder where the saved image should be stored.
        filename (str): name of the file where the resulting plot was saved if
          `base64=False`.
        encoded (str): base-64 encoded representation of the image if it
          `base64=True`.
    """
    def __init__(self, x, y, plot="scatter", subplot_kw=None, gridspec_kw=None,
                 fig_kw=None, plot_kw=None, base64=False, name=None,
                 imgtype=None, index=None, folder=None, save_kw=None):
        import matplotlib.pyplot as plt
        self.x, self.y = x, y
        
        #Create the figure and axes for this plot.
        if fig_kw is None:
            fig_kw = {}
        if "figsize" not in fig_kw:
            fig_kw["figsize"] = (3, 2)
        figure, axes = plt.subplots(subplot_kw=subplot_kw,
                                    gridspec_kw=gridspec_kw,
                                    **fig_kw)

        #Grab the particular plot type that the user wants to make, then
        #generate the figure.
        if plot_kw is None:
            plot_kw = {}
        plotter = getattr(axes, plot)
        points = plotter(x, y, **plot_kw)

        self.folder = folder
        self.name = name
        self.imgtype = imgtype
        self.index = index
        self.base64 = base64
        self.filename, self.savepath = None, None
        self.encoded = None
        
        if save_kw is None:
            save_kw = {}
        if folder is None:
            return
            
        if not base64:
            self.filename = "{}-{}.png".format(name, index)
            savepath = path.join(folder, self.filename)
            plt.savefig(savepath, format="png", **save_kw)
        else:
            import StringIO
            import urllib, base64
            imgdata = StringIO.StringIO()
            plt.savefig(imgdata, format='png', **save_kw)
            imgdata.seek(0)
            b64 = urllib.quote(base64.b64encode(imgdata.buf))
            self.encoded = 'data:image/png;base64,{}'.format(b64)

    @property
    def url(self):
        """Returns the url to place in the `<img>` tag's `src` attribute. If
        :attr:`base64` is True, this will be the full base-64 encoded image. If
        False, then it will be the filename (relative to the interactive plot's
        directory).
        """
        if self.base64:
            return self.encoded
        else:
            return self.filename

def html(data, folder, plot="scatter", subplot_kw=None, gridspec_kw=None,
         fig_kw=None, plot_kw=None, titles=None, ncols=3, font=None):
    """Creates an HTML interactive plot package in the specified directory.

    Args:
        data (dict): keys are point names; values are `dict` with a `location`
          key specifying the 2D tuple `(x, y)` where the point should be
          plotted. The `index` key specifies its order in the plot. Additional
          keys of image types and
          :class:`PointDetailImage` values provide data for plotting auxiliary
          images.
        plot (str): type of plot to perform; one of the plot types on
          :class:`matplotlib.axes._subplots.AxesSubplot`.
        subplot_kw (dict): keywords passed to the
          :meth:`~matplotlib.figure.Figure.add_subplot` call used to create each
          subplot. 
        gridspec_kw (dict): keywords passed to the
          :class:`~matplotlib.gridspec.GridSpec` constructor used to create the
          grid the subplots are placed on.
        fig_kw (dict): keyword arguments are passed to the
          :func:`~matplotlib.pyplot.figure` call.
        plot_kw (dict): keyword arguments passed to the particular plot method.        
        titles (dict): keys are image types in `data`; values are `str` titles
          to display for each of the image types.
        ncols (int): number of columns for the grid of detail images.
        font (dict): key-value pairs for the :func:`matplotlib.rc` configuration
          for `font`.

    Examples:

    This is a simple example of what `data` might look like.

    .. code-block:: python

        data = {
            "gp_2b": {
                "location": (0.2, 1.4),
                "index": 0,
                "potato": PointDetailImage(),
                "EvsV": PointDetailImage()
            },
            "gp_3b": {
                "location": (0.82, 0.75),
                "index": 1,
                "potato": PointDetailImage(),
                "EvsV": PointDetailImage()
            }
        }
    """
    import matplotlib
    if font is None:
        font = {'size'   : 22}
    matplotlib.rc('font', **font)    

    #Sort the data by index.
    sdata = sorted(data.items(), key=(lambda d: d[1]["index"]))
    
    #Next, compile lists of URLs to use for each of the plot types.
    first = next(data.itervalues())
    imgtypes = [k for k in first.keys() if k not in ["location", "index"]]
    images = {k: [] for k in imgtypes}
    names, x, y = [], [], []    
    
    for name, detail in sdata:
        for img in imgtypes:
            images[img].append(detail[img].url)
            
        names.append(name)
        x.append(detail["location"][0])
        y.append(detail["location"][1])

    import matplotlib.pyplot as plt
    import numpy as np
    import mpld3
    if fig_kw is None:
        fig_kw = {}
    figure, axes = plt.subplots(subplot_kw=subplot_kw,
                                gridspec_kw=gridspec_kw,
                                **fig_kw)
    
    #Grab the particular plot type that the user wants to make, then
    #generate the figure. We need to compile the x and y arrays from the
    #location information in each point.
    if plot_kw is None:# pragma: no cover
        plot_kw = {}
    if "s" not in plot_kw:
        plot_kw["s"] = 300
    plotter = getattr(axes, plot)
    points = plotter(np.array(x), np.array(y), **plot_kw)
            
    from matdb.plotting.matd3 import ImagesAtPoint
    plugiap = ImagesAtPoint(points, names, ncols, titles, **images)
    mpld3.plugins.connect(figure, plugiap)

    #Next, copy the JS and CSS dependencies over from the package data into the
    #directory.
    from matdb.utility import relpath, copyonce
    folder = path.abspath(path.expanduser(folder))
    js = relpath("matdb/js/imagepoint.js")
    copyonce(js, path.join(folder, "imagepoint.js"))
    css = relpath("matdb/css/imagepoint.css")
    copyonce(css, path.join(folder, "imagepoint.css"))    

    target = path.join(folder, "index.html")
    with open(target, 'w') as f:
        mpld3.save_html(figure, f)

class ImagesAtPoint(plugins.PluginBase):
    """Plugin for displaying a series of rasterized plots specific to
    the particular point that get's clicked on.

    Args:
        points (matplotlib.collections.PathCollection): list of points returned
          by the matplotlib plotting routine (such as
          :func:`matplotlib.pyplot.plot`, :func:`matplotlib.pyplot.scatter`,
          etc.
        names (list): of names or titles to display for each point in `points`.
        ncols (int): number of columns of images below the main image.
        titles (dict): keys are user-defined plot types in `images`; values are
          `str` titles to use for each of the image types.
        images (dict): keys are user-defined plot types;
          values are a list of URLs to put in the `<img>` tag's `src` attribute
          for each point in `points`.

    Attributes:
        names (list): of names or titles to display for each point in `points`.
        ncols (int): number of columns of images below the main image.
        images (dict): keys are user-defined plot types;
          values are a list of base64-encode URLs (returned by
          :func:`to_base64`) to display for each point in `points`.
    """
    
    JAVASCRIPT = """
    mpld3.register_plugin("image@point", ImagesAtPoint);
    ImagesAtPoint.prototype = Object.create(mpld3.Plugin.prototype);
    ImagesAtPoint.prototype.constructor = ImagesAtPoint;
    ImagesAtPoint.prototype.requiredProps = ["id", "images", "settings"];
    
    function ImagesAtPoint(fig, props){
        mpld3.Plugin.call(this, fig, props);
    };
    
    ImagesAtPoint.prototype.draw = function() {
	var iap = this;
        mpld3_load_lib("imagepoint.js", function() {
           imagepoint(iap);
        });
    };

    """
    settings = ["ncols"]
    def __init__(self, points, names=None, ncols=3, titles=None, **images):
        self.dict_ = {"type": "image@point",
                      "id": utils.get_id(points),
                      "images": images,
                      "settings": {}}

        self.ncols = ncols
        self.names = names
        self.npoints = len(next(images.itervalues()))
        if self.names is None:# pragma: no cover
            self.names = ["Point {}".format(i) for i in range(self.npoints)]
        self.dict_["settings"]["names"] = self.names

        self.titles = titles
        if self.titles is None:# pragma: no cover
            self.titles = {k: k for k in images}
        self.dict_["settings"]["titles"] = self.titles
        
        for attr in self.settings:
            self.dict_["settings"][attr] = getattr(self, attr)
