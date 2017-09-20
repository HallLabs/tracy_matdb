"""Custom plug-ins to create interactive plots with `matplotlib` and
`d3.js` for exploring fits and material properties.
"""
from mpld3 import utils, plugins

def to_base64(fig):
    """Saves the plot drawn on the specified figure to a base64 string and returns
    the URI that could be placed in an <img> tag.

    Args:
        fig (matplotlib.figure.Figure): figure that plots have been drawn to so
          far.
    """
    import StringIO
    import urllib, base64
    imgdata = StringIO.StringIO()
    fig.savefig(imgdata, format='png')
    imgdata.seek(0)
    b64 = urllib.quote(base64.b64encode(imgdata.buf))
    return 'data:image/png;base64,{}'.format(b64)

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
        images (dict): keys are user-defined plot types;
          values are a list of base64-encode URLs (returned by
          :func:`to_base64`) to display for each point in `points`.

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
    def __init__(self, points, names=None, ncols=3, **images):
        self.dict_ = {"type": "image@point",
                      "id": utils.get_id(points),
                      "images": images,
                      "settings": {}}

        self.ncols = ncols
        self.names = names
        self.npoints = len(next(images.itervalues()))
        if self.names is None:
            self.names = ["Point {}".format(i) for i in range(self.npoints)]
        self.dict_["settings"]["names"] = self.names
        
        for attr in self.settings:
            self.dict_["settings"][attr] = getattr(self, attr)
