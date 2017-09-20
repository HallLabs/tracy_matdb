function imagepoint(iap) {
    var jsloads = ['https://cdnjs.cloudflare.com/ajax/libs/salvattore/1.0.9/salvattore.min.js',
		   'https://code.jquery.com/jquery-3.2.1.min.js'];
    var cssloads = ['https://netdna.bootstrapcdn.com/bootstrap/3.1.1/css/bootstrap.min.css',
		    'imagepoint.css'];

    var MATDB_PANEL = "<div class='panel panel-primary col-md-4'><div class='panel-heading'>{{imgtype}}</div><div class='panel-body'>{{content}}</div></div>";

    //Define a function to recursively load javascript and then finally call the
    //last callback when everything is ready.
    var recloadjs = function(index, callback) {
	url = jsloads[index];
	if (index == jsloads.length - 1) {
	    mpld3_load_lib(url, callback);
	} else {
	    mpld3_load_lib(url, function() {
		recloadjs(index + 1, callback);
	    });
	}
    };

    var loadcss = function(url) {
	 var resource = document.createElement('link'); 
	resource.setAttribute("rel", "stylesheet");
	resource.setAttribute("href",url);
	resource.setAttribute("type","text/css");      
	var head = document.getElementsByTagName('head')[0];
	head.appendChild(resource);
    };
    cssloads.map(loadcss);

    var plugindraw = function() {
        var obj = mpld3.get_element(iap.props.id);
        var images = iap.props.images;
        var canvas = iap.fig.canvas;
        var settings = iap.props.settings;
	
        var pdiv = $(iap.fig.root[0][0]).parent();
	var plotbox = $('<div>', {"data-columns": "", id: "grid"});
        pdiv.append(plotbox);

	var rowdivs = [];
	for (kindex in Object.keys(images)) {
	    if (kindex % settings["ncols"] == 0) {
		var rowdiv = $('<div>', {"class": "row"});
		plotbox.append(rowdiv);
		rowdivs.push(rowdiv);
	    } else {
		rowdiv = rowdivs[Math.floor(kindex/settings["ncols"])];
	    }
	    // var cdiv = $('<div>');
	    // rowdiv.append(cdiv);
	};
	
        var imgids = {}
        for (imgtype in images) {
            imgids[imgtype] = "dynimg-" + imgtype;
        };
	
        for (imgtype in images) {
	    var img = $('<img>', {id:imgids[imgtype], width: "100%"});
	    var span = $('<span>', {id:imgids[imgtype] + '-title'});
	    span.text(imgtype);
	    
	    var content = MATDB_PANEL.replace(/\{\{(\w+)\}\}/g, function (match, g1) {
                switch (g1) {
                case 'imgtype':
                    return span[0].outerHTML;
                    break;
                case 'content':
                    return img[0].outerHTML;
                    break;
                }
            });
                                        
            var item = document.createElement('div');
            salvattore['append_elements'](plotbox[0], [item]);
            item.outerHTML = content;
        };

        obj.elements().on("mousedown", function(d, i) {
            for (imgtype in images) {
                var urls = images[imgtype];
                $('#' + imgids[imgtype]).attr("src", urls[i]);
		var title = imgtype + ': ' + settings["names"][i];
		$('#' + imgids[imgtype] + '-title').text(title);
            }
        });
    };
    
    recloadjs(0, plugindraw);
};
