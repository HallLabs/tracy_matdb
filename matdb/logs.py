"""Implements logging for managing statistics and debugging for `matdb`
databases, groups, fitters, analyzers, etc.
"""
import logging as log
import logging.handlers, logging.config
from os import path, mkdir
from matdb.base import debug

_filehandler = log.handlers.RotatingFileHandler
"""logging.handlers.RotatingFileHandler: class, *not* instance to allow for this
module to be named logging and not conflict with built-in logging module.
"""
_tracker_format = '%(asctime)-15s %(message)s'
"""str: format string for logging data and statistics from a tracker/camera
combination.
"""

log.config.dictConfig({
    'version': 1,
    'disable_existing_loggers': False,
})

class Logger(object):
    """Handles the logging of events for a `matdb`.

    Args:
        root (str): path to the root logs directory for this logger.
        identifier (str): name of the log (unique in the application).

    Attributes:
        logger (logging.Logger): logger class from built-in logging module.
    """
    def __init__(self, root, identifier):
        self.logger = log.getLogger(identifier)
        #Set the debugging to appropriate level.
        if debug == True:
            self.logger.setLevel(log.DEBUG)
        elif isinstance(debug, int):
            self.logger.setLevel(debug)
        else:
            self.logger.setLevel(log.INFO)
            
        self.root = path.join(root, "logs")
        if not path.isdir(self.root):
            mkdir(self.root)

        logdict = {
            "maxBytes": 10485760,
            "backupCount": 5
        }
        
        self._debuglog = path.join(self.root, "{}.debug.log".format(uuid))
        """str: path to the log file for debug-level messages.
        """
        self.debughandler = _filehandler(self._debuglog, **logdict)
        self.debughandler.setLevel(log.DEBUG)
        self.logger.addHandler(self.debughandler)
        
        self._log = path.join(self.root, "{}.log".format(uuid))
        """str: path to the file for info-level messages.
        """
        self.infohandler = _filehandler(self._log, **logdict)
        self.infohandler.setLevel(log.INFO)
        self.logger.addHandler(self.infohandler)

        self._errorlog = path.join(self.root, "{}.error.log".format(uuid))
        """str: path to the file for error-level messages.
        """
        self.errorhandler = _filehandler(self._errorlog, **logdict)
        self.errorhandler.setLevel(log.WARNING)
        self.logger.addHandler(self.errorhandler)
        
        formatter = log.Formatter(_tracker_format)
        self.debughandler.setFormatter(formatter)
        self.infohandler.setFormatter(formatter)
        self.errorhandler.setFormatter(formatter)

    def exception(self, location, *args):
        """Logs an exception that occurred at `location`.
        """
        self.logger.debug(location, *args, exc_info=True)
        
    def info(self, message, *args):
        """Logs the specified info message for this tracker.
        Args:
            message (str): formatted message; `args` will be placed in-line.
        """
        self.logger.info(message, *args)
        self.infohandler.flush()

    def debug(self, message, *args):
        """Logs the specified debug message for this tracker.
        Args:
            message (str): formatted message; `args` will be placed in-line.
        """
        self.logger.debug(message, *args)

    def error(self, message, *args):
        """Logs the specified error message for this tracker.
        Args:
            message (str): formatted message; `args` will be placed in-line.
        """
        self.logger.error(message, *args)
        self.errorhandler.flush()

    def warning(self, message, *args):
        """Logs the specified warning message for this tracker.
        Args:
            message (str): formatted message; `args` will be placed in-line.
        """
        self.logger.warning(message, *args)
        self.errorhandler.flush()
