import logging
import sys

# Simple terminal-logger instance 
def SetupLogger(level:str="INFO"):

    custom_logger = logging.getLogger(__name__)
    custom_logger.handlers.clear()

    level = str(level).strip().upper()
    if level not in ["INFO", "DEBUG", "ERROR", "WARNING"]:
        level = "INFO"
    level_code = logging.getLevelName(level)

    # Add terminal output to logger
    sh = logging.StreamHandler(sys.stdout)
    sh.setLevel(level_code)
    custom_logger.addHandler(sh)
    custom_logger.setLevel(level_code)

    # Capture unhandled exceptions
    # https://stackoverflow.com/a/16993115
    def handle_exception(exc_type, exc_value, exc_traceback):
        if issubclass(exc_type, KeyboardInterrupt):
            custom_logger.error("KeyboardInterrupt")
            sys.__excepthook__(exc_type, exc_value, exc_traceback)
            return
        custom_logger.critical("Uncaught exception", exc_info=(exc_type, exc_value, exc_traceback))
    sys.excepthook = handle_exception
    
    return custom_logger
