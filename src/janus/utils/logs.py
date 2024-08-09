import logging
import sys

# Simple terminal-logger instance 
def SetupLogger(level:str="INFO", name:str="FWL"):

    custom_logger = logging.getLogger(name)
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

def GetLogger(name="FWL"):
    # Get logger by name 
    log = logging.getLogger(name)

    # Setup logger if not already 
    if not log.hasHandlers():
        SetupLogger(name)
    
    # Return the logger
    return log
