import logging
import sys

class LoggerSetup:
    """
    A class that duplicates sys.stdout to a logger.
    Any 'print()' call will be sent to both the console (terminal) 
    and a log file.
    """
    def __init__(self, log_file='simulation.log', level=logging.INFO):
        self.terminal = sys.stdout
        self.logger = logging.getLogger('ReactorSimLogger')
        self.logger.setLevel(level)

        # Remove existing handlers to avoid duplicates if run in the same session
        if self.logger.hasHandlers():
            self.logger.handlers.clear()

        # Create file handler
        # 'w' mode overwrites the file each time. Use 'a' to append.
        file_handler = logging.FileHandler(log_file, mode='w')
        file_handler.setLevel(level)
        
        # Set a simple formatter to just log the message itself
        file_handler.setFormatter(logging.Formatter('%(message)s'))
        
        self.logger.addHandler(file_handler)

    def write(self, message):
        """
        This method is called by 'print()'.
        """
        # 1. Write to the original stdout (the console)
        self.terminal.write(message)
        
        # 2. Write to the log file, stripping extra newlines
        if message.strip(): # Avoid logging empty newlines
            self.logger.info(message.rstrip())

    def flush(self):
        """
        This flush method is needed for compatibility with sys.stdout.
        """
        self.terminal.flush()

    def isatty(self):
        """
        Allows the system to check if it's an interactive terminal.
        """
        return self.terminal.isatty()