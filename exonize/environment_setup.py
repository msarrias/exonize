import logging
import os
import shutil
from datetime import datetime


class EnvironmentSetup(object):
    logger: logging.Logger

    def __init__(
            self,
            hard_force: bool,
            soft_force: bool,
            draw_event_multigraphs: bool,
            working_directory: str,
            results_database_path: str,
    ):
        self.__FILE_ONLY_INFO = 9
        self.HARD_FORCE = hard_force
        self.SOFT_FORCE = soft_force
        self.draw_event_multigraphs = draw_event_multigraphs
        self.working_directory = working_directory
        self.results_database_path = results_database_path
        self.setup_environment()

    def configure_logger(self):
        """
        configure_logger is a function that configures the logger.
        INFO level is used for the log file and WARNING and ERROR
        level for the console.
        """
        logging.addLevelName(
            self.__FILE_ONLY_INFO,
            levelName='FILE_ONLY_INFO'
        )
        logging.Logger.file_only_info = self.file_only_info
        self.logger = logging.getLogger(__name__)
        self.logger.setLevel(logging.INFO)

        # Define a filter that allows all messages EXCEPT those at INFO level
        class ExcludeInfoFilter(logging.Filter):
            def filter(self, record):
                return record.levelno != logging.INFO

        # Define file handler for the "FILE_ONLY_INFO" level
        log_file_name = f"exonize_log_{datetime.now():%Y%m%d_%H%M%S}.log"
        log_file_path = os.path.join(self.working_directory, log_file_name)
        file_handler = logging.FileHandler(log_file_path, mode='w')
        file_handler.setLevel(self.__FILE_ONLY_INFO)  # Changed level to custom level
        file_handler.addFilter(ExcludeInfoFilter())  # Added filter
        file_handler.setFormatter(logging.Formatter('%(message)s'))

        # Define console handler for the "INFO" level and above ("DEBUG", "WARNING", "EXCEPTION").
        console_handler = logging.StreamHandler()
        console_handler.setLevel(logging.INFO)
        console_handler.setFormatter(logging.Formatter('[%(levelname)s]: %(message)s'))

        # Add handlers to the logger
        self.logger.addHandler(file_handler)
        self.logger.addHandler(console_handler)

    # noinspection PyProtectedMember
    def file_only_info(
            self,
            message,
            *args,
            **kws
    ):
        if self.logger.isEnabledFor(self.__FILE_ONLY_INFO):
            self.logger._log(self.__FILE_ONLY_INFO, message, args, **kws)

    def setup_environment(self):
        if self.HARD_FORCE:
            if os.path.exists(self.working_directory):
                shutil.rmtree(self.working_directory)
        elif self.SOFT_FORCE:
            if os.path.exists(self.results_database_path):
                os.remove(self.results_database_path)
        os.makedirs(self.working_directory, exist_ok=True)
        if self.draw_event_multigraphs:
            os.makedirs(os.path.join(self.working_directory, 'multigraphs'), exist_ok=True)
        self.configure_logger()
