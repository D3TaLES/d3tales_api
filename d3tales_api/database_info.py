import os
import json
import warnings


def db_info_generator(db_file=None):
    """
    Get information about database connections. This function requires either a db_file argument
    or a DB_INFO_FILE environment variable. This argument or  environment variable should be a path
    to a JSON file containing connection information for the databases. The keys should be the database
    names such as `frontend`, `backend`, `expflow`, and `fireworks`.

    :param db_file: str, path to database info JSON file

    :return: JSON object containing information about database connections
    """
    db_file = db_file or os.environ.get('DB_INFO_FILE') or os.getenv('DB_INFO_FILE')
    if not db_file:
        warnings.warn("Environment variable DB_INFO_FILE not defined. Default database information is in ues. ")
        return {
            "frontend":
                {
                    "host": "mongodb://USERNAME:PASSWORD@DATABASE_IP:DATABASE_PORT/frontend",
                    "database": "ui"
                },
            "backend":
                {
                    "host": "mongodb://USERNAME:PASSWORD@DATABASE_IP:DATABASE_PORT/backend",
                    "database": "backend"
                }
        }
    with open(db_file, "r") as f:
        return json.load(f)


db_info = db_info_generator()
