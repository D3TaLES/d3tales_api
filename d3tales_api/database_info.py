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


def source_groups_generator(group_file=None):
    """
    Get information about database connections. This function requires either a group_file argument
    or a GROUP_FILE environment variable. This argument or  environment variable should be a path
    to a JSON file containing source group information. The keys should be source group names
    and the values should be strings with two digit numbers, e.g., {"Risko": "06"}.

    :param group_file: str, path to database info JSON file

    :return: JSON object containing information about source group codes
    """
    group_file = group_file or os.environ.get('GROUP_FILE') or os.getenv('GROUP_FILE')
    if not group_file:
        print("Environment variable GROUP_FILE not defined. Default group information is in ues. ")
        return {
            "": '00',
            "unknown": '00',
            "Eaton": '01',
            "Robotics": '11',
            "Ganapathysubramanian": '02',
            "Jenkins": '03',
            "Mason": '04',
            "Odom": '05',
            "Odom_Hussein": '05',
            "Odom_Aman": '05',
            "Risko": '06',
            "Sarkar": '07',
            "Shaw": '08',
            "Teague": '09',
            "Csd": '80',
            "Zinc": '81',
            "Risko_Benchmark": '90',
            "Risko_Diheds": '91',
            "Risko_Aman": '92',
            "Risko_Bayesian": '93',
        }
    with open(group_file, "r") as f:
        return json.load(f)


db_info = db_info_generator()
source_groups = source_groups_generator()
