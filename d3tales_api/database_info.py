import os
import json

db_file = os.environ.get('DB_INFO_FILE') or os.getenv('DB_INFO_FILE')
if not db_file:
    raise TypeError("Environment variable DB_INFO_FILE not defined.")
with open(db_file, "r") as f:
    db_info = json.load(f)
