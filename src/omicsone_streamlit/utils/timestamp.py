import os
from datetime import datetime

def build_timestamp():
# Get the current timestamp
    timestamp = datetime.now().strftime("%Y%m%d_%H%M%S")
    return timestamp


