import configparser
import os


SETTINGS_FILE = "settings.ini"
# Get the absolute path of the current script
current_file_path = os.path.abspath(__file__)

# Get the directory of the current script
current_directory = os.path.dirname(current_file_path)


def read_ini():
    """
    Read and display the contents of a .ini file.
    """
    ini_file = os.path.join(current_directory,SETTINGS_FILE)
    config = configparser.ConfigParser()
    config.read(ini_file)
    ini_dict = {section: dict(config[section]) for section in config.sections()}
    return ini_dict

def update_ini(section, key, value):
    """
    Update or add a key-value pair in the specified section of the .ini file.
    """
    ini_file =  os.path.join(current_directory,SETTINGS_FILE)
    config = configparser.ConfigParser()
    config.read(ini_file)
    
    if section not in config:
        config[section] = {}
    config[section][key] = value

    with open(ini_file, 'w') as configfile:
        config.write(configfile)
    print(f"Updated {ini_file}: [{section}] {key} = {value}")


# Load settings on startup
settings = read_ini()
# print(settings)
DATA_DIR = settings["dirs"]["data_dir"]
OUT_DIR = settings["dirs"]["out_dir"]

FASTA_PATH = settings["paths"]["fasta_path"]
CHROM_PATH = settings["paths"]["chrom_path"]
CYTOBAND_PATH = settings["paths"]["cytoband_path"]

USE_CUSTOMIZED_GENE_MAPPING = settings.get("misc", {}).get("use_customized_gene_mapping", "True").lower() in ['true', '1', 'yes']





