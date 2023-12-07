import json

def load_config(file_path):
    """
    Load configuration from a JSON file.

    Parameters:
    file_path (str): Path to the configuration file.

    Returns:
    dict: Configuration data.
    """
    try:
        with open(file_path, 'r') as file:
            config = json.load(file)
        return config
    except FileNotFoundError:
        print(f"Configuration file not found: {file_path}")
        return None
    except json.JSONDecodeError:
        print(f"Error decoding JSON from the file: {file_path}")
        return None