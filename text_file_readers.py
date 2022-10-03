

def readInputFile(input_file_path):

    with open(input_file_path) as f:
        lines = f.readlines()

    clean_lines = [line.strip() for line in lines]
    clean_lines = [line for line in clean_lines if "#" not in line]
    return clean_lines
