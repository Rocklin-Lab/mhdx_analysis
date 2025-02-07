import pandas as pd
import glob
import os
import uuid
import argparse
import shutil

def process_files(input_dir, output_file, temp_dir):
    files = glob.glob(os.path.join(input_dir, "*.json"))
    if not os.path.exists(temp_dir):
        os.makedirs(temp_dir)

    temp_files = []
    for i, file_path in enumerate(files):
        if i % 50 == 0 and i != 0:
            print(f"Processed {i} files. {len(files) - i} files to go...")
            df = pd.concat(temp_files).reset_index(drop=True)
            temp_file_path = os.path.join(temp_dir, f"temp_{i}.json")
            df.to_json(temp_file_path)
            temp_files = [pd.read_json(temp_file_path)]

        df = pd.read_json(file_path)
        df['name'] = os.path.basename(file_path).replace(".json", "")
        temp_files.append(df)

#    if temp_files:  # Final concatenation for remaining files
    df = pd.concat(temp_files).reset_index(drop=True)
    final_temp_file = os.path.join(temp_dir, f"temp_final.json")
    df.to_json(output_file)

    print(f"Finished processing. Output saved to {output_file}")

def main():
    parser = argparse.ArgumentParser(description='Concatenate JSON files from a folder into a single JSON file.')
    parser.add_argument('input_dir', type=str, help='Input directory containing JSON files to concatenate')
    parser.add_argument('output_file', type=str, help='Output JSON file path')
    args = parser.parse_args()

    temp_dir = os.path.join("temp", str(uuid.uuid4()))
    try:
        process_files(args.input_dir, args.output_file, temp_dir)
    finally:
        # Cleanup temporary directory
        if os.path.exists(temp_dir):
            shutil.rmtree(temp_dir)
            print("Temporary files cleaned up.")

if __name__ == "__main__":
    main()
