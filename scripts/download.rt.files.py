import boto3
from botocore.config import Config
from botocore import UNSIGNED
import os
from botocore.exceptions import ClientError

## list subfolders
s3 = boto3.client('s3', config=Config(signature_version=UNSIGNED))

bucket_name = 'gencove-sbir'
subfolder = 'imputation-reference/'
path = f"{bucket_name}/"


if subfolder == '/':
    response = s3.list_objects_v2(Bucket=bucket_name, Delimiter='/')
else:
    response = response = s3.list_objects_v2(Bucket=bucket_name, Prefix=subfolder, Delimiter='/')
    path = f"{path}/{subfolder}"
if 'CommonPrefixes' in response:
    print(f"Folders in {path}:")
    
    for folder in response['CommonPrefixes']:
        print(folder.get('Prefix'))
else:

    print("No folders found.")

## lists files
# Create an S3 client

s3 = boto3.client('s3', config=Config(signature_version=UNSIGNED))

bucket_name = 'gencove-sbir'
subfolder = 'imputation-reference/g1k_v37/chr6/'
response = response = s3.list_objects_v2(Bucket=bucket_name, Prefix=subfolder)

if 'Contents' in response:
    print(f"Files in the {subfolder}:")
    for obj in response['Contents']:
        key = obj['Key']
        # Optionally, skip if the key exactly matches the subfolder (if a folder placeholder exists)
        if key == subfolder:
            continue
        print(key)
else:

    print("No files found in the subfolder.")


## download subfolder
def download_folder(bucket_name, s3_folder, local_dir):
    s3 = boto3.client('s3', config=Config(signature_version=UNSIGNED))
    paginator = s3.get_paginator('list_objects_v2')
    for page in paginator.paginate(Bucket=bucket_name, Prefix=s3_folder):
        for obj in page.get('Contents', []):
            key = obj['Key']
            if not key.endswith('/'):
                local_file_path = os.path.join(local_dir, key[len(s3_folder):].lstrip('/'))
                os.makedirs(os.path.dirname(local_file_path), exist_ok=True)
                try:
                    s3.download_file(bucket_name, key, local_file_path)
                    print(f"Downloaded {key} to {local_file_path}")
                except ClientError as e:
                    print(f"Error downloading {key}: {e}")

 

# Example usage
bucket_name = 'gencove-sbir'
s3_folder = 'imputation-reference/g1k_v37/chunks'
local_dir = 'chunks/'

 

download_folder(bucket_name, s3_folder, local_dir)