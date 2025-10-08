#!/usr/bin/env python3
"""
Download Pre-annotated USPTO Route Data from Figshare

This script downloads annotated synthesis route files from Figshare and
extracts them to the data/uspto_st/ directory.

Usage:
    python scripts/download_uspto_data.py [--output-dir OUTPUT_DIR]
"""

import argparse
import hashlib
import json
import os
import sys
import gzip
import zipfile
import tarfile
from pathlib import Path
from typing import Dict, List, Optional
from urllib.request import urlretrieve
from tqdm import tqdm


class DownloadProgressBar(tqdm):
    """Progress bar for file downloads."""
    
    def update_to(self, b=1, bsize=1, tsize=None):
        """Update progress bar.
        
        Args:
            b: Number of blocks transferred
            bsize: Size of each block (in bytes)
            tsize: Total size (in bytes)
        """
        if tsize is not None:
            self.total = tsize
        self.update(b * bsize - self.n)


def download_file(url: str, output_path: Path, desc: str = "Downloading") -> None:
    """Download a file with progress bar.
    
    Args:
        url: URL to download from
        output_path: Path to save the file
        desc: Description for progress bar
    """
    output_path.parent.mkdir(parents=True, exist_ok=True)
    
    with DownloadProgressBar(unit='B', unit_scale=True, miniters=1, desc=desc) as t:
        urlretrieve(url, str(output_path), reporthook=t.update_to)


def extract_gz(gz_path: Path, output_path: Path) -> None:
    """Extract a gzip file.
    
    Args:
        gz_path: Path to .gz file
        output_path: Path to save extracted file
    """
    print(f"Extracting {gz_path.name}...")
    with gzip.open(gz_path, 'rb') as f_in:
        with open(output_path, 'wb') as f_out:
            f_out.write(f_in.read())


def extract_zip(zip_path: Path, output_dir: Path) -> None:
    """Extract a zip file.
    
    Args:
        zip_path: Path to .zip file
        output_dir: Directory to extract to
    """
    print(f"Extracting {zip_path.name}...")
    with zipfile.ZipFile(zip_path, 'r') as zip_ref:
        zip_ref.extractall(output_dir)


def extract_tar(tar_path: Path, output_dir: Path) -> None:
    """Extract a tar file.
    
    Args:
        tar_path: Path to .tar or .tar.gz file
        output_dir: Directory to extract to
    """
    print(f"Extracting {tar_path.name}...")
    with tarfile.open(tar_path, 'r:*') as tar_ref:
        tar_ref.extractall(output_dir)


def extract_file(file_path: Path, output_dir: Path) -> None:
    """Automatically detect and extract compressed file.
    
    Args:
        file_path: Path to compressed file
        output_dir: Directory to extract to
    """
    if file_path.suffix == '.gz':
        if file_path.stem.endswith('.tar'):
            # .tar.gz file
            extract_tar(file_path, output_dir)
        else:
            # .json.gz or similar
            output_path = output_dir / file_path.stem
            extract_gz(file_path, output_path)
    elif file_path.suffix == '.zip':
        extract_zip(file_path, output_dir)
    elif file_path.suffix in ['.tar', '.tgz']:
        extract_tar(file_path, output_dir)
    else:
        print(f"Warning: Unknown file type {file_path.suffix}, skipping extraction")


def verify_json_file(file_path: Path) -> bool:
    """Verify that a JSON file is valid.
    
    Args:
        file_path: Path to JSON file
        
    Returns:
        True if valid, False otherwise
    """
    try:
        with open(file_path, 'r') as f:
            json.load(f)
        return True
    except (json.JSONDecodeError, IOError) as e:
        print(f"Error validating {file_path}: {e}")
        return False


def get_figshare_file_list(article_id: str) -> List[Dict]:
    """Fetch file list from Figshare API.
    
    Args:
        article_id: Figshare article ID
        
    Returns:
        List of file information dictionaries
    """
    import urllib.request
    import json
    
    api_url = f"https://api.figshare.com/v2/articles/{article_id}"
    
    try:
        with urllib.request.urlopen(api_url) as response:
            data = json.loads(response.read())
            return data.get('files', [])
    except Exception as e:
        print(f"Error fetching Figshare metadata: {e}")
        return []


def main():
    """Main download script."""
    parser = argparse.ArgumentParser(
        description="Download pre-annotated USPTO route data from Figshare",
        formatter_class=argparse.ArgumentDefaultsHelpFormatter
    )
    parser.add_argument(
        "--output-dir",
        type=Path,
        default=Path("data/uspto_st"),
        help="Directory to save the downloaded data"
    )
    parser.add_argument(
        "--keep-compressed",
        action="store_true",
        help="Keep compressed files after extraction"
    )
    parser.add_argument(
        "--article-id",
        type=str,
        default="30146374",
        help="Figshare article ID"
    )
    
    args = parser.parse_args()
    
    # Create output directory
    args.output_dir.mkdir(parents=True, exist_ok=True)
    temp_dir = args.output_dir / "temp_downloads"
    temp_dir.mkdir(exist_ok=True)
    
    print("=" * 70)
    print("USPTO Annotated Routes Data Downloader")
    print("=" * 70)
    print(f"\nTarget directory: {args.output_dir.resolve()}")
    print(f"Figshare article: https://figshare.com/articles/{args.article_id}")
    print()
    
    # Fetch file list from Figshare
    print("Fetching file list from Figshare API...")
    files = get_figshare_file_list(args.article_id)
    
    if not files:
        print("Error: Could not fetch file list from Figshare.")
        print("\nPlease download manually from:")
        print(f"https://figshare.com/articles/{args.article_id}")
        sys.exit(1)
    
    print(f"Found {len(files)} files to download\n")
    
    # Download and extract each file
    downloaded_files = []
    for i, file_info in enumerate(files, 1):
        file_name = file_info['name']
        file_url = file_info['download_url']
        file_size = file_info.get('size', 0)
        
        print(f"\n[{i}/{len(files)}] Processing: {file_name}")
        print(f"Size: {file_size / (1024**3):.2f} GB")
        
        # Download file
        download_path = temp_dir / file_name
        if download_path.exists():
            print(f"File already downloaded, skipping download...")
        else:
            try:
                download_file(file_url, download_path, desc=f"Downloading {file_name}")
                print(f"✓ Downloaded to {download_path}")
            except Exception as e:
                print(f"✗ Error downloading {file_name}: {e}")
                continue
        
        # Extract file
        try:
            extract_file(download_path, args.output_dir)
            print(f"✓ Extracted to {args.output_dir}")
            downloaded_files.append(file_name)
            
            # Optionally remove compressed file
            if not args.keep_compressed:
                download_path.unlink()
                print(f"✓ Removed compressed file")
        except Exception as e:
            print(f"✗ Error extracting {file_name}: {e}")
            continue
    
    # Clean up temp directory
    try:
        if temp_dir.exists() and not any(temp_dir.iterdir()):
            temp_dir.rmdir()
    except Exception:
        pass
    
    # Verify extracted JSON files
    print("\n" + "=" * 70)
    print("Verifying extracted files...")
    print("=" * 70)
    
    json_files = list(args.output_dir.glob("*.json"))
    valid_files = 0
    
    for json_file in json_files:
        print(f"\nValidating {json_file.name}...", end=" ")
        if verify_json_file(json_file):
            print("✓ Valid")
            valid_files += 1
            
            # Get route count
            try:
                with open(json_file, 'r') as f:
                    data = json.load(f)
                    if isinstance(data, list):
                        print(f"  Contains {len(data)} routes")
            except Exception:
                pass
        else:
            print("✗ Invalid or corrupted")
    
    # Summary
    print("\n" + "=" * 70)
    print("Download Summary")
    print("=" * 70)
    print(f"Downloaded: {len(downloaded_files)} files")
    print(f"Validated: {valid_files} JSON files")
    print(f"Location: {args.output_dir.resolve()}")
    
    if valid_files > 0:
        print("\n✓ Success! You can now use these routes for retrieval and analysis.")
        print("\nExample usage:")
        print(f"  synth-strategy retrieve \\")
        print(f"    --query 'your query' \\")
        print(f"    --route-db {args.output_dir} \\")
        print(f"    --metadata-db data/function_metadata_database.json")
    else:
        print("\n✗ No valid files found. Please check the downloads and try again.")
        sys.exit(1)


if __name__ == "__main__":
    main()
