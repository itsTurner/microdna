# MicroDNA detection
Detect and quantify microDNA occurences in alignment data.

# Command Line Interface
**Dependencies listed in `requirements.txt`. To install, run:**
```bash
$ pip install -r requirements.txt
```

## Usage
```
usage: src/explore.py [-h] -r REFERENCE [-t THRESHOLD] [-o OUT_FILE]

options:
  -h, --help            show this help message and exit
  -r REFERENCE, --reference REFERENCE
                        Path to BAM file
  -t THRESHOLD, --threshold THRESHOLD
                        Score threshold to include circle
  -o OUT_FILE, --out_file OUT_FILE
                        Path to log file
```

## Example
```shell
$ python3 src/explore.py \
    -r data/SRR413984.sorted.NC_000001.10.bam \
    -t 1.0 \
    -o docs/results/circles.txt
```
[Results from above run](docs/results/circles.txt)