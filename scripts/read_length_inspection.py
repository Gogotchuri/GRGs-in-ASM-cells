"""
Helper script to inspect FASTQ headers. And get read lengths.
We need this to set the most accurate START sjdb_overhang argument
**It is 63**
"""
import gzip
from pathlib import Path

def inspect_fastq_header(gz_path: str, n: int = 20):
	"""Print first n FASTQ headers to check read length."""
	avg = 0
	c = 0
	with gzip.open(gz_path, 'rt') as f:
		for i, line in enumerate(f):
			if i % 4 == 1:  # Sequence line
				avg += line.strip().__len__()
				c += 1
			if c >= n:
				break
	return avg / c

# Check lengths (first 5) for every FASTQ file
fastq_files = Path("data/raw").glob("*.fastq.gz")
avg = 0
c = 0
for fastq_path in fastq_files:
	avg += inspect_fastq_header(fastq_path.absolute().__str__())
	c += 1
print(f"Average read length: {avg / c:.0f}")
