# %%
import sys

# %%
import yaml
import pysam
import pyBigWig

# %%
import numpy as np
import tensorflow as tf

# %%
import click
from tqdm import tqdm
import pandas as pd
import numpy as np
from scipy.stats import poisson
from pathos.multiprocessing import ProcessingPool as PathosPool

# %%
def nan_to_zero(x):
    """Replaces nan's with zeros."""
    return 0 if math.isnan(x) else x

# %%
class BigWig:
    def __init__(self, bigWigFile) -> None:
        self._bigWig = bigWigFile
    
    def values(self, chrom, start, end):
        # lazy-load bigWig
        if not isinstance(self._bigWig, pyBigWig.pyBigWig):
            self._bigWig = pyBigWig.open(self._bigWig)
        
        try:
            return self._bigWig.values(chrom, start, end)
        except RuntimeError as e:
            print('WARNING: ' + str(e), file=sys.stderr)
            print(f'WARNING: ({chrom}\t{start}\t{end})', file=sys.stderr)
            return [0.0]*(end-start)
    
    
# %%        
class Track:
    def __init__(self, bigWigPlus, bigWigMinus) -> None:
        self._bigWigPlus = BigWig(bigWigPlus)
        self._bigWigMinus = BigWig(bigWigMinus)

    def profile(self, chrom, start, end, strand, reverse=True):
        """Return the signal profile for a given sample. 

        Args:
            chrom  (str): Chromosome (chr1, chr2, ...)
            start  (int): 0-based start position
            end    (int): 0-based end position
            strand (str): Strand ('+' or '-')

        Raises:
            ValueError: [description]

        Returns:
            numpy.array: Numpy array of shape (end-start, 1)
        """

        if strand == '+':
            bigWig = self._bigWigPlus
        elif strand == '-':
            bigWig = self._bigWigMinus
        else:
            raise ValueError(f'Unspected strand: {strand}')

        profile = bigWig.values(str(chrom), start, end)
        profile = [nan_to_zero(c) for c in profile]

        if strand == '-' and reverse:
            profile = list(reversed(profile))

        profile = np.array(profile, dtype=np.float32)

        return profile
    
    def window(self, chrom, center, strand='+', reverse=False, size=201):
        start, end = center - int(size/2), center + int(size/2) + (size % 2)
        return self.profile(chrom, start, end, strand, reverse=reverse)

    
# %%
def poisson_test(observed_counts, mu):
    """Perform a poisson test (H1: Greater) on the observed counts, given an event rate mu."""
    
    return 1 - np.sum(poisson.pmf(range(observed_counts), mu=mu))

# %%
def poisson_count_threshold(mu, pvalue = 0.05):
    """Get the count threshold for a certain p-value and mu."""
    
    p = 1.0
    c = 0
    while p > pvalue:
        p -= poisson.pmf(c, mu=mu)
        c += 1
    return c

# %%
def scan_region(row, track, tile_size = 50, min_pval = 0.05, min_count = 5, min_height=2):
        profile = track.profile(row.chrom, int(row.start), int(row.stop), row.strand, reverse=False)
        if np.sum(profile) <= min_count:
            return
        
        # expected counts in window given the mean counts in the region/transcript
        mu = np.mean(profile) * tile_size
        
        # compute the count threshold for the given p-value and take the most restrictive (i.e. max) threshold
        threshold = max(min_count, poisson_count_threshold(mu, pvalue=min_pval))
        
        # slide over region/profile and find enriched sites
        i = 0
        while i < (len(profile) - tile_size + 1):
            # within-chromosome start/stop positions
            start_i, stop_i = int(row.start)+i, int(row.start)+i+tile_size
            
            # counts in window
            profile_i = profile[i:(i + tile_size)]
            counts_i = int(np.sum(profile_i))
            
            # skip this region if counts are below threshold and minimum peak height is not reached
            if counts_i < threshold or np.max(profile_i) < min_height:
                i += 1
            else:
                # compute p-value and return 'enriched' site
                yield row.chrom, start_i, stop_i, row.name, row.score, row.strand, counts_i, poisson_test(counts_i, mu)
                
                # slide over half the window size (to prevent excessive number of sites in a region)
                i += int(tile_size/2)

# %%
def scan_region_factory(track, tile_size = 50, min_pval = 0.05, min_count = 5, min_height=2):
    """Factory function to wrap scan_region."""
    return lambda row: [site for site in scan_region(row, track, tile_size, min_pval, min_count, min_height)]

# %%
@click.command()
@click.argument('bed')
@click.option('--bigWigPlus')
@click.option('--bigWigMinus')
@click.option('-o', '--output')
@click.option('--tile-size', type=int, default=50)
@click.option('--min-pval', type=float, default=0.01)
@click.option('--min-count', type=int, default=8)
@click.option('--min-height', type=int, default=2)
@click.option('-t', '--threads', type=int, default=1)
def main(bed, bigwigplus, bigwigminus, output, tile_size, min_pval, min_count, min_height, threads):
    track = Track(bigwigplus, bigwigminus)
    bed_df = pd.read_csv(bed, sep='\t', header=None)
    bed_df.columns = ['chrom', 'start', 'stop', 'name', 'score', 'strand']
    
    with PathosPool(threads) as p:
        result = p.uimap(scan_region_factory(track, tile_size, min_pval, min_count, min_height), tqdm(bed_df.itertuples(), total=len(bed_df)))
    
    with open(output, 'w') as f:
        for sites in result:
            for site in sites:
                print('\t'.join(map(str, site)), file=f)

# %%
if __name__ == '__main__':
    main()