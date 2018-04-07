# Off Target Filtering Tools for CRISPR

Computational tools for identifying CRISPR targets and offtargets.

# Preparing an off-target index for the human genome

First, download hg38.fa.gz, the human genome.

The next few steps take about 2 minutes and produce ~4 GB txt file
contaning all CRISPR 20-mer guides in the human genome,
converted to forward direction.

    cd crispr_sites

    make

    gzip -dc generated_files/untracked/hg38.fa.gz | ./crispr_sites > human_targets.txt

The next step loads the off-target guides into memory to facilitate
offtarget filtering.

    cd offtarget
    go build
    HOST=file://`pwd`/../crispr_sites/human_targets.txt ./offtarget

You may have to wait about a minute for the offtarget command to
print that it is ready to receive connections.


# Filtering a batch of targets against the index

Store a list of targets you wish to filter in 

    batch_filter/all_targets.txt

Open another terminal window and in it run

    cd batch_filter
    ./batch_filter.py
    
The targets from your list that match against the human targets
are listed in 

    batch_filter/off_targets.txt

# crispr_sites unit tests

    cd crispr_sites

    make tests
