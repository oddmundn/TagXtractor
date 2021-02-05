# TagXtractor
Software for extraction of molecular tag from sequencing data

Isolate molecular tags, move them from within the sequenced read to the
header region, and remove the spacer region. The ion version is adapted for
Ion Torrent chemistry and a simplified iDES system in Stavanger.
This version writes the resulting fastq file to stdout for piping
It can also read from stdin if using the infile file name -
Current version also has an option to trim the ends of the reads
