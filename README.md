# samstrip
This command-line tool removes data from [SAM files](https://en.wikipedia.org/wiki/SAM_(file_format))
which does not contain any information about alignment location or quality.

By default, the SEQ and QUAL fields are emptied, and all auxiliary fields except the NM field is removed.
This can often reduce the size of BAM files by 90%.

## Usage
`samstrip` operates on SAM files only. To work on BAM files, use `samtools` to convert from BAM to SAM
(and optionally back again).
`samstrip` reads the SAM file from stdin and prints the stripped file to stdout.

Example usage:
```
$ # Stripping an exiting BAM file:
$ samtools view -h file.bam | samstrip | samtools view -b - > stripped.bam

$ # Stripping a SAM file while creating it:
$ minimap2 -ax sr ref.fa fw.fq rv.fq | samstrip | samtools view -b - > file.bam
```

### Options
* `--noheader`: Allow input to not have any headers. By default, an input without SAM headers will cause `samstrip` to error.
  The default behaviour is because it's easy to forget passing `-h` to `samtools view`.
  Headers are required for conversion to BAM files, and you should nearly always store your data in BAM format.
* `--keep`: Pass a list of SAM auxiliary tags to keep in output. All others are removed. Defaults to "NM".
Examples:
```
$ cat in | samstrip > out # default: keep tag 'NM' only
$ cat in | samstrip --keep NM AS rl > out # keep tags 'NM', 'AS', 'rl'
$ cat in | samstrip --keep > out # do not keep any tags
```
* `--remove`: Same as `--keep`, but only remove the given tags. If an empty list of tags is passed,
  don't remove any tags. Incompatible with `--keep`.