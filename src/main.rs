use clap::Parser;
use core::str;
use std::io::{self, BufRead, BufReader, BufWriter, Write};

// Iterates over tab-separated fields in u8 slice. This uses memchr for efficiency.
struct FieldIterator<'a> {
    // None if the iterator has no more elements to output,
    // because the previous iteration emitted the remainder of the slice.
    // This allows us to distinguish between 'the iterator is done' and
    // 'there is one more element, it just happens to be empty'
    bytes: Option<&'a [u8]>,
}

impl<'a> FieldIterator<'a> {
    fn new(s: &'a [u8]) -> Self {
        Self { bytes: Some(s) }
    }
}

impl<'a> Iterator for FieldIterator<'a> {
    type Item = &'a [u8];

    fn next(&mut self) -> Option<Self::Item> {
        if let Some(mem) = self.bytes {
            if let Some(next_pos) = memchr::memchr(b'\t', mem) {
                // Safety: Since the memchr found a tab at next_pos, next_pos is a valid index
                // into mem. Therefore, by the behaviour of `split_at_unchecked`,
                // new can't be empty, and therefore this cannot be out of bounds
                unsafe {
                    let (res, new) = mem.split_at_unchecked(next_pos);
                    self.bytes = Some(new.get_unchecked(1..new.len()));
                    Some(res)
                }
            } else {
                self.bytes.take()
            }
        } else {
            None
        }
    }
}

/// Strip a SAM file of non-alignment information. See --help for more info
#[derive(Parser)]
#[command(author, version, about, long_about = HELP_MESSAGE)]
struct Cli {
    /// List of aux tags to keep in file (default: NM)
    #[arg(long, num_args(0..))]
    keep: Option<Vec<String>>,

    /// List of aux tags to remove (incompatible with --keep)
    #[arg(long, num_args(0..))]
    remove: Option<Vec<String>>,

    /// Allow input without SAM header
    #[arg(long, default_value_t = false)]
    noheader: bool,
}

enum Tags {
    KeepAll,
    RemoveAll,
    Keep(Vec<u16>),
    Remove(Vec<u16>),
}

impl Tags {
    fn new(keep: Option<Vec<String>>, remove: Option<Vec<String>>) -> Self {
        match (keep, remove) {
            (Some(_), Some(_)) => {
                eprintln!("Error: Cannot pass both arguments --keep and --remove.");
                std::process::exit(1);
            }
            (None, None) => Self::Keep(convert_tags(vec!["NM".to_string()])),
            (Some(k), None) => {
                if k.is_empty() {
                    Self::RemoveAll
                } else {
                    Self::Keep(convert_tags(k))
                }
            }
            (None, Some(r)) => {
                if r.is_empty() {
                    Self::KeepAll
                } else {
                    Self::Remove(convert_tags(r))
                }
            }
        }
    }

    fn add_fields<W: Write>(&self, io: &mut W, line: &[u8]) {
        let mut need_newline = true;
        match self {
            Self::RemoveAll => (),
            Self::KeepAll => {
                unwrap_pipe(io.write_all(b"\t"));
                unwrap_pipe(io.write_all(line));
                need_newline = false;
            }
            Self::Keep(k) => {
                for mem in FieldIterator::new(line).filter(|m| starts_with_tag(m, k)) {
                    unwrap_pipe(io.write_all(b"\t"));
                    unwrap_pipe(io.write_all(mem));
                    need_newline = mem.last().map_or(false, |&b| b != b'\n');
                }
            }
            Self::Remove(r) => {
                for mem in FieldIterator::new(line).filter(|m| !starts_with_tag(m, r)) {
                    unwrap_pipe(io.write_all(b"\t"));
                    unwrap_pipe(io.write_all(mem));
                    need_newline = mem.last().map_or(false, |&b| b != b'\n');
                }
            }
        }
        if need_newline {
            unwrap_pipe(io.write_all(b"\n"));
        }
    }
}

fn main() {
    let cli = Cli::parse();
    let tags = Tags::new(cli.keep, cli.remove);
    // Buffering the writer leads to a large speedup. Buffering the reader also gives
    // a small speedup
    let mut stdin = BufReader::new(io::stdin().lock());
    let mut stdout = BufWriter::new(io::stdout().lock());
    let mut line = write_header(&mut stdout, &mut stdin, cli.noheader);
    while !line.is_empty() {
        // Write the line from the start to the 9th tab, as the first 9 fields
        // should be copied over unchanged
        let ninth = memchr::memchr_iter(b'\t', &line).nth(8).unwrap_or_else(|| {
            eprintln!("Error: In SAM alignment line, did not see all required fields");
            std::process::exit(1)
        });
        unwrap_pipe(stdout.write_all(&line[..=ninth]));
        // The 10th and 11th field are SEQ and QUAL - we replace them with *
        unwrap_pipe(stdout.write_all(b"*\t*"));
        if let Some(aux_mem) = memchr::memchr_iter(b'\t', &line[ninth + 1..])
            .nth(1)
            .map(|eleventh| &line[eleventh + ninth + 2..])
        {
            tags.add_fields(&mut stdout, aux_mem);
        } else {
            unwrap_pipe(stdout.write_all(b"\n"));
        }
        line.clear();
        stdin.read_until(b'\n', &mut line).unwrap();
    }
}

// Returns the first non-header line (or an empty Vec, if no such line)
fn write_header<O: Write, I: BufRead>(out: &mut O, inp: &mut I, noheader: bool) -> Vec<u8> {
    let first_byte = inp.fill_buf().unwrap().first().copied();
    // Unless --noheader is passed, error if the first byte exists but is not a @
    if !noheader && first_byte.is_some_and(|b| b != b'@') {
        eprintln!(
        "Error: First SAM line did not start with a @, indicating a missing header. \
        \nDid you remember to pass in the whole SAM file? \
        If producing the SAM file using `samtools view`, remember the `-h` flag to include the header."
    );
        std::process::exit(1);
    }
    let mut line: Vec<u8> = Vec::new();
    let mut seen_versions: Vec<u32> = Vec::new();
    while inp.read_until(b'\n', &mut line).unwrap() > 0 {
        if !line.starts_with(b"@") {
            break;
        };
        // For PG headers, we need a unique ID field. We write them in the format "samstrip.X",
        // where X is some u32. Since it must be unique, we need to keep track of what versions
        // we've already seen in the file, so we can pick the first version not in the file
        if line.starts_with(b"@PG\t") {
            if let Some(ver) = FieldIterator::new(&line)
                .find(|mem| mem.starts_with(b"ID:samstrip."))
                .and_then(|mem| str::from_utf8(&mem[b"ID:samstrip.".len()..]).ok())
                .and_then(|s| s.parse::<u32>().ok())
            {
                seen_versions.push(ver)
            }
        }
        unwrap_pipe(out.write_all(&line));
        line.clear();
    }
    // The idea here is that after sorting and dedupping, seen_versions must correspond contain the
    // number 0..seen_versions.len(), UNLESS a version is missing.
    // So, the first unseen version is: The first element in seen_versions that is not equal to its
    // index, or if no such element exist, it's the length of `seen_versions` (i.e. the next element
    // that would have been added to the vec)
    seen_versions.sort_unstable();
    // Technically this should already be deduplicated, because samstrip will never write
    // a PG ID that already exists. Nonetheless, we do it anyway for safety's sake.
    seen_versions.dedup();
    let first_unseen_version = seen_versions
        .iter()
        .enumerate()
        .position(|(i, &v)| i as u32 != v)
        .unwrap_or(seen_versions.len()) as u32;
    if first_byte.is_some_and(|b| b == b'@') {
        write_pg_header(out, first_unseen_version);
    }
    line
}

// This is useful for consumers of the stripped SAM file, because they probably wonder
// where all the fields went. The PG header will tell them that this SAM file was processed
// with samstrip, and which version, and what the program "samstrip" even is.
fn write_pg_header<O: Write>(out: &mut O, id: u32) {
    let command_line = std::env::args()
        .map(|arg| arg.as_str().to_string())
        .collect::<Vec<_>>()
        .join(" ");
    let description = "See https://github.com/jakobnissen/samstrip";
    let line = format!(
        "@PG\tID:samstrip.{}\tPN:samstrip\tVN:{}\tCL:{}\tDS:{}\n",
        id,
        env!("CARGO_PKG_VERSION"),
        command_line,
        description,
    );
    unwrap_pipe(out.write_all(line.as_bytes()));
}

// SAM aux tags are two ASCII alphanumerical symbols. We encode them in a u16, so we can check
// for their presence more quickly.
fn convert_tags(tags: Vec<String>) -> Vec<u16> {
    tags.iter().map(|tag| {
        match tag.as_bytes() {
            [x @ (b'A'..=b'Z' | b'a'..=b'z'), y @ (b'A'..=b'Z' | b'a'..=b'z' | b'0'..=b'9')] => {
                (*x as u16) | ((*y as u16) << 8)
            },
            _ => {eprintln!("Error: Tag argument \"{}\" does not conform to SAM specified regex [A-Za-z][A-Za-z0-9]", tag);
            std::process::exit(1)},
        }
    }).collect()
}

fn starts_with_tag(mem: &[u8], tags: &[u16]) -> bool {
    if mem.len() < 2 {
        eprintln!("Error: Could not parse SAM auxiliary field: Less than two bytes in length");
        std::process::exit(1)
    }
    // Safety: We read first two elements and just checked there are at least two elements.
    // Note that this works correctly for both big- and little endian platforms,
    // but it's more efficient for little-endian ones.
    let u = unsafe { (*mem.get_unchecked(0) as u16) | ((*mem.get_unchecked(1) as u16) << 8) };
    // We could use tags.contain(&u), but this implementation SIMDs.
    // That probably won't matter in 99% of cases, but the user could pass a long list of tags,
    // in which case it would.
    // SIMD'ing does not degrade performance measurably, because the extra branches related to whether
    // the SIMD loop should kick in is constant across the program's run, and therefore will be almost
    // perfectly predicted by the branch predictor.
    tags.iter().fold(false, |acc, &t| acc | (t == u))
}

// This is necessary because stdout may decide to stop allowing more input. E.g. if
// the user does `cat file.sam | samstrip | head -10`, then `head` will close the pipe
// after 10 lines, which will make writing to the pipe return a BrokenPipe error.
// This particular failure more is not really a failure, but indicates samstrip has no
// more output to produce, so we just exit the program immediately.
fn unwrap_pipe<T>(x: std::io::Result<T>) -> T {
    match x {
        Ok(r) => r,
        Err(ref err) => {
            if err.kind() == std::io::ErrorKind::BrokenPipe {
                std::process::exit(0)
            } else {
                x.unwrap()
            }
        }
    }
}

const HELP_MESSAGE: &str =
    "Reads a SAM file from stdin, and prints the equivalent stripped file to stdout.
A stripped file has the SEQ and QUAL fields removed, and all auxiliary fields.
Barring any aligner-specific auxiliary fields, a stripped SAM file contain the same
alignment information as a full file, but takes up less disk space.

The program optionally takes `--keep`, a list of auxiliary tags to keep in the output.
This defaults to 'NM', which many tools assume is always present.
Examples:
`cat in | samstrip > out` - default: keep tag 'NM' only
`cat in | samstrip --keep NM AS rl > out` - keep tags 'NM', 'AS', 'rl'
`cat in | samstrip --keep > out` - do not keep any tags

Similarly, the `--remove` option only removes the given tags. If no tags are passed,
all tags are kept.

Example usage:
Stripping an exiting BAM file:
`samtools view -h file.bam | samstrip | samtools view -b - > stripped.bam`

Stripping a BAM file while creating it:
`minimap2 -ax sr ref.fa fw.fq rv.fq | samstrip | samtools view -b - > file.bam`";
