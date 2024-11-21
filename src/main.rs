use clap::Parser;
use core::str;
use std::io::{self, BufRead, BufReader, BufWriter, Write};

struct FieldIterator<'a> {
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
                let (res, new) = mem.split_at(next_pos);
                // Safety: New is guaranteed to be nonempty by the semantics of split_at
                unsafe { self.bytes = Some(new.get_unchecked(1..new.len())) }
                Some(res)
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
    /// List of auxiliary tags to keep in file (default: NM)
    #[arg(long, num_args(0..))]
    keep: Option<Vec<String>>,

    /// Allow input without SAM header
    #[arg(long, default_value_t = false)]
    noheader: bool,
}

fn main() {
    let cli = Cli::parse();
    // Basically, if user didn't pass it in, default to NM, if user explicitly passed
    // an empty list of tags, return None.
    let keep_tags = {
        let strtags = cli.keep.unwrap_or(vec!["NM".to_string()]);
        if strtags.is_empty() {
            None
        } else {
            Some(convert_tags(strtags))
        }
    };
    let mut stdin = BufReader::new(io::stdin().lock());
    let mut stdout = BufWriter::new(io::stdout().lock());
    let mut line = write_header(&mut stdout, &mut stdin, cli.noheader);
    while !line.is_empty() {
        // Write the line from the start to the 9th tab, as the first 9 fields
        // should be copied over unchanged
        let start = memchr::memchr_iter(b'\t', &line).nth(8).unwrap_or_else(|| {
            eprintln!("Error: In SAM alignment line, did not see all required fields");
            std::process::exit(1)
        });
        unwrap_pipe(stdout.write_all(&line[..=start]));
        // The 10th and 11th field are SEQ and QUAL - we replace them with *
        unwrap_pipe(stdout.write_all(b"*\t*"));
        let mut needs_newline = true;
        if let Some(ref tags) = keep_tags {
            // Loop over the last fields - skip the first two SEQ and QUAL fields which we replaced above,
            // and then look for the NM:i field which should always be present. If we find it, copy it
            // to the out buffer.
            let fields = FieldIterator::new(&line[start + 1..]);
            for field in fields.skip(2).filter(|f| starts_with_tag(f, tags)) {
                unwrap_pipe(stdout.write_all(b"\t"));
                unwrap_pipe(stdout.write_all(field));
                needs_newline = *field.last().unwrap() != b'\n';
            }
        }
        if needs_newline {
            unwrap_pipe(stdout.write_all(b"\n"));
        }
        line.clear();
        stdin.read_until(b'\n', &mut line).unwrap();
    }
}

fn write_header<O: Write, I: BufRead>(out: &mut O, inp: &mut I, noheader: bool) -> Vec<u8> {
    let first_byte = inp.fill_buf().unwrap().first().copied();
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
        if line.starts_with(b"@PG\t") {
            if let Some(ver) = FieldIterator::new(&line)
                .find(|sl| sl.starts_with(b"ID"))
                .and_then(|mem| str::from_utf8(mem).ok())
                .and_then(|s| s.strip_prefix("ID:samstrip."))
                .and_then(|s| s.parse::<u32>().ok())
            {
                seen_versions.push(ver)
            }
        }
        unwrap_pipe(out.write_all(&line));
        line.clear();
    }
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
    // Safety: We read first two elements and just checked there are at least two
    let u = unsafe { (*mem.get_unchecked(0) as u16) | ((*mem.get_unchecked(1) as u16) << 8) };
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

The program optionally takes `-k`, a list of auxiliary tags to keep in the output.
This defaults to 'NM', which many tools assume is always present.
Examples:
`cat in | samstrip > out` - default: keep tag 'NM' only
`cat in | samstrip --keep NM AS rl > out` - keep tags 'NM', 'AS', 'rl'
`cat in | samstrip --keep > out` - do not keep any tags

Example usage:
Stripping an exiting BAM file:
`samtools view -h file.bam | samstrip | samtools view -b - > stripped.bam`

Stripping a BAM file while creating it:
`minimap2 -ax sr ref.fa fw.fq rv.fq | samstrip | samtools view -b - > file.bam`";
