//! Export records through Helicase for the benchmark's independent validator.
//! Each stdout row is <record name><TAB><DNA>; no FASTA parsing in Python.
use helicase::{input::*, *};
use std::io::{self, BufWriter, Write};

const CONFIG: Config = ParserOptions::default().config();

fn main() -> io::Result<()> {
    let path = std::env::args_os().nth(1).expect("FASTA path required");
    let input = FileInput::new(path)?;
    let mut parser = FastxParser::<CONFIG>::from_input(input)?;
    let mut out = BufWriter::new(io::stdout().lock());
    while parser.next().is_some() {
        let name = parser
            .get_header()
            .split(|b| b.is_ascii_whitespace())
            .next()
            .unwrap_or(&[]);
        out.write_all(name)?;
        out.write_all(b"\t")?;
        out.write_all(parser.get_dna_string())?;
        out.write_all(b"\n")?;
    }
    out.flush()
}
