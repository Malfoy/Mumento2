use hashes::HashableSequence;
use helicase::{input::*, *};
use std::{io::Read, path::Path};

// The bucketing queue consumes ASCII and packs it itself. Reading packed DNA
// here would immediately unpack every base, only to pack it again downstream.
const CONFIG: Config = ParserOptions::default().config();

#[derive(Copy, Clone)]
pub enum DnaSequencesFileType {
    FASTA,
    FASTQ,
    GFA,
    BINARY,
}

#[derive(Copy, Clone)]
pub struct DnaSequence<'a, S: HashableSequence + 'a> {
    pub ident_data: &'a [u8],
    pub seq: S,
    pub format: DnaSequencesFileType,
}

pub struct SequencesReader;
impl SequencesReader {
    pub fn new() -> Self {
        Self
    }

    pub fn process_file_extended<F: for<'a> FnMut(DnaSequence<'a, &'a [u8]>)>(
        &mut self,
        source: impl AsRef<Path>,
        func: F,
        line_split_copyback: Option<usize>,
        copy_ident: bool,
        remove_file: bool,
    ) {
        let data = FileInput::new(source.as_ref()).expect("Cannot open input sequence file");
        if data.first_byte() == 0 {
            if remove_file {
                std::fs::remove_file(source).unwrap();
            }
            return;
        }
        let mut parser = FastxParser::<CONFIG>::from_input(data)
            .unwrap_or_else(|e| panic!("Cannot parse {}: {e}", source.as_ref().display()));
        Self::process_parser(&mut parser, func, line_split_copyback, copy_ident);
        drop(parser);
        if remove_file {
            std::fs::remove_file(source).unwrap();
        }
    }

    pub fn process_reader_extended(
        &mut self,
        reader: &mut dyn Read,
        _name: &Path,
        func: impl for<'a> FnMut(DnaSequence<'a, &'a [u8]>),
        line_split_copyback: Option<usize>,
        copy_ident: bool,
    ) {
        let data = ReaderInput::new(reader);
        if data.first_byte() == 0 {
            return;
        }
        let mut parser =
            FastxParser::<CONFIG>::from_input(data).expect("Cannot parse sequence stream");
        Self::process_parser(&mut parser, func, line_split_copyback, copy_ident);
    }

    fn process_parser<P: HelicaseParser + Iterator>(
        parser: &mut P,
        mut func: impl for<'a> FnMut(DnaSequence<'a, &'a [u8]>),
        copyback: Option<usize>,
        copy_ident: bool,
    ) {
        let mut normalized = Vec::new();
        let mut identifier = Vec::new();
        while parser.next().is_some() {
            identifier.clear();
            let format = match parser.format() {
                helicase::parser::Format::Fasta => DnaSequencesFileType::FASTA,
                helicase::parser::Format::Fastq => DnaSequencesFileType::FASTQ,
            };
            if copy_ident {
                identifier.push(if matches!(format, DnaSequencesFileType::FASTA) {
                    b'>'
                } else {
                    b'@'
                });
                identifier.extend_from_slice(parser.get_header());
            }
            let sequence = parser.get_dna_string();
            if sequence.is_empty() {
                continue;
            }
            let overlap = copyback.unwrap_or(0);
            let size = if copyback.is_some() {
                ::config::DEFAULT_OUTPUT_BUFFER_SIZE
                    .max(overlap.saturating_mul(2))
                    .max(overlap + 1)
            } else {
                sequence.len()
            };
            let mut start: usize = 0;
            loop {
                let end = start.saturating_add(size).min(sequence.len());
                normalized.clear();
                normalized.extend(sequence[start..end].iter().map(|&base| {
                    let upper = base.to_ascii_uppercase();
                    match upper {
                        b'A' | b'C' | b'G' | b'T' => upper,
                        _ => b'N',
                    }
                }));
                func(DnaSequence {
                    ident_data: &identifier,
                    seq: &normalized,
                    format,
                });
                if end == sequence.len() {
                    break;
                }
                start = end - overlap;
            }
        }
    }
    pub fn archive_file_type(path: &Path) -> Option<DnaSequencesFileType> {
        let name = path.file_name()?.to_str()?;
        let name = crate::sequences_stream::tar::strip_compression_suffix(name);
        match Path::new(name).extension()?.to_str()? {
            "fq" | "fastq" => Some(DnaSequencesFileType::FASTQ),
            "fa" | "fasta" | "fna" | "ffn" => Some(DnaSequencesFileType::FASTA),
            _ => None,
        }
    }

    pub fn file_type(path: &Path) -> Option<DnaSequencesFileType> {
        let mut name = path.file_name()?.to_str()?;
        while let Some((stem, extension)) = name.rsplit_once('.') {
            match extension {
                "fq" | "fastq" => return Some(DnaSequencesFileType::FASTQ),
                "fa" | "fasta" | "fna" | "ffn" => return Some(DnaSequencesFileType::FASTA),
                _ => name = stem,
            }
        }
        None
    }
}
