//! Seekable scratch storage: independent zstd frames, with a private in-memory
//! seek index. The complete file is a valid concatenated zstd stream.
use std::{
    fs::File,
    io::{self, Read, Seek, SeekFrom, Write},
    path::{Path, PathBuf},
};

const BLOCK_BYTES: usize = 64 * 1024;

#[derive(Clone)]
pub(crate) struct Store {
    path: PathBuf,
    frames: Vec<(u64, usize)>,
}

pub(crate) struct StoreWriter {
    file: File,
    store: Store,
    buffer: Vec<u8>,
    compressor: zstd::bulk::Compressor<'static>,
    compressed: Vec<u8>,
}

impl StoreWriter {
    pub fn new(path: &Path) -> io::Result<Self> {
        Ok(Self {
            file: File::create(path)?,
            store: Store {
                path: path.into(),
                frames: Vec::new(),
            },
            buffer: Vec::with_capacity(BLOCK_BYTES),
            compressor: zstd::bulk::Compressor::new(1)?,
            compressed: Vec::with_capacity(zstd::zstd_safe::compress_bound(BLOCK_BYTES)),
        })
    }

    fn frame(&mut self) -> io::Result<()> {
        self.compressed.clear();
        self.compressor
            .compress_to_buffer(&self.buffer, &mut self.compressed)?;
        let offset = self.file.stream_position()?;
        self.file.write_all(&self.compressed)?;
        self.store.frames.push((offset, self.compressed.len()));
        self.buffer.clear();
        Ok(())
    }

    pub fn finish(mut self) -> io::Result<Store> {
        if !self.buffer.is_empty() || self.store.frames.is_empty() {
            self.frame()?;
        }
        self.file.flush()?;
        Ok(self.store)
    }
}

impl Write for StoreWriter {
    fn write(&mut self, mut bytes: &[u8]) -> io::Result<usize> {
        let len = bytes.len();
        while !bytes.is_empty() {
            let take = (BLOCK_BYTES - self.buffer.len()).min(bytes.len());
            self.buffer.extend_from_slice(&bytes[..take]);
            bytes = &bytes[take..];
            if self.buffer.len() == BLOCK_BYTES {
                self.frame()?;
            }
        }
        Ok(len)
    }

    fn flush(&mut self) -> io::Result<()> {
        self.file.flush()
    }
}

pub(crate) struct StoreReader<'a> {
    file: File,
    store: &'a Store,
    cached: usize,
    decoded: Vec<u8>,
    compressed: Vec<u8>,
    decompressor: zstd::bulk::Decompressor<'static>,
}

impl Store {
    pub fn open(&self) -> io::Result<StoreReader<'_>> {
        Ok(StoreReader {
            file: File::open(&self.path)?,
            store: self,
            cached: usize::MAX,
            decoded: Vec::with_capacity(BLOCK_BYTES),
            compressed: Vec::new(),
            decompressor: zstd::bulk::Decompressor::new()?,
        })
    }
}

impl StoreReader<'_> {
    fn load(&mut self, block: usize) -> io::Result<()> {
        if self.cached == block {
            return Ok(());
        }
        let &(offset, len) = self.store.frames.get(block).ok_or_else(|| {
            io::Error::new(io::ErrorKind::UnexpectedEof, "scratch offset out of bounds")
        })?;
        self.compressed.resize(len, 0);
        self.file.seek(SeekFrom::Start(offset))?;
        self.file.read_exact(&mut self.compressed)?;
        self.decoded.clear();
        self.decompressor
            .decompress_to_buffer(&self.compressed, &mut self.decoded)?;
        self.cached = block;
        Ok(())
    }

    pub fn byte(&mut self, offset: usize) -> io::Result<u8> {
        self.load(offset / BLOCK_BYTES)?;
        self.decoded
            .get(offset % BLOCK_BYTES)
            .copied()
            .ok_or_else(|| {
                io::Error::new(io::ErrorKind::UnexpectedEof, "scratch offset out of bounds")
            })
    }

    pub fn record(&mut self, offset: usize) -> io::Result<[u8; 16]> {
        debug_assert!(offset.is_multiple_of(16));
        self.load(offset / BLOCK_BYTES)?;
        Ok(
            self.decoded[offset % BLOCK_BYTES..offset % BLOCK_BYTES + 16]
                .try_into()
                .unwrap(),
        )
    }
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn random_access_and_standard_zstd_decoding() -> io::Result<()> {
        let temp = tempfile::tempdir()?;
        let path = temp.path().join("scratch.zst");
        let expected: Vec<u8> = (0..BLOCK_BYTES * 3 + 123)
            .map(|i| (i % 251) as u8)
            .collect();
        let mut writer = StoreWriter::new(&path)?;
        writer.write_all(&expected)?;
        let store = writer.finish()?;
        assert_eq!(zstd::stream::decode_all(File::open(path)?)?, expected);
        let mut reader = store.open()?;
        for i in [BLOCK_BYTES * 2, 0, BLOCK_BYTES - 1, expected.len() - 1, 17] {
            assert_eq!(reader.byte(i)?, expected[i]);
        }
        Ok(())
    }
}
