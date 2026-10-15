use std::{
    fs::{File, OpenOptions},
    io::{self, Read, Write},
    path::Path,
    sync::{
        Arc,
        atomic::{AtomicU64, Ordering},
    },
};

#[derive(Clone)]
pub struct ConcurrentFileWriter {
    file_and_offset: Arc<(File, AtomicU64)>,
}

impl ConcurrentFileWriter {
    pub fn create(path: impl AsRef<Path>) -> io::Result<Self> {
        let mut file = OpenOptions::new()
            .create(true)
            .write(true)
            .read(true)
            .truncate(true)
            .open(path)?;

        let empty = zstd::bulk::compress(&[], 1)?;
        file.write_all(&empty)?;
        Ok(Self {
            file_and_offset: Arc::new((file, AtomicU64::new(empty.len() as u64))),
        })
    }

    pub fn append_existing(path: impl AsRef<Path>) -> io::Result<Self> {
        let file = OpenOptions::new()
            .create(true)
            .write(true)
            .read(true)
            .open(path)?;

        let len = file.metadata()?.len();

        Ok(Self {
            file_and_offset: Arc::new((file, AtomicU64::new(len))),
        })
    }

    /// Atomically reserves space, writes `buf`, and returns the starting offset.
    pub fn write(&self, buf: &[u8]) -> io::Result<u64> {
        let compressed = zstd::bulk::compress(buf, 1)?;
        let len = u64::try_from(compressed.len())
            .map_err(|_| io::Error::new(io::ErrorKind::InvalidInput, "buffer too large"))?;

        let offset = self.file_and_offset.1.fetch_add(len, Ordering::Relaxed);

        write_all_at(&self.file_and_offset.0, &compressed, offset)?;

        Ok(offset)
    }

    /// Resizes `buf` to `len` and reads exactly `len` bytes from `offset`.
    ///
    /// Existing contents are overwritten from the start of the vector.
    pub fn read_all_at(&self, buf: &mut Vec<u8>, offset: u64, len: usize) -> io::Result<()> {
        struct AtReader<'a> {
            file: &'a File,
            offset: u64,
        }
        impl Read for AtReader<'_> {
            fn read(&mut self, out: &mut [u8]) -> io::Result<usize> {
                let n = read_at(self.file, out, self.offset)?;
                self.offset += n as u64;
                Ok(n)
            }
        }
        buf.clear();
        let reader = AtReader {
            file: &self.file_and_offset.0,
            offset,
        };
        zstd::stream::read::Decoder::new(reader)?
            .single_frame()
            .take(len as u64)
            .read_to_end(buf)?;
        if buf.len() != len {
            return Err(io::Error::new(
                io::ErrorKind::InvalidData,
                "indirect zstd record length mismatch",
            ));
        }
        Ok(())
    }

    pub fn file(&self) -> &File {
        &self.file_and_offset.0
    }
}

#[cfg(unix)]
fn write_all_at(file: &File, mut buf: &[u8], mut offset: u64) -> io::Result<()> {
    use std::os::unix::fs::FileExt;

    while !buf.is_empty() {
        let written = file.write_at(buf, offset)?;

        if written == 0 {
            return Err(io::Error::new(
                io::ErrorKind::WriteZero,
                "failed to write whole buffer",
            ));
        }

        buf = &buf[written..];
        offset += written as u64;
    }

    Ok(())
}

#[cfg(windows)]
fn write_all_at(file: &File, mut buf: &[u8], mut offset: u64) -> io::Result<()> {
    use std::os::windows::fs::FileExt;

    while !buf.is_empty() {
        let written = file.seek_write(buf, offset)?;

        if written == 0 {
            return Err(io::Error::new(
                io::ErrorKind::WriteZero,
                "failed to write whole buffer",
            ));
        }

        buf = &buf[written..];
        offset += written as u64;
    }

    Ok(())
}

#[cfg(unix)]
fn read_at(file: &File, buf: &mut [u8], offset: u64) -> io::Result<usize> {
    use std::os::unix::fs::FileExt;
    file.read_at(buf, offset)
}

#[cfg(windows)]
fn read_at(file: &File, buf: &mut [u8], offset: u64) -> io::Result<usize> {
    use std::os::windows::fs::FileExt;
    file.seek_read(buf, offset)
}
