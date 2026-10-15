//! Compatibility writer API, now backed by seekable zstd frames as well.
use crate::buckets::writers::compressed_binary_writer::{
    CompressedBinaryWriter, CompressedCheckpointSize, CompressionLevelInfo,
};
use crate::buckets::LockFreeBucket;
use crate::memory_data_size::MemoryDataSize;
use crate::memory_fs::file::{internal::MemoryFileMode, reader::FileRangeReference};
use crate::utils::memory_size_to_log2;
use bincode::Encode;
use std::path::{Path, PathBuf};

pub const LOCK_FREE_BUCKET_MAGIC: &[u8; 16] = b"PLAIN_INTR_BKT_M";
#[derive(Clone)]
pub struct LockFreeCheckpointSize(u8);
impl LockFreeCheckpointSize {
    pub const fn new_from_size(size: MemoryDataSize) -> Self {
        Self(memory_size_to_log2(size))
    }
    pub const fn new_from_log2(value: u8) -> Self {
        Self(value)
    }
}
pub struct LockFreeBinaryWriter {
    inner: CompressedBinaryWriter,
}
impl LockFreeBinaryWriter {
    pub const CHECKPOINT_SIZE_UNLIMITED: LockFreeCheckpointSize =
        LockFreeCheckpointSize::new_from_log2(62);
}
impl LockFreeBucket for LockFreeBinaryWriter {
    type InitData = (MemoryFileMode, LockFreeCheckpointSize);
    fn new_serialized_data_format(
        path: &Path,
        (mode, checkpoint): &Self::InitData,
        index: usize,
        format: &[u8],
    ) -> Self {
        Self {
            inner: CompressedBinaryWriter::new_serialized_data_format(
                path,
                &(
                    *mode,
                    CompressedCheckpointSize::new_from_log2(checkpoint.0),
                    CompressionLevelInfo {
                        fast_disk: 1,
                        slow_disk: 3,
                    },
                ),
                index,
                format,
            ),
        }
    }
    fn set_checkpoint_data<T: Encode>(&self, data: Option<&T>, range: Option<FileRangeReference>) {
        self.inner.set_checkpoint_data(data, range);
    }
    fn get_bucket_size(&self) -> u64 {
        self.inner.get_bucket_size()
    }
    fn write_data(&self, bytes: &[u8]) -> u64 {
        self.inner.write_data(bytes)
    }
    fn get_path(&self) -> PathBuf {
        self.inner.get_path()
    }
    fn finalize(self) {
        self.inner.finalize();
    }
}
